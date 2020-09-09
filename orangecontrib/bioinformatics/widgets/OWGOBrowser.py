""" GO Enrichment Analysis """
import gc
import sys
import math
import logging
import operator
import itertools
import webbrowser
from types import SimpleNamespace
from typing import Dict, List, Tuple
from functools import reduce
from collections import defaultdict
from concurrent.futures import Future

import numpy
from requests.exceptions import ConnectTimeout, ConnectionError, RequestException

from AnyQt.QtGui import QBrush
from AnyQt.QtCore import Qt, QSize, QTimer, Signal, QThread
from AnyQt.QtCore import pyqtSlot as Slot
from AnyQt.QtWidgets import (
    QMenu,
    QLabel,
    QDialog,
    QCheckBox,
    QSplitter,
    QTreeView,
    QTreeWidget,
    QVBoxLayout,
    QItemDelegate,
    QTreeWidgetItem,
)

import Orange.data
from Orange.widgets import gui, widget, settings
from Orange.widgets.utils.concurrent import Task, FutureWatcher, ThreadExecutor, methodinvoke

from orangecontrib.bioinformatics import go
from orangecontrib.bioinformatics.ncbi import gene, taxonomy
from orangecontrib.bioinformatics.utils import statistics, serverfiles
from orangecontrib.bioinformatics.go.config import DOMAIN, FILENAME_ONTOLOGY, FILENAME_ANNOTATION
from orangecontrib.bioinformatics.widgets.utils.data import TableAnnotation, check_table_annotation


class EnsureDownloaded(Task):
    advance = Signal()
    progress = Signal(float)

    def __init__(self, files_list, parent=None):
        Task.__init__(self, parent)
        self.files_list = files_list

    def run(self):
        nfiles = len(self.files_list)
        count = 100 * nfiles
        counter = itertools.count()

        def advance():
            self.advance.emit()
            self.progress.emit(100.0 * next(counter) / count)

        for domain, filename in self.files_list:
            ensure_downloaded(domain, filename, advance=advance)


def ensure_downloaded(domain, filename, advance=None):
    serverfiles.localpath_download(domain, filename, callback=advance)


def isstring(var):
    return isinstance(var, Orange.data.StringVariable)


class TreeNode(object):
    def __init__(self, value, children):
        self.value = value
        self.children = children


class GOTreeWidget(QTreeWidget):
    def __init__(self, parent=None):
        super().__init__(parent=parent)
        self._currMenu = QMenu()
        self._currAction = self._currMenu.addAction("View term on AmiGO website")

    def contextMenuEvent(self, event):
        super().contextMenuEvent(event)

        def browser_action(_term):
            if isinstance(_term, go.Term):
                _term = _term.id
            webbrowser.open(f'http://amigo.geneontology.org/amigo/term/{_term}')

        term = self.itemAt(event.pos()).term
        self._currAction.triggered.connect(lambda: browser_action(term))
        self._currMenu.popup(event.globalPos())


class State:
    #: Ready to run
    Ready = 1
    #: Downloading datasets needed to run.
    Downloading = 2
    #: Running the enrichment task
    Running = 4
    #: The current executing task is stale. Need to reschedule another update
    #: once this one completes.
    #: Only applicable with Downloading and Running
    Stale = 8


class OWGOBrowser(widget.OWWidget):
    name = "GO Browser"
    description = "Enrichment analysis for Gene Ontology terms."
    icon = "../widgets/icons/OWGOBrowser.svg"
    priority = 60

    inputs = [
        ("Cluster Data", Orange.data.Table, "set_dataset", widget.Single + widget.Default),
        ("Reference Data", Orange.data.Table, "set_reference_dataset"),
    ]

    outputs = [("Data on Selected Genes", Orange.data.Table), ("Enrichment Report", Orange.data.Table)]

    settingsHandler = settings.DomainContextHandler()

    gene_attr_index = settings.ContextSetting(0)
    use_attr_names = settings.ContextSetting(False)
    use_reference_dataset = settings.Setting(False)
    aspect_index = settings.Setting(0)
    use_evidence_type = settings.Setting({et: True for et in go.evidence_types_ordered})
    filter_by_num_of_instances = settings.Setting(False)
    min_num_of_instances = settings.Setting(1)
    filter_by_p_value = settings.Setting(True)
    max_p_value = settings.Setting(0.2)
    filter_by_p_value_nofdr = settings.Setting(False)
    max_p_value_no_fdr = settings.Setting(0.01)
    prob_func = settings.Setting(0)
    selection_direct_annotation = settings.Setting(0)
    selection_disjoint = settings.Setting(0)

    class Error(widget.OWWidget.Error):
        serverfiles_unavailable = widget.Msg(
            'Can not locate annotation files, ' 'please check your connection and try again.'
        )

    def __init__(self, parent=None):
        super().__init__(self, parent)

        self.input_data = None
        self.gene_info = None
        self.ref_data = None
        self.ontology = None
        self.annotations = None
        self.loaded_annotation_code = None
        self.treeStructRootKey = None
        self.probFunctions = [statistics.Binomial(), statistics.Hypergeometric()]
        self.selectedTerms = []

        self.selectionChanging = 0
        self.__state = State.Ready
        self.__scheduletimer = QTimer(self, singleShot=True)
        self.__scheduletimer.timeout.connect(self.__update)

        #############
        # GUI
        #############
        self.tabs = gui.tabWidget(self.controlArea)
        # Input tab
        self.inputTab = gui.createTabPage(self.tabs, "Input")
        box = gui.widgetBox(self.inputTab, "Info")
        self.infoLabel = gui.widgetLabel(box, "No data on input\n")

        gui.button(
            box,
            self,
            "Ontology/Annotation Info",
            callback=self.show_info,
            tooltip="Show information on loaded ontology and annotations",
        )

        self.referenceRadioBox = gui.radioButtonsInBox(
            self.inputTab,
            self,
            "use_reference_dataset",
            ["Entire genome", "Reference set (input)"],
            tooltips=["Use entire genome for reference", "Use genes from Referece Examples input signal as reference"],
            box="Reference",
            callback=self.__invalidate,
        )

        self.referenceRadioBox.buttons[1].setDisabled(True)
        gui.radioButtonsInBox(
            self.inputTab,
            self,
            "aspect_index",
            ["Biological process", "Cellular component", "Molecular function"],
            box="Aspect",
            callback=self.__invalidate,
        )

        # Filter tab
        self.filterTab = gui.createTabPage(self.tabs, "Filter")
        box = gui.widgetBox(self.filterTab, "Filter GO Term Nodes")
        gui.checkBox(
            box,
            self,
            "filter_by_num_of_instances",
            "Genes",
            callback=self.filter_and_display_graph,
            tooltip="Filter by number of input genes mapped to a term",
        )
        ibox = gui.indentedBox(box)
        gui.spin(
            ibox,
            self,
            'min_num_of_instances',
            1,
            100,
            step=1,
            label='#:',
            labelWidth=15,
            callback=self.filter_and_display_graph,
            callbackOnReturn=True,
            tooltip="Min. number of input genes mapped to a term",
        )

        gui.checkBox(
            box,
            self,
            "filter_by_p_value_nofdr",
            "p-value",
            callback=self.filter_and_display_graph,
            tooltip="Filter by term p-value",
        )

        gui.doubleSpin(
            gui.indentedBox(box),
            self,
            'max_p_value_no_fdr',
            1e-8,
            1,
            step=1e-8,
            label='p:',
            labelWidth=15,
            callback=self.filter_and_display_graph,
            callbackOnReturn=True,
            tooltip="Max term p-value",
        )

        # use filter_by_p_value for FDR, as it was the default in prior versions
        gui.checkBox(
            box, self, "filter_by_p_value", "FDR", callback=self.filter_and_display_graph, tooltip="Filter by term FDR"
        )
        gui.doubleSpin(
            gui.indentedBox(box),
            self,
            'max_p_value',
            1e-8,
            1,
            step=1e-8,
            label='p:',
            labelWidth=15,
            callback=self.filter_and_display_graph,
            callbackOnReturn=True,
            tooltip="Max term p-value",
        )

        box = gui.widgetBox(box, "Significance test")

        gui.radioButtonsInBox(
            box,
            self,
            "prob_func",
            ["Binomial", "Hypergeometric"],
            tooltips=["Use binomial distribution test", "Use hypergeometric distribution test"],
            callback=self.__invalidate,
        )  # TODO: only update the p values
        box = gui.widgetBox(self.filterTab, "Evidence codes in annotation", addSpace=True)
        self.evidenceCheckBoxDict = {}
        for etype in go.evidence_types_ordered:
            ecb = QCheckBox(etype, toolTip=go.evidence_types[etype], checked=self.use_evidence_type[etype])
            ecb.toggled.connect(self.__on_evidence_changed)
            box.layout().addWidget(ecb)
            self.evidenceCheckBoxDict[etype] = ecb

        # Select tab
        self.selectTab = gui.createTabPage(self.tabs, "Select")
        box = gui.radioButtonsInBox(
            self.selectTab,
            self,
            "selection_direct_annotation",
            ["Directly or Indirectly", "Directly"],
            box="Annotated genes",
            callback=self.example_selection,
        )

        box = gui.widgetBox(self.selectTab, "Output", addSpace=True)
        gui.radioButtonsInBox(
            box,
            self,
            "selection_disjoint",
            btnLabels=["All selected genes", "Term-specific genes", "Common term genes"],
            tooltips=[
                "Outputs genes annotated to all selected GO terms",
                "Outputs genes that appear in only one of selected GO terms",
                "Outputs genes common to all selected GO terms",
            ],
            callback=self.example_selection,
        )

        # ListView for DAG, and table for significant GOIDs
        self.DAGcolumns = ['GO term', 'Cluster', 'Reference', 'p-value', 'FDR', 'Genes', 'Enrichment']

        self.splitter = QSplitter(Qt.Vertical, self.mainArea)
        self.mainArea.layout().addWidget(self.splitter)

        # list view
        self.listView = GOTreeWidget(self.splitter)
        self.listView.setSelectionMode(QTreeView.ExtendedSelection)
        self.listView.setAllColumnsShowFocus(1)
        self.listView.setColumnCount(len(self.DAGcolumns))
        self.listView.setHeaderLabels(self.DAGcolumns)

        self.listView.header().setSectionsClickable(True)
        self.listView.header().setSortIndicatorShown(True)
        self.listView.header().setSortIndicator(self.DAGcolumns.index('p-value'), Qt.AscendingOrder)
        self.listView.setSortingEnabled(True)
        self.listView.setItemDelegateForColumn(6, EnrichmentColumnItemDelegate(self))
        self.listView.setRootIsDecorated(True)

        self.listView.itemSelectionChanged.connect(self.view_selection_changed)

        # table of significant GO terms
        self.sigTerms = QTreeWidget(self.splitter)
        self.sigTerms.setColumnCount(len(self.DAGcolumns))
        self.sigTerms.setHeaderLabels(self.DAGcolumns)
        self.sigTerms.setSortingEnabled(True)
        self.sigTerms.setSelectionMode(QTreeView.ExtendedSelection)
        self.sigTerms.header().setSortIndicator(self.DAGcolumns.index('p-value'), Qt.AscendingOrder)
        self.sigTerms.setItemDelegateForColumn(6, EnrichmentColumnItemDelegate(self))

        self.sigTerms.itemSelectionChanged.connect(self.table_selection_changed)

        self.sigTableTermsSorted = []
        self.graph = {}
        self.originalGraph = None

        self.inputTab.layout().addStretch(1)
        self.filterTab.layout().addStretch(1)
        self.selectTab.layout().addStretch(1)

        class AnnotationSlot(SimpleNamespace):
            taxid = ...  # type: str
            name = ...  # type: str
            filename = ...  # type:str

            @staticmethod
            def parse_tax_id(f_name):
                return f_name.split('.')[0]

        try:
            remote_files = serverfiles.ServerFiles().listfiles(DOMAIN)
        except (ConnectTimeout, RequestException, ConnectionError):
            # TODO: Warn user about failed connection to the remote server
            remote_files = []

        self.available_annotations = [
            AnnotationSlot(
                taxid=AnnotationSlot.parse_tax_id(annotation_file),
                name=taxonomy.common_taxid_to_name(AnnotationSlot.parse_tax_id(annotation_file)),
                filename=FILENAME_ANNOTATION.format(AnnotationSlot.parse_tax_id(annotation_file)),
            )
            for _, annotation_file in set(remote_files + serverfiles.listfiles(DOMAIN))
            if annotation_file != FILENAME_ONTOLOGY
        ]
        self._executor = ThreadExecutor()

    def sizeHint(self):
        return QSize(1000, 700)

    def __on_evidence_changed(self):
        for etype, cb in self.evidenceCheckBoxDict.items():
            self.use_evidence_type[etype] = cb.isChecked()
        self.__invalidate()

    def clear(self):
        self.infoLabel.setText("No data on input\n")
        self.warning(0)
        self.warning(1)
        self.clear_graph()

        self.send("Data on Selected Genes", None)
        self.send("Enrichment Report", None)

    @check_table_annotation
    def set_dataset(self, data=None):
        self.closeContext()
        self.clear()
        self.Error.clear()
        # print(data)
        if data:
            self.input_data = data
            self.tax_id = data.attributes[TableAnnotation.tax_id]
            self.use_attr_names = data.attributes[TableAnnotation.gene_as_attr_name]
            if self.use_attr_names:
                self.gene_id_attribute = data.attributes[TableAnnotation.gene_id_attribute]
            else:
                self.gene_id_column = data.attributes[TableAnnotation.gene_id_column]
            self.annotation_index = None

            _c2i = {a.taxid: i for i, a in enumerate(self.available_annotations)}
            try:
                self.annotation_index = _c2i[self.tax_id]
            except KeyError:
                self.Error.serverfiles_unavailable()
                # raise ValueError('Taxonomy {} not supported.'.format(self.tax_id))
                return

            self.gene_info = gene.GeneInfo(self.tax_id)

            self.__invalidate()

    @check_table_annotation
    def set_reference_dataset(self, data=None):
        self.Error.clear()
        if data:
            self.ref_data = data
            self.ref_tax_id = data.attributes[TableAnnotation.tax_id]
            self.ref_use_attr_names = data.attributes[TableAnnotation.gene_as_attr_name]
            if self.ref_use_attr_names:
                self.ref_gene_id_attribute = data.attributes[TableAnnotation.gene_id_attribute]
            else:
                self.ref_gene_id_column = data.attributes[TableAnnotation.gene_id_column]

        self.referenceRadioBox.buttons[1].setDisabled(not bool(data))
        self.referenceRadioBox.buttons[1].setText("Reference set")
        if self.input_data is not None and self.use_reference_dataset:
            self.use_reference_dataset = 0 if not data else 1
            self.__invalidate()

    @Slot()
    def __invalidate(self):
        # Invalidate the current results or pending task and schedule an
        # update.
        self.__scheduletimer.start()
        if self.__state != State.Ready:
            self.__state |= State.Stale

        self.set_graph({})
        self.ref_genes = None
        self.input_genes = None

    def __invalidate_annotations(self):
        self.annotations = None
        self.loaded_annotation_code = None
        if self.input_data:
            self.infoLabel.setText("...\n")
        self.__invalidate()

    @Slot()
    def __update(self):
        self.__scheduletimer.stop()
        if self.input_data is None:
            return

        if self.__state & State.Running:
            self.__state |= State.Stale
        elif self.__state & State.Downloading:
            self.__state |= State.Stale
        elif self.__state & State.Ready:
            if self.__ensure_data():
                self.load()
                self.enrichment()
            else:
                assert self.__state & State.Downloading
                assert self.isBlocking()

    def __get_ref_genes(self):
        self.ref_genes = []

        if self.ref_use_attr_names:
            for variable in self.input_data.domain.attributes:
                self.ref_genes.append(str(variable.attributes.get(self.ref_gene_id_attribute, '?')))
        else:
            genes, _ = self.ref_data.get_column_view(self.ref_gene_id_column)
            self.ref_genes = [str(g) for g in genes]

    def __get_input_genes(self):
        self.input_genes = []

        if self.use_attr_names:
            for variable in self.input_data.domain.attributes:
                self.input_genes.append(str(variable.attributes.get(self.gene_id_attribute, '?')))
        else:
            genes, _ = self.input_data.get_column_view(self.gene_id_column)
            self.input_genes = [str(g) for g in genes]

    def filter_annotated_genes(self, genes):
        matchedgenes = self.annotations.get_gene_names_translator(genes).values()
        return matchedgenes, [gene for gene in genes if gene not in matchedgenes]

    def __start_download(self, files_list):
        # type: (List[Tuple[str, str]]) -> None
        task = EnsureDownloaded(files_list)
        task.progress.connect(self._progress_bar_set)

        f = self._executor.submit(task)
        fw = FutureWatcher(f, self)
        fw.finished.connect(self.__download_finish)
        fw.finished.connect(fw.deleteLater)
        fw.resultReady.connect(self.__invalidate)

        self.progressBarInit(processEvents=None)
        self.setBlocking(True)
        self.setStatusMessage("Downloading")
        self.__state = State.Downloading

    @Slot(Future)
    def __download_finish(self, result):
        # type: (Future[None]) -> None
        assert QThread.currentThread() is self.thread()
        assert result.done()
        self.setBlocking(False)
        self.setStatusMessage("")
        self.progressBarFinished(processEvents=False)
        try:
            result.result()
        except ConnectTimeout:
            logging.getLogger(__name__).error("Error:")
            self.error(
                2,
                "Internet connection error, unable to load data. "
                + "Check connection and create a new GO Browser widget.",
            )
        except RequestException as err:
            logging.getLogger(__name__).error("Error:")
            self.error(2, "Internet error:\n" + str(err))
        except BaseException as err:
            logging.getLogger(__name__).error("Error:")
            self.error(2, "Error:\n" + str(err))
            raise
        else:
            self.error(2)
        finally:
            self.__state = State.Ready

    def __ensure_data(self):
        # Ensure that all required database (ontology and annotations for
        # the current selected organism are present. If not start a download in
        # the background. Return True if all dbs are present and false
        # otherwise
        assert self.__state == State.Ready
        annotation = self.available_annotations[self.annotation_index]
        go_files = [fname for domain, fname in serverfiles.listfiles(DOMAIN)]
        files = []

        if annotation.filename not in go_files:
            files.append(("go", annotation.filename))

        if FILENAME_ONTOLOGY not in go_files:
            files.append((DOMAIN, FILENAME_ONTOLOGY))
        if files:
            self.__start_download(files)
            assert self.__state == State.Downloading
            return False
        else:
            return True

    def load(self):
        a = self.available_annotations[self.annotation_index]

        if self.ontology is None:
            self.ontology = go.Ontology()

        if a.taxid != self.loaded_annotation_code:
            self.annotations = None
            gc.collect()  # Force run garbage collection
            self.annotations = go.Annotations(a.taxid)
            self.loaded_annotation_code = a.taxid
            count = defaultdict(int)
            gene_sets = defaultdict(set)

            for anno in self.annotations.annotations:
                count[anno.evidence] += 1
                gene_sets[anno.evidence].add(anno.gene_id)
            for etype in go.evidence_types_ordered:
                ecb = self.evidenceCheckBoxDict[etype]
                ecb.setEnabled(bool(count[etype]))
                ecb.setText(etype + ": %i annots(%i genes)" % (count[etype], len(gene_sets[etype])))

    def enrichment(self):
        assert self.input_data is not None
        assert self.__state == State.Ready

        if not self.annotations.ontology:
            self.annotations.ontology = self.ontology

        self.error(1)
        self.warning([0, 1])

        self.__get_input_genes()
        self.input_genes = set(self.input_genes)
        self.known_input_genes = self.annotations.get_genes_with_known_annotation(self.input_genes)

        # self.clusterGenes = clusterGenes = self.annotations.map_to_ncbi_id(self.input_genes).values()

        self.infoLabel.setText(
            "%i unique genes on input\n%i (%.1f%%) genes with known annotations"
            % (
                len(self.input_genes),
                len(self.known_input_genes),
                100.0 * len(self.known_input_genes) / len(self.input_genes) if len(self.input_genes) else 0.0,
            )
        )

        if not self.use_reference_dataset or self.ref_data is None:
            self.information(2)
            self.information(1)
            self.ref_genes = set(self.gene_info.keys())

        elif self.ref_data is not None:
            self.__get_ref_genes()
            self.ref_genes = set(self.ref_genes)

            ref_count = len(self.ref_genes)
            if ref_count == 0:
                self.ref_genes = self.annotations.genes()
                self.referenceRadioBox.buttons[1].setText("Reference set")
                self.referenceRadioBox.buttons[1].setDisabled(True)
                self.information(
                    2, "Unable to extract gene names from reference dataset. " "Using entire genome for reference"
                )
                self.use_reference_dataset = 0
            else:
                self.referenceRadioBox.buttons[1].setText("Reference set ({} genes)".format(ref_count))
                self.referenceRadioBox.buttons[1].setDisabled(False)
                self.information(2)
        else:
            self.use_reference_dataset = 0
            self.ref_genes = []

        if not self.ref_genes:
            self.error(1, "No valid reference set")
            return {}

        evidences = []
        for etype in go.evidence_types_ordered:
            if self.use_evidence_type[etype]:
                evidences.append(etype)
        aspect = ['Process', 'Component', 'Function'][self.aspect_index]

        self.progressBarInit(processEvents=False)
        self.setBlocking(True)
        self.__state = State.Running

        if self.input_genes:
            f = self._executor.submit(
                self.annotations.get_enriched_terms,
                self.input_genes,
                self.ref_genes,
                evidences,
                aspect=aspect,
                prob=self.probFunctions[self.prob_func],
                use_fdr=False,
                progress_callback=methodinvoke(self, "_progress_bar_set", (float,)),
            )
            fw = FutureWatcher(f, parent=self)
            fw.done.connect(self.__on_enrichment_done)
            fw.done.connect(fw.deleteLater)
            return
        else:
            f = Future()
            f.set_result({})
            self.__on_enrichment_done(f)

    def __on_enrichment_done(self, results):
        # type: (Future[Dict[str, tuple]]) -> None
        self.progressBarFinished(processEvents=False)
        self.setBlocking(False)
        self.setStatusMessage("")
        if self.__state & State.Stale:
            self.__state = State.Ready
            self.__invalidate()
            return

        self.__state = State.Ready
        try:
            results = results.result()  # type: Dict[str, tuple]
        except Exception as ex:
            results = {}
            error = str(ex)
            self.error(1, error)

        if results:
            terms = list(results.items())
            fdr_vals = statistics.FDR([d[1] for _, d in terms])
            terms = [(key, d + (fdr,)) for (key, d), fdr in zip(terms, fdr_vals)]
            terms = dict(terms)

        else:
            terms = {}

        self.terms = terms

        if not self.terms:
            self.warning(0, "No enriched terms found.")
        else:
            self.warning(0)

        self.treeStructDict = {}
        ids = self.terms.keys()

        self.treeStructRootKey = None

        parents = {}
        for _id in ids:
            parents[_id] = {term for _, term in self.ontology[_id].related}

        children = {}
        for term in self.terms:
            children[term] = {id for id in ids if term in parents[id]}

        for term in self.terms:
            self.treeStructDict[term] = TreeNode(self.terms[term], children[term])
            if not self.ontology[term].related and not getattr(self.ontology[term], "is_obsolete", False):
                self.treeStructRootKey = term

        self.set_graph(terms)
        self._update_enrichment_report_output()
        self.commit()

    def _update_enrichment_report_output(self):
        terms = sorted(self.terms.items(), key=lambda item: item[1][1])
        # Create and send the enrichemnt report table.
        terms_domain = Orange.data.Domain(
            [],
            [],
            # All is meta!
            [
                Orange.data.StringVariable("GO Term Id"),
                Orange.data.StringVariable("GO Term Name"),
                Orange.data.ContinuousVariable("Cluster Frequency"),
                Orange.data.ContinuousVariable("Genes in Cluster", number_of_decimals=0),
                Orange.data.ContinuousVariable("Reference Frequency"),
                Orange.data.ContinuousVariable("Genes in Reference", number_of_decimals=0),
                Orange.data.ContinuousVariable("p-value"),
                Orange.data.ContinuousVariable("FDR"),
                Orange.data.ContinuousVariable("Enrichment"),
                Orange.data.StringVariable("Genes"),
            ],
        )

        terms = [
            [
                t_id,
                self.ontology[t_id].name,
                len(genes) / len(self.input_genes),
                len(genes),
                r_count / len(self.ref_genes),
                r_count,
                p_value,
                fdr,
                len(genes) / len(self.input_genes) * len(self.ref_genes) / r_count,
                ",".join(genes),
            ]
            for t_id, (genes, p_value, r_count, fdr) in terms
            if genes and r_count
        ]

        if terms:
            x = numpy.empty((len(terms), 0))
            m = numpy.array(terms, dtype=object)
            terms_table = Orange.data.Table.from_numpy(terms_domain, x, metas=m)
        else:
            terms_table = None
        self.send("Enrichment Report", terms_table)

    @Slot(float)
    def _progress_bar_set(self, value):
        assert QThread.currentThread() is self.thread()
        self.progressBarSet(value, processEvents=None)

    @Slot()
    def _progress_bar_finish(self):
        assert QThread.currentThread() is self.thread()
        self.progressBarFinished(processEvents=None)

    def filter_graph(self, graph):
        if self.filter_by_p_value_nofdr:
            graph = go.filter_by_p_value(graph, self.max_p_value_no_fdr)
        if self.filter_by_p_value:  # FDR
            graph = dict(filter(lambda item: item[1][3] <= self.max_p_value, graph.items()))
        if self.filter_by_num_of_instances:
            graph = dict(filter(lambda item: len(item[1][0]) >= self.min_num_of_instances, graph.items()))
        return graph

    def filter_and_display_graph(self):
        if self.input_data and self.originalGraph is not None:
            self.graph = self.filter_graph(self.originalGraph)
            if self.originalGraph and not self.graph:
                self.warning(1, "All found terms were filtered out.")
            else:
                self.warning(1)
            self.clear_graph()
            self.display_graph()

    def set_graph(self, graph=None):
        self.originalGraph = graph
        if graph:
            self.filter_and_display_graph()
        else:
            self.graph = {}
            self.clear_graph()

    def clear_graph(self):
        self.listView.clear()
        self.listViewItems = []
        self.sigTerms.clear()

    def display_graph(self):
        from_parent_dict = {}
        self.termListViewItemDict = {}
        self.listViewItems = []

        def enrichment(t):
            try:
                return len(t[0]) / t[2] * (len(self.ref_genes) / len(self.input_genes))
            except ZeroDivisionError:
                # TODO: find out why this happens
                return 0

        max_fold_enrichment = max([enrichment(term) for term in self.graph.values()] or [1])

        def add_node(term, parent, parent_display_node):
            if (parent, term) in from_parent_dict:
                return
            if term in self.graph:
                display_node = GOTreeWidgetItem(
                    self.ontology[term],
                    self.graph[term],
                    len(self.input_genes),
                    len(self.ref_genes),
                    max_fold_enrichment,
                    parent_display_node,
                )
                display_node.goId = term
                self.listViewItems.append(display_node)
                if term in self.termListViewItemDict:
                    self.termListViewItemDict[term].append(display_node)
                else:
                    self.termListViewItemDict[term] = [display_node]
                from_parent_dict[(parent, term)] = True
                parent = term
            else:
                display_node = parent_display_node

            for c in self.treeStructDict[term].children:
                add_node(c, parent, display_node)

        if self.treeStructDict:
            add_node(self.treeStructRootKey, None, self.listView)

        terms = self.graph.items()
        terms = sorted(terms, key=lambda item: item[1][1])
        self.sigTableTermsSorted = [t[0] for t in terms]

        self.sigTerms.clear()
        for i, (t_id, (genes, p_value, ref_count, fdr)) in enumerate(terms):
            item = GOTreeWidgetItem(
                self.ontology[t_id],
                (genes, p_value, ref_count, fdr),
                len(self.input_genes),
                len(self.ref_genes),
                max_fold_enrichment,
                self.sigTerms,
            )
            item.goId = t_id

        self.listView.expandAll()
        for i in range(5):
            self.listView.resizeColumnToContents(i)
            self.sigTerms.resizeColumnToContents(i)
        self.sigTerms.resizeColumnToContents(6)
        width = min(self.listView.columnWidth(0), 350)
        self.listView.setColumnWidth(0, width)
        self.sigTerms.setColumnWidth(0, width)

    def view_selection_changed(self):
        if self.selectionChanging:
            return

        self.selectionChanging = 1
        self.selectedTerms = []
        selected = self.listView.selectedItems()
        self.selectedTerms = list({lvi.term.id for lvi in selected})
        self.example_selection()
        self.selectionChanging = 0

    def table_selection_changed(self):
        if self.selectionChanging:
            return

        self.selectionChanging = 1
        self.selectedTerms = []
        selected_ids = {self.sigTerms.itemFromIndex(index).goId for index in self.sigTerms.selectedIndexes()}

        for i in range(self.sigTerms.topLevelItemCount()):
            item = self.sigTerms.topLevelItem(i)
            selected = item.goId in selected_ids
            term = item.goId

            if selected:
                self.selectedTerms.append(term)

            for lvi in self.termListViewItemDict[term]:
                try:
                    lvi.setSelected(selected)
                    if selected:
                        lvi.setExpanded(True)
                except RuntimeError:  # Underlying C/C++ object deleted
                    pass
        self.selectionChanging = 0
        self.example_selection()

    def example_selection(self):
        self.commit()

    def commit(self):
        if self.input_data is None or self.originalGraph is None or self.annotations is None:
            return
        if self.__state & State.Stale:
            return

        terms = set(self.selectedTerms)
        genes = reduce(operator.ior, (set(self.graph[term][0]) for term in terms), set())

        evidences = []
        for etype in go.evidence_types_ordered:
            if self.use_evidence_type[etype]:
                evidences.append(etype)

        all_terms = self.annotations.get_annotated_terms(
            genes, direct_annotation_only=self.selection_direct_annotation, evidence_codes=evidences
        )

        if self.selection_disjoint > 0:
            count = defaultdict(int)
            for term in self.selectedTerms:
                for g in all_terms.get(term, []):
                    count[g] += 1
            ccount = 1 if self.selection_disjoint == 1 else len(self.selectedTerms)
            selected_genes = [gene for gene, c in count.items() if c == ccount and gene in genes]
        else:
            selected_genes = reduce(operator.ior, (set(all_terms.get(term, [])) for term in self.selectedTerms), set())

        if self.use_attr_names:
            selected = [
                column
                for column in self.input_data.domain.attributes
                if self.gene_id_attribute in column.attributes
                and str(column.attributes[self.gene_id_attribute]) in set(selected_genes)
            ]

            domain = Orange.data.Domain(selected, self.input_data.domain.class_vars, self.input_data.domain.metas)
            new_data = self.input_data.from_table(domain, self.input_data)
            self.send("Data on Selected Genes", new_data)

        else:
            selected_rows = []
            for row_index, row in enumerate(self.input_data):
                gene_in_row = str(row[self.gene_id_column])
                if gene_in_row in self.input_genes and gene_in_row in selected_genes:
                    selected_rows.append(row_index)

                if selected_rows:
                    selected = self.input_data[selected_rows]
                else:
                    selected = None

                self.send("Data on Selected Genes", selected)

    def show_info(self):
        dialog = QDialog(self)
        dialog.setModal(False)
        dialog.setLayout(QVBoxLayout())
        label = QLabel(dialog)
        label.setText("Ontology:\n" + self.ontology.header if self.ontology else "Ontology not loaded!")
        dialog.layout().addWidget(label)

        label = QLabel(dialog)
        label.setText(
            "Annotations:\n" + self.annotations.header.replace("!", "")
            if self.annotations
            else "Annotations not loaded!"
        )
        dialog.layout().addWidget(label)
        dialog.show()

    def onDeleteWidget(self):
        """Called before the widget is removed from the canvas."""
        self.annotations = None
        self.ontology = None
        gc.collect()  # Force collection


def fmtp(score):
    return "%0.5f" % score if score > 10e-4 else "%0.1e" % score


def fmtpdet(score):
    return "%0.9f" % score if score > 10e-4 else "%0.5e" % score


class GOTreeWidgetItem(QTreeWidgetItem):
    def __init__(self, term, enrichment_result, n_cluster_genes, n_ref_genes, max_fold_enrichment, parent):
        super().__init__(parent)
        self.term = term
        self.enrichmentResult = enrichment_result
        self.nClusterGenes = n_cluster_genes
        self.nRefGenes = n_ref_genes
        self.maxFoldEnrichment = max_fold_enrichment

        querymapped, pvalue, refmappedcount, fdr = enrichment_result

        querymappedcount = len(querymapped)
        if refmappedcount > 0 and n_ref_genes > 0 and n_cluster_genes > 0:
            enrichment = (querymappedcount / refmappedcount) * (n_ref_genes / n_cluster_genes)
        else:
            enrichment = numpy.nan

        self.enrichment = enrichment

        self.setText(0, term.name)

        fmt = "%" + str(-int(math.log(max(n_cluster_genes, 1)))) + "i (%.2f%%)"
        self.setText(1, fmt % (querymappedcount, 100.0 * querymappedcount / (n_cluster_genes or 1)))

        fmt = "%" + str(-int(math.log(max(n_ref_genes, 1)))) + "i (%.2f%%)"
        self.setText(2, fmt % (refmappedcount, 100.0 * refmappedcount / (n_ref_genes or 1)))

        self.setText(3, fmtp(pvalue))
        self.setToolTip(3, fmtpdet(pvalue))
        self.setText(4, fmtp(fdr))  # FDR
        self.setToolTip(4, fmtpdet(fdr))
        self.setText(5, ", ".join(querymapped))
        self.setText(6, "%.2f" % enrichment)
        self.setToolTip(6, "%.2f" % enrichment)
        self.setToolTip(0, "<p>" + term.__repr__()[6:].strip().replace("\n", "<br>"))
        self.sortByData = [
            term.name,
            querymappedcount,
            refmappedcount,
            pvalue,
            fdr,
            ", ".join(querymapped),
            enrichment,
        ]

    def data(self, col, role):
        if role == Qt.UserRole:
            if self.maxFoldEnrichment > 0:
                return self.enrichment / self.maxFoldEnrichment
            else:
                return numpy.nan
        else:
            return super().data(col, role)

    def __lt__(self, other):
        col = self.treeWidget().sortColumn()
        return self.sortByData[col] < other.sortByData[col]


class EnrichmentColumnItemDelegate(QItemDelegate):
    def paint(self, painter, option, index):
        self.drawBackground(painter, option, index)
        value = index.data(Qt.UserRole)
        if isinstance(value, float) and numpy.isfinite(value):
            painter.save()
            painter.setBrush(QBrush(Qt.blue, Qt.SolidPattern))
            painter.drawRect(
                option.rect.x(), option.rect.y(), int(value * (option.rect.width() - 1)), option.rect.height() - 1
            )
            painter.restore()
        else:
            super().paint(painter, option, index)


if __name__ == "__main__":

    def test_main(argv=sys.argv):
        from AnyQt.QtWidgets import QApplication

        app = QApplication(list(argv))
        argv = app.arguments()
        if len(argv) > 1:
            data = Orange.data.Table(argv[1])
        else:
            data = None

        w = OWGOBrowser()
        w.show()
        w.raise_()
        w.setDataset(data)
        w.handleNewSignals()
        rval = app.exec_()
        w.setDataset(None)
        w.handleNewSignals()
        w.saveSettings()
        w.onDeleteWidget()
        return rval

    sys.exit(test_main())
