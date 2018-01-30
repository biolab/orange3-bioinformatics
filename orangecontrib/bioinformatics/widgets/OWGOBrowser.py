""" GO Enrichment Analysis """
import sys
import math
import gc
import webbrowser
import operator
import numpy
import Orange.data


from requests.exceptions import ConnectTimeout
from collections import defaultdict
from functools import reduce

from AnyQt.QtWidgets import (
    QTreeWidget, QTreeWidgetItem, QTreeView, QMenu, QCheckBox, QSplitter,
    QDialog, QVBoxLayout, QLabel, QItemDelegate
)
from AnyQt.QtGui import QBrush
from AnyQt.QtCore import Qt, QSize

from Orange.widgets import widget, gui, settings
from Orange.widgets.utils.datacaching import data_hints
from Orange.widgets.utils.concurrent import ThreadExecutor

from orangecontrib.bioinformatics import go
from orangecontrib.bioinformatics.ncbi import taxonomy
from orangecontrib.bioinformatics.utils import serverfiles, statistics
from orangecontrib.bioinformatics.widgets.utils.data import GENE_NAME, TAX_ID
from orangecontrib.bioinformatics.go.config import DOMAIN, FILENAME_ONTOLOGY


def isstring(var):
    return isinstance(var, Orange.data.StringVariable)


def listAvailable():
    taxids = taxonomy.common_taxids()
    essential = [(taxonomy.name(taxid), 'gene_association.{}'.format(taxid)) for taxid in taxids
                 if (DOMAIN, 'gene_association.{}'.format(taxid)) in serverfiles.ServerFiles().listfiles(DOMAIN)]
    return dict(essential)


class TreeNode(object):
    def __init__(self, value, children):
        self.value = value
        self.children = children


class GOTreeWidget(QTreeWidget):
    def contextMenuEvent(self, event):
        super().contextMenuEvent(event)
        term = self.itemAt(event.pos()).term
        self._currMenu = QMenu()
        self._currAction = self._currMenu.addAction("View term on AmiGO website")
        self._currAction.triggered.connect(lambda: self.BrowserAction(term))
        self._currMenu.popup(event.globalPos())

    def BrowserAction(self, term):
        if isinstance(term, go.Term):
            term = term.id
        webbrowser.open("http://amigo.geneontology.org/cgi-bin/amigo/term-details.cgi?term="+term)


class OWGOBrowser(widget.OWWidget):
    name = "GO Browser"
    description = "Enrichment analysis for Gene Ontology terms."
    icon = "../widgets/icons/OWGOBrowser.svg"
    priority = 7

    inputs = [("Cluster Data", Orange.data.Table,
               "setDataset", widget.Single + widget.Default),
              ("Reference Data", Orange.data.Table,
               "setReferenceDataset")]

    outputs = [("Data on Selected Genes", Orange.data.Table),
               ("Enrichment Report", Orange.data.Table)]

    settingsHandler = settings.DomainContextHandler()

    annotationIndex = settings.ContextSetting(0)
    geneAttrIndex = settings.ContextSetting(0)
    useAttrNames = settings.ContextSetting(False)
    useReferenceDataset = settings.Setting(False)
    aspectIndex = settings.Setting(0)

    useEvidenceType = settings.Setting(
        {et: True for et in go.evidenceTypesOrdered})

    filterByNumOfInstances = settings.Setting(False)
    minNumOfInstances = settings.Setting(1)
    filterByPValue = settings.Setting(True)
    maxPValue = settings.Setting(0.2)
    filterByPValue_nofdr = settings.Setting(False)
    maxPValue_nofdr = settings.Setting(0.01)
    probFunc = settings.Setting(0)

    selectionDirectAnnotation = settings.Setting(0)
    selectionDisjoint = settings.Setting(0)
    selectionAddTermAsClass = settings.Setting(0)

    Ready, Initializing, Running = 0, 1, 2

    def __init__(self, parent=None):
        super().__init__(self, parent)

        self.clusterDataset = None
        self.referenceDataset = None
        self.ontology = None
        self.annotations = None
        self.loadedAnnotationCode = "---"
        self.treeStructRootKey = None
        self.probFunctions = [statistics.Binomial(), statistics.Hypergeometric()]
        self.selectedTerms = []

        self.selectionChanging = 0
        self.__state = OWGOBrowser.Initializing

        self.annotationCodes = []

        #############
        # GUI
        #############
        self.tabs = gui.tabWidget(self.controlArea)
        # Input tab
        self.inputTab = gui.createTabPage(self.tabs, "Input")
        box = gui.widgetBox(self.inputTab, "Info")
        self.infoLabel = gui.widgetLabel(box, "No data on input\n")

        gui.button(box, self, "Ontology/Annotation Info",
                   callback=self.ShowInfo,
                   tooltip="Show information on loaded ontology and annotations")

        box = gui.widgetBox(self.inputTab, "Organism")
        self.annotationComboBox = gui.comboBox(
            box, self, "annotationIndex", items=self.annotationCodes,
            callback=self._updateEnrichment, tooltip="Select organism")

        genebox = gui.widgetBox(self.inputTab, "Gene Names")
        self.geneAttrIndexCombo = gui.comboBox(
            genebox, self, "geneAttrIndex", callback=self._updateEnrichment,
            tooltip="Use this attribute to extract gene names from input data")
        self.geneAttrIndexCombo.setDisabled(self.useAttrNames)

        self.useAttrNames_checkbox = gui.checkBox(genebox, self, "useAttrNames", "Use column names",
                                                  tooltip="Use column names for gene names",
                                                  callback=self._updateEnrichment)
        self.useAttrNames_checkbox.toggled[bool].connect(self.geneAttrIndexCombo.setDisabled)

        self.referenceRadioBox = gui.radioButtonsInBox(
            self.inputTab, self, "useReferenceDataset",
            ["Entire genome", "Reference set (input)"],
            tooltips=["Use entire genome for reference",
                      "Use genes from Referece Examples input signal as reference"],
            box="Reference", callback=self._updateEnrichment)

        self.referenceRadioBox.buttons[1].setDisabled(True)
        gui.radioButtonsInBox(
            self.inputTab, self, "aspectIndex",
            ["Biological process", "Cellular component", "Molecular function"],
            box="Aspect", callback=self._updateEnrichment)

        # Filter tab
        self.filterTab = gui.createTabPage(self.tabs, "Filter")
        box = gui.widgetBox(self.filterTab, "Filter GO Term Nodes")
        gui.checkBox(box, self, "filterByNumOfInstances", "Genes",
                     callback=self.FilterAndDisplayGraph, 
                     tooltip="Filter by number of input genes mapped to a term")
        ibox = gui.indentedBox(box)
        gui.spin(ibox, self, 'minNumOfInstances', 1, 100,
                 step=1, label='#:', labelWidth=15,
                 callback=self.FilterAndDisplayGraph,
                 callbackOnReturn=True,
                 tooltip="Min. number of input genes mapped to a term")

        gui.checkBox(box, self, "filterByPValue_nofdr", "p-value",
                     callback=self.FilterAndDisplayGraph,
                     tooltip="Filter by term p-value")

        gui.doubleSpin(gui.indentedBox(box), self, 'maxPValue_nofdr', 1e-8, 1,
                       step=1e-8,  label='p:', labelWidth=15,
                       callback=self.FilterAndDisplayGraph,
                       callbackOnReturn=True,
                       tooltip="Max term p-value")

        # use filterByPValue for FDR, as it was the default in prior versions
        gui.checkBox(box, self, "filterByPValue", "FDR",
                     callback=self.FilterAndDisplayGraph,
                     tooltip="Filter by term FDR")
        gui.doubleSpin(gui.indentedBox(box), self, 'maxPValue', 1e-8, 1,
                       step=1e-8,  label='p:', labelWidth=15,
                       callback=self.FilterAndDisplayGraph,
                       callbackOnReturn=True,
                       tooltip="Max term p-value")

        box = gui.widgetBox(box, "Significance test")

        gui.radioButtonsInBox(box, self, "probFunc", ["Binomial", "Hypergeometric"],
                              tooltips=["Use binomial distribution test",
                                        "Use hypergeometric distribution test"],
                              callback=self._updateEnrichment)
        box = gui.widgetBox(self.filterTab, "Evidence codes in annotation",
                              addSpace=True)
        self.evidenceCheckBoxDict = {}
        for etype in go.evidenceTypesOrdered:
            ecb = QCheckBox(
                etype, toolTip=go.evidenceTypes[etype],
                checked=self.useEvidenceType[etype])
            ecb.toggled.connect(self.__on_evidenceChanged)
            box.layout().addWidget(ecb)
            self.evidenceCheckBoxDict[etype] = ecb

        # Select tab
        self.selectTab = gui.createTabPage(self.tabs, "Select")
        box = gui.radioButtonsInBox(
            self.selectTab, self, "selectionDirectAnnotation",
            ["Directly or Indirectly", "Directly"],
            box="Annotated genes",
            callback=self.ExampleSelection)

        box = gui.widgetBox(self.selectTab, "Output", addSpace=True)
        gui.radioButtonsInBox(
            box, self, "selectionDisjoint",
            btnLabels=["All selected genes",
                       "Term-specific genes",
                       "Common term genes"],
            tooltips=["Outputs genes annotated to all selected GO terms",
                      "Outputs genes that appear in only one of selected GO terms", 
                      "Outputs genes common to all selected GO terms"],
            callback=[self.ExampleSelection,
                      self.UpdateAddClassButton])

        self.addClassCB = gui.checkBox(
            box, self, "selectionAddTermAsClass", "Add GO Term as class",
            callback=self.ExampleSelection)

        # ListView for DAG, and table for significant GOIDs
        self.DAGcolumns = ['GO term', 'Cluster', 'Reference', 'p-value',
                           'FDR', 'Genes', 'Enrichment']

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
        self.listView.setSortingEnabled(True)
        self.listView.setItemDelegateForColumn(
            6, EnrichmentColumnItemDelegate(self))
        self.listView.setRootIsDecorated(True)

        self.listView.itemSelectionChanged.connect(self.ViewSelectionChanged)

        # table of significant GO terms
        self.sigTerms = QTreeWidget(self.splitter)
        self.sigTerms.setColumnCount(len(self.DAGcolumns))
        self.sigTerms.setHeaderLabels(self.DAGcolumns)
        self.sigTerms.setSortingEnabled(True)
        self.sigTerms.setSelectionMode(QTreeView.ExtendedSelection)
        self.sigTerms.setItemDelegateForColumn(
            6, EnrichmentColumnItemDelegate(self))

        self.sigTerms.itemSelectionChanged.connect(self.TableSelectionChanged)

        self.sigTableTermsSorted = []
        self.graph = {}

        self.inputTab.layout().addStretch(1)
        self.filterTab.layout().addStretch(1)
        self.selectTab.layout().addStretch(1)

        self.__initialize_finish()

    def sizeHint(self):
        return QSize(1000, 700)

    def __initialize_finish(self):
        self.setBlocking(False)

        try:
            self.annotationFiles = listAvailable()
        except ConnectTimeout:
            self.error(2, "Internet connection error, unable to load data. " + \
                          "Check connection and create a new GO Browser widget.")
            self.filterTab.setEnabled(False)
            self.inputTab.setEnabled(False)
            self.selectTab.setEnabled(False)
            self.listView.setEnabled(False)
            self.sigTerms.setEnabled(False)
        else:
            self.annotationCodes = sorted(self.annotationFiles.keys())
            self.annotationComboBox.clear()
            self.annotationComboBox.addItems(self.annotationCodes)
            self.annotationComboBox.setCurrentIndex(self.annotationIndex)
            self.__state = OWGOBrowser.Ready

    def __on_evidenceChanged(self):
        for etype, cb in self.evidenceCheckBoxDict.items():
            self.useEvidenceType[etype] = cb.isChecked()
        self._updateEnrichment()

    def clear(self):
        self.infoLabel.setText("No data on input\n")
        self.warning(0)
        self.warning(1)
        self.geneAttrIndexCombo.clear()
        self.ClearGraph()

        self.send("Data on Selected Genes", None)
        self.send("Enrichment Report", None)

    def setDataset(self, data=None):
        if self.__state == OWGOBrowser.Initializing:
            self.__initialize_finish()

        self.clear()
        self.clusterDataset = data

        if data is not None:
            domain = data.domain
            allvars = domain.variables + domain.metas
            self.candidateGeneAttrs = [var for var in allvars if isstring(var)]

            self.geneAttrIndexCombo.clear()
            for var in self.candidateGeneAttrs:
                self.geneAttrIndexCombo.addItem(*gui.attributeItem(var))
            taxid = data_hints.get_hint(data, TAX_ID, '')

            if taxid is not None:
                filename = "gene_association.%s" % taxid
                if filename in self.annotationFiles.values():
                    self.annotationIndex = [i for i, name in enumerate(self.annotationCodes)
                                            if self.annotationFiles[name] == filename].pop()

            self.useAttrNames = data_hints.get_hint(data, GENE_NAME, default=self.useAttrNames)

            self.geneAttrIndex = min(self.geneAttrIndex, len(self.candidateGeneAttrs) - 1)
            if len(self.candidateGeneAttrs) == 0:
                self.useAttrNames = True
                self.geneAttrIndex = -1
            elif self.geneAttrIndex < len(self.candidateGeneAttrs):
                self.geneAttrIndex = len(self.candidateGeneAttrs) - 1

            self._updateEnrichment()

    def setReferenceDataset(self, data=None):
        self.referenceDataset = data
        self.referenceRadioBox.buttons[1].setDisabled(not bool(data))
        self.referenceRadioBox.buttons[1].setText("Reference set")
        if self.clusterDataset is not None and self.useReferenceDataset:
            self.useReferenceDataset = 0 if not data else 1
            graph = self.Enrichment()
            self.SetGraph(graph)
        elif self.clusterDataset:
            self.__updateReferenceSetButton()

    def handleNewSignals(self):
        super().handleNewSignals()

    def _updateEnrichment(self):
        if self.clusterDataset is not None and \
                self.__state == OWGOBrowser.Ready:
            pb = gui.ProgressBar(self, 100)
            self.Load(pb=pb)
            graph = self.Enrichment(pb=pb)
            self.SetGraph(graph)

    def __updateReferenceSetButton(self):
        allgenes, refgenes = None, []
        if self.referenceDataset:
            try:
                allgenes = self.genesFromTable(self.referenceDataset)
            except Exception:
                allgenes = []
        self.referenceRadioBox.buttons[1].setDisabled(not bool(allgenes))
        self.referenceRadioBox.buttons[1].setText("Reference set " + ("(%i genes, %i matched)" % (len(allgenes), len(refgenes)) if allgenes and refgenes else ""))

    def genesFromTable(self, data):
        if self.useAttrNames:
            genes = [v.name for v in data.domain.variables]
        else:
            attr = self.candidateGeneAttrs[min(self.geneAttrIndex, len(self.candidateGeneAttrs) - 1)]
            genes = [str(ex[attr]) for ex in data if not numpy.isnan(ex[attr])]
            if any("," in gene for gene in genes):
                self.information(0, "Separators detected in gene names. Assuming multiple genes per example.")
                genes = reduce(operator.iadd, (genes.split(",") for genes in genes), [])
        return genes

    def FilterAnnotatedGenes(self, genes):
        matchedgenes = self.annotations.get_gene_names_translator(genes).values()
        return matchedgenes, [gene for gene in genes if gene not in matchedgenes]

    def Load(self, pb=None):

        if self.__state == OWGOBrowser.Ready:
            go_files, tax_files = serverfiles.listfiles("GO"), serverfiles.listfiles("Taxonomy")
            calls = []
            pb, finish = (gui.ProgressBar(self, 0), True) if pb is None else (pb, False)
            count = 0
            if not tax_files:
                calls.append((taxonomy.DOMAIN, taxonomy.FILENAME))
                count += 1
            org = self.annotationCodes[min(self.annotationIndex, len(self.annotationCodes)-1)]
            if org != self.loadedAnnotationCode:
                count += 1
                if self.annotationFiles[org] not in go_files:
                    calls.append(("GO", self.annotationFiles[org]))
                    count += 1

            if FILENAME_ONTOLOGY not in go_files:
                calls.append((DOMAIN, FILENAME_ONTOLOGY))
                count += 1
            if not self.ontology:
                count += 1
            pb.iter += count * 100

            for args in calls:
                serverfiles.localpath_download(*args, **dict(callback=pb.advance))

            i = len(calls)
            if not self.ontology:
                self.ontology = go.Ontology(progress_callback=lambda value: pb.advance())
                i += 1

            if org != self.loadedAnnotationCode:
                self.annotations = None
                gc.collect()  # Force run garbage collection
                code = self.annotationFiles[org].split('.')[-1]
                self.annotations = go.Annotations(code, progress_callback=lambda value: pb.advance())
                i += 1
                self.loadedAnnotationCode = org
                count = defaultdict(int)
                geneSets = defaultdict(set)

                for anno in self.annotations.annotations:
                    count[anno.evidence] += 1
                    geneSets[anno.evidence].add(anno.gene_id)
                for etype in go.evidenceTypesOrdered:
                    ecb = self.evidenceCheckBoxDict[etype]
                    ecb.setEnabled(bool(count[etype]))
                    ecb.setText(etype + ": %i annots(%i genes)" % (count[etype], len(geneSets[etype])))
            if finish:
                pb.finish()

    def Enrichment(self, pb=None):
        assert self.clusterDataset is not None

        pb = gui.ProgressBar(self, 100) if pb is None else pb
        if not self.annotations.ontology:
            self.annotations.ontology = self.ontology

        self.error(1)
        self.warning([0, 1])

        if self.useAttrNames:
            clusterGenes = [v.name for v in self.clusterDataset.domain.attributes]
            self.information(0)
        elif 0 <= self.geneAttrIndex < len(self.candidateGeneAttrs):
            geneAttr = self.candidateGeneAttrs[self.geneAttrIndex]
            clusterGenes = [str(ex[geneAttr]) for ex in self.clusterDataset
                            if not numpy.isnan(ex[geneAttr])]
            if any("," in gene for gene in clusterGenes):
                self.information(0, "Separators detected in cluster gene names. Assuming multiple genes per example.")
                clusterGenes = reduce(operator.iadd, (genes.split(",") for genes in clusterGenes), [])
            else:
                self.information(0)
        else:
            self.error(1, "Failed to extract gene names from input dataset!")
            return {}

        genesSetCount = len(set(clusterGenes))

        self.clusterGenes = clusterGenes = self.annotations.map_to_ncbi_id(clusterGenes).values()

        self.infoLabel.setText("%i unique genes on input\n%i (%.1f%%) genes with known annotations" %
                               (genesSetCount, len(clusterGenes), 100.0*len(clusterGenes)/genesSetCount if genesSetCount else 0.0))

        referenceGenes = None
        if not self.useReferenceDataset or self.referenceDataset is None:
            self.information(2)
            self.information(1)
            referenceGenes = self.annotations.genes()

        elif self.referenceDataset is not None:
            if self.useAttrNames:
                referenceGenes = [v.name for v in self.referenceDataset.domain.attributes]
                self.information(1)
            elif geneAttr in (self.referenceDataset.domain.variables +
                              self.referenceDataset.domain.metas):
                referenceGenes = [str(ex[geneAttr]) for ex in self.referenceDataset
                                  if not numpy.isnan(ex[geneAttr])]
                if any("," in gene for gene in clusterGenes):
                    self.information(1, "Separators detected in reference gene names. Assuming multiple genes per example.")
                    referenceGenes = reduce(operator.iadd, (genes.split(",") for genes in referenceGenes), [])
                else:
                    self.information(1)
            else:
                self.information(1)
                referenceGenes = None

            if referenceGenes is None:
                referenceGenes = list(self.annotations.genes())
                self.referenceRadioBox.buttons[1].setText("Reference set")
                self.referenceRadioBox.buttons[1].setDisabled(True)
                self.information(2, "Unable to extract gene names from reference dataset. Using entire genome for reference")
                self.useReferenceDataset = 0
            else:
                refc = len(referenceGenes)
                # referenceGenes = self.annotations.get_gene_names_translator(referenceGenes).values()
                self.referenceRadioBox.buttons[1].setText("Reference set (%i genes, %i matched)" % (refc, len(referenceGenes)))
                self.referenceRadioBox.buttons[1].setDisabled(False)
                self.information(2)
        else:
            self.useReferenceDataset = 0
        if not referenceGenes:
            self.error(1, "No valid reference set")
            return {}

        self.referenceGenes = referenceGenes
        evidences = []
        for etype in go.evidenceTypesOrdered:
            if self.useEvidenceType[etype]:
                evidences.append(etype)
        aspect = ['Process', 'Component', 'Function'][self.aspectIndex]

        if clusterGenes:
            self.terms = terms = self.annotations.get_enriched_terms(
                clusterGenes, referenceGenes, evidences, aspect=aspect,
                prob=self.probFunctions[self.probFunc], use_fdr=False,
                progress_callback=lambda value: pb.advance())
            ids = []
            pvals = []
            for i, d in self.terms.items():
                ids.append(i)
                pvals.append(d[1])
            for i, fdr in zip(ids, statistics.FDR(pvals)):  # save FDR as the last part of the tuple
                terms[i] = tuple(list(terms[i]) + [fdr])
        else:
            self.terms = terms = {}
        if not self.terms:
            self.warning(0, "No enriched terms found.")
        else:
            self.warning(0)

        pb.finish()
        self.treeStructDict = {}
        ids = self.terms.keys()

        self.treeStructRootKey = None

        parents = {}
        for id in ids:
            parents[id] = set([term for _, term in self.ontology[id].related])

        children = {}
        for term in self.terms:
            children[term] = set([id for id in ids if term in parents[id]])

        for term in self.terms:
            self.treeStructDict[term] = TreeNode(self.terms[term], children[term])
            if not self.ontology[term].related and not getattr(self.ontology[term], "is_obsolete", False):
                self.treeStructRootKey = term

        # print(terms)
        return terms

    def FilterGraph(self, graph):
        if self.filterByPValue_nofdr:
            graph = go.filterByPValue(graph, self.maxPValue_nofdr)
        if self.filterByPValue:  # FDR
            graph = dict(filter(lambda item: item[1][3] <= self.maxPValue, graph.items()))
        if self.filterByNumOfInstances:
            graph = dict(filter(lambda item: len(item[1][0]) >= self.minNumOfInstances, graph.items()))
        return graph

    def FilterAndDisplayGraph(self):
        if self.clusterDataset:
            self.graph = self.FilterGraph(self.originalGraph)
            if self.originalGraph and not self.graph:
                self.warning(1, "All found terms were filtered out.")
            else:
                self.warning(1)
            self.ClearGraph()
            self.DisplayGraph()

    def SetGraph(self, graph=None):
        self.originalGraph = graph
        if graph:
            self.FilterAndDisplayGraph()
        else:
            self.graph = {}
            self.ClearGraph()

    def ClearGraph(self):
        self.listView.clear()
        self.listViewItems=[]
        self.sigTerms.clear()

    def DisplayGraph(self):
        fromParentDict = {}
        self.termListViewItemDict = {}
        self.listViewItems = []
        enrichment = lambda t: len(t[0]) / t[2] * (len(self.referenceGenes) / len(self.clusterGenes))
        maxFoldEnrichment = max([enrichment(term) for term in self.graph.values()] or [1])

        def addNode(term, parent, parentDisplayNode):
            if (parent, term) in fromParentDict:
                return
            if term in self.graph:
                displayNode = GOTreeWidgetItem(self.ontology[term], self.graph[term], len(self.clusterGenes), len(self.referenceGenes), maxFoldEnrichment, parentDisplayNode)
                displayNode.goId = term
                self.listViewItems.append(displayNode)
                if term in self.termListViewItemDict:
                    self.termListViewItemDict[term].append(displayNode)
                else:
                    self.termListViewItemDict[term] = [displayNode]
                fromParentDict[(parent, term)] = True
                parent = term
            else:
                displayNode = parentDisplayNode

            for c in self.treeStructDict[term].children:
                addNode(c, parent, displayNode)

        if self.treeStructDict:
            addNode(self.treeStructRootKey, None, self.listView)

        terms = self.graph.items()
        terms = sorted(terms, key=lambda item: item[1][1])
        self.sigTableTermsSorted = [t[0] for t in terms]

        self.sigTerms.clear()
        for i, (t_id, (genes, p_value, refCount, fdr)) in enumerate(terms):
            item = GOTreeWidgetItem(self.ontology[t_id],
                                    (genes, p_value, refCount, fdr),
                                    len(self.clusterGenes),
                                    len(self.referenceGenes),
                                    maxFoldEnrichment,
                                    self.sigTerms)
            item.goId = t_id

        self.listView.expandAll()
        for i in range(5):
            self.listView.resizeColumnToContents(i)
            self.sigTerms.resizeColumnToContents(i)
        self.sigTerms.resizeColumnToContents(6)
        width = min(self.listView.columnWidth(0), 350)
        self.listView.setColumnWidth(0, width)
        self.sigTerms.setColumnWidth(0, width)

        # Create and send the enrichemnt report table.
        termsDomain = Orange.data.Domain(
            [], [],
            # All is meta!
            [Orange.data.StringVariable("GO Term Id"),
             Orange.data.StringVariable("GO Term Name"),
             Orange.data.ContinuousVariable("Cluster Frequency"),
             Orange.data.ContinuousVariable("Genes in Cluster", number_of_decimals=0),
             Orange.data.ContinuousVariable("Reference Frequency"),
             Orange.data.ContinuousVariable("Genes in Reference", number_of_decimals=0),
             Orange.data.ContinuousVariable("p-value"),
             Orange.data.ContinuousVariable("FDR"),
             Orange.data.ContinuousVariable("Enrichment"),
             Orange.data.StringVariable("Genes")])

        terms = [[t_id,
                  self.ontology[t_id].name,
                  len(genes) / len(self.clusterGenes),
                  len(genes),
                  r_count / len(self.referenceGenes),
                  r_count,
                  p_value,
                  fdr,
                  len(genes) / len(self.clusterGenes) * \
                  len(self.referenceGenes) / r_count,
                  ",".join(genes)
                  ]
                 for t_id, (genes, p_value, r_count, fdr) in terms]

        if terms:
            X = numpy.empty((len(terms), 0))
            M = numpy.array(terms, dtype=object)
            termsTable = Orange.data.Table.from_numpy(termsDomain, X, metas=M)
        else:
            termsTable = Orange.data.Table(termsDomain)
        self.send("Enrichment Report", termsTable)

    def ViewSelectionChanged(self):
        if self.selectionChanging:
            return

        self.selectionChanging = 1
        self.selectedTerms = []
        selected = self.listView.selectedItems()
        self.selectedTerms = list(set([lvi.term.id for lvi in selected]))
        self.ExampleSelection()
        self.selectionChanging = 0

    def TableSelectionChanged(self):
        if self.selectionChanging:
            return

        self.selectionChanging = 1
        self.selectedTerms = []
        selectedIds = set([self.sigTerms.itemFromIndex(index).goId for index in self.sigTerms.selectedIndexes()])

        for i in range(self.sigTerms.topLevelItemCount()):
            item = self.sigTerms.topLevelItem(i)
            selected = item.goId in selectedIds
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

        self.ExampleSelection()
        self.selectionChanging = 0

    def UpdateAddClassButton(self):
        self.addClassCB.setEnabled(self.selectionDisjoint == 1)

    def ExampleSelection(self):
        self.commit()

    def commit(self):
        if self.clusterDataset is None:
            return

        terms = set(self.selectedTerms)
        genes = reduce(operator.ior,
                       (set(self.graph[term][0]) for term in terms), set())

        evidences = []
        for etype in go.evidenceTypesOrdered:
            if self.useEvidenceType[etype]:
                #  if getattr(self, "useEvidence" + etype):
                evidences.append(etype)
        allTerms = self.annotations.get_annotated_terms(
            genes, direct_annotation_only=self.selectionDirectAnnotation,
            evidence_codes=evidences)

        if self.selectionDisjoint > 0:
            count = defaultdict(int)
            for term in self.selectedTerms:
                for g in allTerms.get(term, []):
                    count[g] += 1
            ccount = 1 if self.selectionDisjoint == 1 else len(self.selectedTerms)
            selectedGenes = [gene for gene, c in count.items()
                             if c == ccount and gene in genes]
        else:
            selectedGenes = reduce(
                operator.ior,
                (set(allTerms.get(term, [])) for term in self.selectedTerms),
                set())

        if self.useAttrNames:
            vars = [self.clusterDataset.domain[gene]
                    for gene in set(selectedGenes)]
            domain = Orange.data.Domain(
                vars, self.clusterDataset.domain.class_vars,
                self.clusterDataset.domain.metas)
            newdata = self.clusterDataset.from_table(domain, self.clusterDataset)

            self.send("Data on Selected Genes", newdata)

        elif self.candidateGeneAttrs:
            selectedExamples = []

            geneAttr = self.candidateGeneAttrs[min(self.geneAttrIndex, len(self.candidateGeneAttrs)-1)]

            if self.selectionDisjoint == 1:
                goVar = Orange.data.DiscreteVariable(
                    "GO Term", values=list(self.selectedTerms))
                newDomain = Orange.data.Domain(
                    self.clusterDataset.domain.variables, goVar,
                    self.clusterDataset.domain.metas)
                goColumn = []
            for i, ex in enumerate(self.clusterDataset):
                if not numpy.isnan(ex[geneAttr]) and any(gene in selectedGenes for gene in str(ex[geneAttr]).split(",")):
                    if self.selectionDisjoint == 1 and self.selectionAddTermAsClass:
                        terms = filter(lambda term: any(gene in self.graph[term][0]
                                                        for gene in str(ex[geneAttr]).split(",")), self.selectedTerms)
                        term = sorted(terms)[0]
                        goColumn.append(goVar.values.index(term))
                    selectedExamples.append(i)

            if selectedExamples:
                selectedExamples = self.clusterDataset[selectedExamples]
                if self.selectionDisjoint == 1 and self.selectionAddTermAsClass:
                    selectedExamples = Orange.data.Table.from_table(newDomain, selectedExamples)
                    view, issparse = selectedExamples.get_column_view(goVar)
                    assert not issparse
                    view[:] = goColumn
            else:
                selectedExamples = None

            self.send("Data on Selected Genes", selectedExamples)

    def ShowInfo(self):
        dialog = QDialog(self)
        dialog.setModal(False)
        dialog.setLayout(QVBoxLayout())
        label = QLabel(dialog)
        label.setText("Ontology:\n" + self.ontology.header
                      if self.ontology else "Ontology not loaded!")
        dialog.layout().addWidget(label)

        label = QLabel(dialog)
        label.setText("Annotations:\n" + self.annotations.header.replace("!", "")
                      if self.annotations else "Annotations not loaded!")
        dialog.layout().addWidget(label)
        dialog.show()

    def onDeleteWidget(self):
        """Called before the widget is removed from the canvas.
        """
        self.annotations = None
        self.ontology = None
        gc.collect()  # Force collection


fmtp = lambda score: "%0.5f" % score if score > 10e-4 else "%0.1e" % score
fmtpdet = lambda score: "%0.9f" % score if score > 10e-4 else "%0.5e" % score


class GOTreeWidgetItem(QTreeWidgetItem):
    def __init__(self, term, enrichmentResult, nClusterGenes, nRefGenes,
                 maxFoldEnrichment, parent):
        super().__init__(parent)
        self.term = term
        self.enrichmentResult = enrichmentResult
        self.nClusterGenes = nClusterGenes
        self.nRefGenes = nRefGenes
        self.maxFoldEnrichment = maxFoldEnrichment

        querymapped, pvalue, refmappedcount, fdr = enrichmentResult

        querymappedcount = len(querymapped)
        if refmappedcount > 0 and nRefGenes > 0 and nClusterGenes > 0:
            enrichment = (querymappedcount / refmappedcount) * (nRefGenes / nClusterGenes)
        else:
            enrichment = numpy.nan

        self.enrichment = enrichment  # = lambda t: len(t[0]) / t[2] * (nRefGenes / (nClusterGenes or 1))

        self.setText(0, term.name)

        fmt = "%" + str(-int(math.log(max(nClusterGenes, 1)))) + "i (%.2f%%)"
        self.setText(1, fmt % (querymappedcount,
                               100.0 * querymappedcount / (nClusterGenes or 1)))

        fmt = "%" + str(-int(math.log(max(nRefGenes, 1)))) + "i (%.2f%%)"
        self.setText(2, fmt % (refmappedcount,
                               100.0 * refmappedcount / (nRefGenes or 1)))

        self.setText(3, fmtp(pvalue))
        self.setToolTip(3, fmtpdet(pvalue))
        self.setText(4, fmtp(fdr))  # FDR
        self.setToolTip(4, fmtpdet(fdr))
        self.setText(5, ", ".join(querymapped))
        self.setText(6, "%.2f" % (enrichment))
        self.setToolTip(6, "%.2f" % (enrichment))
        self.setToolTip(0, "<p>" + term.__repr__()[6:].strip().replace("\n", "<br>"))
        self.sortByData = [term.name, querymappedcount, refmappedcount,
                           pvalue, fdr, ", ".join(querymapped), enrichment]

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
            painter.drawRect(option.rect.x(), option.rect.y(),
                             int(value * (option.rect.width() - 1)),
                             option.rect.height() - 1)
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
