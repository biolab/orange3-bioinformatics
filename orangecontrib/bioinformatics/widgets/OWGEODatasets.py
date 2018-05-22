""" Gene Expression Omnibus datasets widget """
import glob
import os
import sys
import numpy
import urllib.request as urlrequest


from collections import defaultdict
from functools import partial, lru_cache

from AnyQt.QtWidgets import (
    QLineEdit, QSplitter, QTreeView, QTreeWidget, QTreeWidgetItem,
)
from AnyQt.QtCore import (
    Qt, QThread, QCoreApplication, QSortFilterProxyModel, QItemSelectionModel, Slot
)
from AnyQt.QtGui import QStandardItemModel, QStandardItem

from Orange.data import Table, Domain, StringVariable
from Orange.widgets import gui
from Orange.widgets.widget import OWWidget
from Orange.widgets.settings import Setting
from Orange.widgets.utils.concurrent import ThreadExecutor, Task, methodinvoke

from orangecontrib.bioinformatics.utils import serverfiles
from orangecontrib.bioinformatics.geo.utils import gds_ensure_downloaded
from orangecontrib.bioinformatics.geo.dataset import GDS, GDSInfo, DOMAIN
from orangecontrib.bioinformatics.widgets.utils.gui import TokenListCompleter
from orangecontrib.bioinformatics.widgets.utils.data import (
    GENE_AS_ATTRIBUTE_NAME, TAX_ID, GENE_ID_COLUMN, GENE_ID_ATTRIBUTE
)
from orangecontrib.bioinformatics.ncbi.gene import NCBI_ID


TextFilterRole = next(gui.OrangeUserRole)


class MySortFilterProxyModel(QSortFilterProxyModel):
    def __init__(self, parent=None):
        QSortFilterProxyModel.__init__(self, parent)
        self._filter_strings = []
        self._cache = {}
        self._cache_fixed = {}
        self._cache_prefix = {}
        self._row_text = {}

        # Create a cached version of _filteredRows
        self._filteredRows = lru_cache(100)(self._filteredRows)

    def setSourceModel(self, model):
        """ Set the source model for the filter.
        """
        self._filter_strings = []
        self._cache = {}
        self._cache_fixed = {}
        self._cache_prefix = {}
        self._row_text = {}
        QSortFilterProxyModel.setSourceModel(self, model)

    def addFilterFixedString(self, string, invalidate=True):
        """ Add `string` filter to the list of filters. If invalidate is
        True the filter cache will be recomputed.
        """
        self._filter_strings.append(string)
        all_rows = range(self.sourceModel().rowCount())
        row_text = [self.rowFilterText(row) for row in all_rows]
        self._cache[string] = [string in text for text in row_text]
        if invalidate:
            self.updateCached()
            self.invalidateFilter()

    def removeFilterFixedString(self, index=-1, invalidate=True):
        """ Remove the `index`-th filter string. If invalidate is True the
        filter cache will be recomputed.
        """
        string = self._filter_strings.pop(index)
        del self._cache[string]
        if invalidate:
            self.updateCached()
            self.invalidate()

    def setFilterFixedStrings(self, strings):
        """Set a list of string to be the new filters.
        """
        to_remove = set(self._filter_strings) - set(strings)
        to_add = set(strings) - set(self._filter_strings)
        for str in to_remove:
            self.removeFilterFixedString(
                self._filter_strings.index(str),
                invalidate=False)

        for str in to_add:
            self.addFilterFixedString(str, invalidate=False)
        self.updateCached()
        self.invalidate()

    def _filteredRows(self, filter_strings):
        """Return a dictionary mapping row indexes to True False values.

        .. note:: This helper function is wrapped in the __init__ method.

        """
        all_rows = range(self.sourceModel().rowCount())
        cache = self._cache
        return dict([(row, all([cache[str][row] for str in filter_strings]))
                     for row in all_rows])

    def updateCached(self):
        """Update the combined filter cache.
        """
        self._cache_fixed = self._filteredRows(
            tuple(sorted(self._filter_strings)))

    def setFilterFixedString(self, string):
        """Should this raise an error? It is not being used.
        """
        QSortFilterProxyModel.setFilterFixedString(self, string)

    def rowFilterText(self, row):
        """Return text for `row` to filter on.
        """
        f_role = self.filterRole()
        f_column = self.filterKeyColumn()
        s_model = self.sourceModel()
        data = s_model.data(s_model.index(row, f_column), f_role)
        return str(data)

    def filterAcceptsRow(self, row, parent):
        return self._cache_fixed.get(row, True)

    def lessThan(self, left, right):
        # TODO: Remove fixed column handling
        if left.column() == 1 and right.column():
            left_gds = str(left.data(Qt.DisplayRole))
            right_gds = str(right.data(Qt.DisplayRole))
            left_gds = left_gds.lstrip("GDS")
            right_gds = right_gds.lstrip("GDS")
            try:
                return int(left_gds) < int(right_gds)
            except ValueError:
                pass
        return QSortFilterProxyModel.lessThan(self, left, right)


def childiter(item):
    """ Iterate over the children of an QTreeWidgetItem instance.
    """
    for i in range(item.childCount()):
        yield item.child(i)


class OWGEODatasets(OWWidget):
    name = "GEO Data Sets"
    description = "Access to Gene Expression Omnibus data sets."
    icon = "icons/OWGEODatasets.svg"
    priority = 2

    inputs = []
    outputs = [("Expression Data", Table)]

    settingsList = ["outputRows", "mergeSpots", "gdsSelectionStates",
                    "splitterSettings", "currentGds", "autoCommit",
                    "datasetNames"]

    outputRows = Setting(True)
    mergeSpots = Setting(True)
    gdsSelectionStates = Setting({})
    currentGds = Setting(None)
    datasetNames = Setting({})
    splitterSettings = Setting(
        (b'\x00\x00\x00\xff\x00\x00\x00\x00\x00\x00\x00\x02\x00\x00\x01\xea\x00\x00\x00\xd7\x01\x00\x00\x00\x07\x01\x00\x00\x00\x02',
         b'\x00\x00\x00\xff\x00\x00\x00\x00\x00\x00\x00\x02\x00\x00\x01\xb5\x00\x00\x02\x10\x01\x00\x00\x00\x07\x01\x00\x00\x00\x01')
    )

    autoCommit = Setting(False)

    def __init__(self, parent=None, signalManager=None, name=" GEO Data Sets"):
        OWWidget.__init__(self, parent, signalManager, name)

        self.selectionChanged = False
        self.filterString = ""
        self.datasetName = ""

        ## GUI
        box = gui.widgetBox(self.controlArea, "Info", addSpace=True)
        self.infoBox = gui.widgetLabel(box, "Initializing\n\n")

        box = gui.widgetBox(self.controlArea, "Output", addSpace=True)
        gui.radioButtonsInBox(box, self, "outputRows",
                              ["Genes in rows", "Samples in rows"], "Rows",
                              callback=self.commitIf)

        gui.checkBox(box, self, "mergeSpots", "Merge spots of same gene",
                     callback=self.commitIf)

        gui.separator(box)
        self.nameEdit = gui.lineEdit(
            box, self, "datasetName", "Data set name",
            tooltip="Override the default output data set name",
            callback=self.onNameEdited
        )
        self.nameEdit.setPlaceholderText("")

        if sys.version_info < (3, ):
            box = gui.widgetBox(self.controlArea, "Commit", addSpace=True)
            self.commitButton = gui.button(
                box, self, "Commit", callback=self.commit)
            cb = gui.checkBox(box, self, "autoCommit", "Commit on any change")
            gui.setStopper(self, self.commitButton, cb, "selectionChanged",
                           self.commit)
        else:
            gui.auto_commit(self.controlArea, self, "autoCommit", "Commit",
                            box="Commit")
            self.commitIf = self.commit

        gui.rubber(self.controlArea)

        gui.widgetLabel(self.mainArea, "Filter")
        self.filterLineEdit = QLineEdit(
            textChanged=self.filter
        )
        self.completer = TokenListCompleter(
            self, caseSensitivity=Qt.CaseInsensitive
        )
        self.filterLineEdit.setCompleter(self.completer)

        self.mainArea.layout().addWidget(self.filterLineEdit)

        splitter = QSplitter(Qt.Vertical, self.mainArea)
        self.mainArea.layout().addWidget(splitter)
        self.treeWidget = QTreeView(splitter)

        self.treeWidget.setSelectionMode(QTreeView.SingleSelection)
        self.treeWidget.setRootIsDecorated(False)
        self.treeWidget.setSortingEnabled(True)
        self.treeWidget.setAlternatingRowColors(True)
        self.treeWidget.setUniformRowHeights(True)
        self.treeWidget.setEditTriggers(QTreeView.NoEditTriggers)

        linkdelegate = gui.LinkStyledItemDelegate(self.treeWidget)
        self.treeWidget.setItemDelegateForColumn(1, linkdelegate)
        self.treeWidget.setItemDelegateForColumn(8, linkdelegate)
        self.treeWidget.setItemDelegateForColumn(
            0, gui.IndicatorItemDelegate(self.treeWidget,
                                         role=Qt.DisplayRole))

        proxyModel = MySortFilterProxyModel(self.treeWidget)
        self.treeWidget.setModel(proxyModel)
        self.treeWidget.selectionModel().selectionChanged.connect(
            self.updateSelection
        )
        self.treeWidget.viewport().setMouseTracking(True)

        splitterH = QSplitter(Qt.Horizontal, splitter)

        box = gui.widgetBox(splitterH, "Description")
        self.infoGDS = gui.widgetLabel(box, "")
        self.infoGDS.setWordWrap(True)
        gui.rubber(box)

        box = gui.widgetBox(splitterH, "Sample Annotations")
        self.annotationsTree = QTreeWidget(box)
        self.annotationsTree.setHeaderLabels(
            ["Type (Sample annotations)", "Sample count"]
        )
        self.annotationsTree.setRootIsDecorated(True)
        box.layout().addWidget(self.annotationsTree)
        self.annotationsTree.itemChanged.connect(
            self.annotationSelectionChanged
        )
        self._annotationsUpdating = False
        self.splitters = splitter, splitterH

        for sp, setting in zip(self.splitters, self.splitterSettings):
            sp.splitterMoved.connect(self.splitterMoved)
            sp.restoreState(setting)

        self.searchKeys = ["dataset_id", "title", "platform_organism",
                           "description"]

        self.gds = []
        self.gds_info = None

        self.resize(1000, 600)

        self.setBlocking(True)
        self.setEnabled(False)
        self.progressBarInit()

        self._executor = ThreadExecutor()

        func = partial(get_gds_model,
                       methodinvoke(self, "_setProgress", (float,)))
        self._inittask = Task(function=func)
        self._inittask.finished.connect(self._initializemodel)
        self._executor.submit(self._inittask)

        self._datatask = None

    @Slot(float)
    def _setProgress(self, value):
        self.progressBarValue = value

    def _initializemodel(self):
        assert self.thread() is QThread.currentThread()
        model, self.gds_info, self.gds = self._inittask.result()
        model.setParent(self)

        proxy = self.treeWidget.model()
        proxy.setFilterKeyColumn(0)
        proxy.setFilterRole(TextFilterRole)
        proxy.setFilterCaseSensitivity(False)
        proxy.setFilterFixedString(self.filterString)

        proxy.setSourceModel(model)
        proxy.sort(0, Qt.DescendingOrder)

        self.progressBarFinished()
        self.setBlocking(False)
        self.setEnabled(True)

        filter_items = " ".join(
            gds[key] for gds in self.gds for key in self.searchKeys
        )
        tr_chars = ",.:;!?(){}[]_-+\\|/%#@$^&*<>~`"
        tr_table = str.maketrans(tr_chars, " " * len(tr_chars))
        filter_items = filter_items.translate(tr_table)

        filter_items = sorted(set(filter_items.split(" ")))
        filter_items = [item for item in filter_items if len(item) > 3]

        self.completer.setTokenList(filter_items)

        if self.currentGds:
            current_id = self.currentGds["dataset_id"]
            gdss = [(i, proxy.data(proxy.index(i, 1), Qt.DisplayRole))
                    for i in range(proxy.rowCount())]
            current = [i for i, data in gdss if data and data == current_id]
            if current:
                current_index = proxy.index(current[0], 0)
                self.treeWidget.selectionModel().select(
                    current_index,
                    QItemSelectionModel.Select | QItemSelectionModel.Rows
                )
                self.treeWidget.scrollTo(
                    current_index, QTreeView.PositionAtCenter)

        for i in range(8):
            self.treeWidget.resizeColumnToContents(i)

        self.treeWidget.setColumnWidth(
            1, min(self.treeWidget.columnWidth(1), 300))
        self.treeWidget.setColumnWidth(
            2, min(self.treeWidget.columnWidth(2), 200))

        self.updateInfo()

    def updateInfo(self):
        gds_info = self.gds_info
        text = ("%i datasets\n%i datasets cached\n" %
                (len(gds_info),
                 len(glob.glob(serverfiles.localpath("GEO") + "/GDS*"))))
        filtered = self.treeWidget.model().rowCount()
        if len(self.gds) != filtered:
            text += ("%i after filtering") % filtered
        self.infoBox.setText(text)

    def updateSelection(self, *args):
        current = self.treeWidget.selectedIndexes()
        mapToSource = self.treeWidget.model().mapToSource
        current = [mapToSource(index).row() for index in current]
        if current:
            self.currentGds = self.gds[current[0]]
            self.setAnnotations(self.currentGds)
            self.infoGDS.setText(self.currentGds.get("description", ""))
            self.nameEdit.setPlaceholderText(self.currentGds["title"])
            self.datasetName = \
                self.datasetNames.get(self.currentGds["dataset_id"], "")
        else:
            self.currentGds = None
            self.nameEdit.setPlaceholderText("")
            self.datasetName = ""

        self.commitIf()

    def setAnnotations(self, gds):
        self._annotationsUpdating = True
        self.annotationsTree.clear()

        annotations = defaultdict(set)
        subsetscount = {}
        for desc in gds["subsets"]:
            annotations[desc["type"]].add(desc["description"])
            subsetscount[desc["description"]] = str(len(desc["sample_id"]))

        for type, subsets in annotations.items():
            key = (gds["dataset_id"], type)
            subsetItem = QTreeWidgetItem(self.annotationsTree, [type])
            subsetItem.setFlags(subsetItem.flags() | Qt.ItemIsUserCheckable |
                                Qt.ItemIsTristate)
            subsetItem.setCheckState(
                0, self.gdsSelectionStates.get(key, Qt.Checked)
            )
            subsetItem.key = key
            for subset in subsets:
                key = (gds["dataset_id"], type, subset)
                item = QTreeWidgetItem(
                    subsetItem, [subset, subsetscount.get(subset, "")]
                )
                item.setFlags(item.flags() | Qt.ItemIsUserCheckable)
                item.setCheckState(
                    0, self.gdsSelectionStates.get(key, Qt.Checked)
                )
                item.key = key
        self._annotationsUpdating = False
        self.annotationsTree.expandAll()
        for i in range(self.annotationsTree.columnCount()):
            self.annotationsTree.resizeColumnToContents(i)

    def annotationSelectionChanged(self, item, column):
        if self._annotationsUpdating:
            return
        for i in range(self.annotationsTree.topLevelItemCount()):
            item = self.annotationsTree.topLevelItem(i)
            self.gdsSelectionStates[item.key] = item.checkState(0)
            for j in range(item.childCount()):
                child = item.child(j)
                self.gdsSelectionStates[child.key] = child.checkState(0)

    def filter(self):
        filter_string = self.filterLineEdit.text()
        proxyModel = self.treeWidget.model()
        if proxyModel:
            strings = filter_string.lower().strip().split()
            proxyModel.setFilterFixedStrings(strings)
            self.updateInfo()

    def selectedSamples(self):
        """
        Return the currently selected sample annotations.

        The return value is a list of selected (sample type, sample value)
        tuples.

        .. note:: if some Sample annotation type has no selected values.
                  this method will return all values for it.

        """
        samples = []
        unused_types = []
        used_types = []
        for stype in childiter(self.annotationsTree.invisibleRootItem()):
            selected_values = []
            all_values = []
            for sval in childiter(stype):
                value = (str(stype.text(0)), str(sval.text(0)))
                if self.gdsSelectionStates.get(sval.key, True):
                    selected_values.append(value)
                all_values.append(value)
            if selected_values:
                samples.extend(selected_values)
                used_types.append(str(stype.text(0)))
            else:
                # If no sample of sample type is selected we don't filter
                # on it.
                samples.extend(all_values)
                unused_types.append(str(stype.text(0)))

        return samples, used_types

    def commitIf(self):
        if self.autoCommit:
            self.commit()
        else:
            self.selectionChanged = True

    @Slot(int, int)
    def progressCompleted(self, value, total):
        if total > 0:
            self.progressBarSet(100. * value / total, processEvents=False)
        else:
            pass
            # TODO: report 'indeterminate progress'

    def commit(self):
        if self.currentGds:
            self.error(0)
            sample_type = None
            self.progressBarInit(processEvents=None)

            _, groups = self.selectedSamples()
            if len(groups) == 1 and self.outputRows:
                sample_type = groups[0]

            self.setEnabled(False)
            self.setBlocking(True)

            progress = methodinvoke(self, "progressCompleted", (int, int))

            def get_data(gds_id, report_genes, transpose, sample_type, title):
                gds_ensure_downloaded(gds_id, progress)
                gds = GDS(gds_id)
                data = gds.get_data(
                    report_genes=report_genes, transpose=transpose,
                    sample_type=sample_type
                )
                data.name = title
                return data

            get_data = partial(
                get_data, self.currentGds["dataset_id"],
                report_genes=self.mergeSpots,
                transpose=self.outputRows,
                sample_type=sample_type,
                title=self.datasetName or self.currentGds["title"]
            )
            self._datatask = Task(function=get_data)
            self._datatask.finished.connect(self._on_dataready)
            self._executor.submit(self._datatask)

    def _on_dataready(self):
        self.setEnabled(True)
        self.setBlocking(False)
        self.progressBarFinished(processEvents=False)

        try:
            data = self._datatask.result()
        except urlrequest.URLError as error:
            self.error(0, ("Error while connecting to the NCBI ftp server! "
                           "'%s'" % error))
            sys.excepthook(type(error), error, getattr(error, "__traceback__"))
            return
        finally:
            self._datatask = None

        data_name = data.name
        samples, _ = self.selectedSamples()

        self.warning(0)
        message = None
        from orangecontrib.bioinformatics.ncbi.gene import GeneMatcher

        gene_matcher = GeneMatcher(self.currentGds.get('taxid', ''))

        if self.outputRows:
            def samplesinst(ex):
                out = []
                for meta in data.domain.metas:
                    out.append((meta.name, ex[meta].value))

                if data.domain.class_var.name != 'class':
                    out.append((data.domain.class_var.name,
                                ex[data.domain.class_var].value))

                return out
            samples = set(samples)
            mask = [samples.issuperset(samplesinst(ex)) for ex in data]
            data = data[numpy.array(mask, dtype=bool)]
            gene_matcher.match_table_attributes(data)
            if len(data) == 0:
                message = "No samples with selected sample annotations."
        else:
            samples = set(samples)
            domain = Domain(
                [attr for attr in data.domain.attributes
                 if samples.issuperset(attr.attributes.items())],
                data.domain.class_var,
                data.domain.metas
            )
#             domain.addmetas(data.domain.getmetas())

            if len(domain.attributes) == 0:
                message = "No samples with selected sample annotations."
            stypes = set(s[0] for s in samples)
            for attr in domain.attributes:
                attr.attributes = dict(
                    (key, value) for key, value in attr.attributes.items()
                    if key in stypes
                )

            data = Table(domain, data)

            if 'gene' in data.domain:
                gene_column = data.domain['gene']
                gene_names = data.get_column_view(gene_column)[0]
                gene_matcher.genes = gene_names
                gene_matcher.run_matcher()

                domain_ids = Domain([], metas=[StringVariable(NCBI_ID)])
                data_ids = [[str(gene.ncbi_id) if gene.ncbi_id else '?'] for gene in gene_matcher.genes]
                table_ids = Table(domain_ids, data_ids)

                data = Table.concatenate([data, table_ids])

        if message is not None:
            self.warning(0, message)

        data.attributes[TAX_ID] = self.currentGds.get('taxid', '')
        data.attributes[GENE_AS_ATTRIBUTE_NAME] = bool(self.outputRows)

        if not bool(self.outputRows):
            data.attributes[GENE_ID_COLUMN] = NCBI_ID
        else:
            data.attributes[GENE_ID_ATTRIBUTE] = NCBI_ID

        data.name = data_name
        self.send("Expression Data", data)

        model = self.treeWidget.model().sourceModel()
        row = self.gds.index(self.currentGds)

        model.setData(model.index(row, 0),  " ", Qt.DisplayRole)

        self.updateInfo()
        self.selectionChanged = False


    def splitterMoved(self, *args):
        self.splitterSettings = [bytes(sp.saveState()) for sp in self.splitters]

    def send_report(self):
        self.report_items("GEO Dataset",
                          [("ID", self.currentGds['dataset_id']), ("Title", self.currentGds['title']),
                           ("Organism", self.currentGds['sample_organism'])])
        self.report_items("Data", [("Samples", self.currentGds['sample_count']),
                                   ("Features", self.currentGds['feature_count']),
                                   ("Genes", self.currentGds['gene_count'])])
        self.report_name("Sample annotations")
        subsets = defaultdict(list)
        for subset in self.currentGds['subsets']:
            subsets[subset['type']].append((subset['description'], len(subset['sample_id'])))
        self.report_html += "<ul>"
        for type in subsets:
            self.report_html += "<b>" + type + ":</b></br>"
            for desc, count in subsets[type]:
                self.report_html += 9 * "&nbsp" + "<b>{}:</b> {}</br>".format(desc, count)
        self.report_html += "</ul>"

    def onDeleteWidget(self):
        if self._inittask:
            self._inittask.future().cancel()
            self._inittask.finished.disconnect(self._initializemodel)
        if self._datatask:
            self._datatask.future().cancel()
            self._datatask.finished.disconnect(self._on_dataready)
        self._executor.shutdown(wait=False)

        super(OWGEODatasets, self).onDeleteWidget()

    def onNameEdited(self):
        if self.currentGds:
            gds_id = self.currentGds["dataset_id"]
            self.datasetNames[gds_id] = self.nameEdit.text()
            self.commitIf()


def get_gds_model(progress=lambda val: None):
    """
    Initialize and return a GDS datasets model.

    :param progress: A progress callback.
    :rval tuple:
        A tuple of (QStandardItemModel, GDSInfo, [GDS])

    .. note::
        The returned QStandardItemModel's thread affinity is set to
        the GUI thread.

    """
    progress(1)
    info = GDSInfo()
    search_keys = ["dataset_id", "title", "platform_organism", "description"]
    cache_dir = serverfiles.localpath(DOMAIN)
    gds_link = "http://www.ncbi.nlm.nih.gov/sites/GDSbrowser?acc={0}"
    pm_link = "http://www.ncbi.nlm.nih.gov/pubmed/{0}"
    gds_list = []

    def is_cached(gds):
        return os.path.exists(os.path.join(cache_dir, gds["dataset_id"]) +
                              ".soft.gz")

    def item(displayvalue, item_values={}):
        item = QStandardItem()
        item.setData(displayvalue, Qt.DisplayRole)
        for role, value in item_values.items():
            item.setData(value, role)
        return item

    def gds_to_row(gds):
        #: Text for easier full search.
        search_text = " | ".join([gds.get(key, "").lower()
                                  for key in search_keys])
        row = [
            item(" " if is_cached(gds) else "",
                 {TextFilterRole: search_text}),
            item(gds["dataset_id"],
                 {gui.LinkRole: gds_link.format(gds["dataset_id"])}),
            item(gds["title"]),
            item(gds["platform_organism"]),
            item(len(gds["samples"])),
            item(gds["feature_count"]),
            item(gds["gene_count"]),
            item(len(gds["subsets"])),
            item(gds.get("pubmed_id", ""),
                 {gui.LinkRole: pm_link.format(gds["pubmed_id"])
                            if gds.get("pubmed_id")
                            else None})
        ]
        return row

    model = QStandardItemModel()
    model.setHorizontalHeaderLabels(
        ["", "ID", "Title", "Organism", "Samples", "Features",
         "Genes", "Subsets", "PubMedID"]
    )
    progress(20)
    for gds in info.values():
        model.appendRow(gds_to_row(gds))

        gds_list.append(gds)

    progress(50)

    if QThread.currentThread() is not QCoreApplication.instance().thread():
        model.moveToThread(QCoreApplication.instance().thread())
    return model, info, gds_list


if __name__ == "__main__":
    def main_test():
        from AnyQt.QtWidgets import QApplication
        app = QApplication([])
        w = OWGEODatasets()
        w.show()
        w.raise_()
        r = app.exec_()
        w.saveSettings()
        return r

    sys.exit(main_test())
