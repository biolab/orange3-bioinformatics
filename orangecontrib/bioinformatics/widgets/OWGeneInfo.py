""" Gene Info Widget Display gene summary information from NCBI Gene database. """
import sys
import math
import Orange
import numpy as np

from itertools import chain
from collections import defaultdict, OrderedDict
from functools import partial, lru_cache

from AnyQt.QtWidgets import QTreeView
from AnyQt.QtGui import QFont, QColor
from AnyQt.QtCore import (
    Qt, QSize, QThread, QAbstractItemModel, QSortFilterProxyModel,
    QModelIndex, QItemSelection, QItemSelectionModel, Slot
)

from Orange.widgets import widget, gui, settings
from Orange.widgets.utils.concurrent import ThreadExecutor, Task, methodinvoke
from Orange.widgets.utils.signals import Output, Input

from orangecontrib.bioinformatics.ncbi import gene, taxonomy
from orangecontrib.bioinformatics.utils import serverfiles
from orangecontrib.bioinformatics.widgets.utils.data import TAX_ID, GENE_AS_ATTRIBUTE_NAME
from orangecontrib.bioinformatics.widgets.utils.data import append_columns


def ensure_downloaded(domain, filename, advance=None):
    serverfiles.localpath_download(domain, filename, callback=advance)


class TreeModel(QAbstractItemModel):

    def __init__(self, data, header, parent):
        QAbstractItemModel.__init__(self, parent)
        self._data = data
        self._dataDict = {}
        self._header = header
        self._roleData = {Qt.DisplayRole: self._data}
        self._roleData = partial(
            defaultdict,
            partial(defaultdict,
                    partial(defaultdict, lambda: None)))(self._roleData)

        self.setColumnLinks()

    def setColumnLinks(self):
        font = QFont()
        font.setUnderline(True)

        for i, link in enumerate(self._data):
            self._roleData[gui.LinkRole][i][0] = gene.NCBI_DETAIL_LINK.format(self._roleData[Qt.DisplayRole][i][0])
            self._roleData[Qt.FontRole][i][0] = font
            self._roleData[Qt.ForegroundRole][i][0] = QColor(Qt.blue)

            self._roleData[gui.LinkRole][i][1] = gene.ENSEMBL_DETAIL_LINK.format(self._roleData[Qt.DisplayRole][i][1])
            self._roleData[Qt.FontRole][i][1] = font
            self._roleData[Qt.ForegroundRole][i][1] = QColor(Qt.blue)

    def setRoleData(self, role, row, col, data):
        self._roleData[role][row][col] = data

    def data(self, index, role=Qt.DisplayRole):
        row, col = index.row(), index.column()
        return self._roleData[role][row][col]

    def index(self, row, col, parent=QModelIndex()):
        return self.createIndex(row, col, 0)

    def parent(self, index):
        return QModelIndex()

    def rowCount(self, index=QModelIndex()):
        if index.isValid():
            return 0
        else:
            return len(self._data)

    def columnCount(self, index=QModelIndex()):
        return len(self._header)

    def headerData(self, section, orientation, role=Qt.DisplayRole):
        if role == Qt.DisplayRole:
            return self._header[section]
        return None


@lru_cache(maxsize=2)
def get_ncbi_info(taxid):
    return gene.GeneInfo(taxid)


HEADER_SCHEMA = {label: index for index, label in enumerate(gene.GENE_INFO_HEADER_LABELS)}


def ncbi_info(taxid, genes, advance=None):
    ensure_downloaded(gene.DOMAIN, gene.FILENAME, advance)
    info = get_ncbi_info(taxid)

    gene_matcher = gene.GeneMatcher(str(taxid))
    gene_matcher.genes = genes
    gene_matcher.run_matcher()

    map_input_to_ensembl = OrderedDict({input_name: '' for input_name in genes})
    ret = []

    for gene_obj in gene_matcher.get_known_genes():
        if gene_obj.ncbi_id:
            gi = info.get_gene_by_id(gene_obj.ncbi_id)

            ensemble = ''
            if 'Ensembl' in gi.db_refs:
                ensemble = gi.db_refs['Ensembl']
                map_input_to_ensembl[gene_obj.input_name] = ensemble

            ret.append((gene_obj.input_name,
                        [str(gi.gene_id),
                         str(ensemble),
                         gi.symbol + " (%s)" % gene_obj.input_name if gene_obj.input_name != gi.symbol else gi.symbol,
                         gi.locus_tag or "",
                         gi.chromosome or "",
                         gi.description or "",
                         gi.synonyms,
                         gi.symbol_from_nomenclature_authority or ""])
                       )
        else:
            ret.append(None)

    return map_input_to_ensembl, ret


class OWGeneInfo(widget.OWWidget):
    name = "Gene Info"
    description = "Displays gene information from NCBI and other sources."
    icon = "../widgets/icons/OWGeneInfo.svg"
    priority = 5

    class Inputs:
        data = Input("Data", Orange.data.Table)

    class Outputs:
        selected_genes = Output("Selected Genes", Orange.data.Table)
        data = Output("Data", Orange.data.Table)

    settingsHandler = settings.DomainContextHandler()

    organism_index = settings.ContextSetting(0)
    taxid = settings.ContextSetting("9606")

    gene_attr = settings.ContextSetting(0)

    auto_commit = settings.Setting(False)
    search_string = settings.Setting("")

    useAttr = settings.ContextSetting(False)
    useAltSource = settings.ContextSetting(False)

    def __init__(self, parent=None, ):
        super().__init__(self, parent)

        self.selectionChangedFlag = False

        self.__initialized = False
        self.initfuture = None
        self.itemsfuture = None

        self.map_input_to_ensembl = None
        self.infoLabel = gui.widgetLabel(
            gui.widgetBox(self.controlArea, "Info", addSpace=True),
            "Initializing\n"
        )

        self.organisms = None
        self.organismBox = gui.widgetBox(
            self.controlArea, "Organism", addSpace=True)

        self.organismComboBox = gui.comboBox(
            self.organismBox, self, "organism_index",
            callback=self._onSelectedOrganismChanged)

        box = gui.widgetBox(self.controlArea, "Gene names", addSpace=True)
        self.geneAttrComboBox = gui.comboBox(
            box, self, "gene_attr",
            "Gene attribute", callback=self.updateInfoItems
        )
        self.geneAttrComboBox.setEnabled(not self.useAttr)

        self.geneAttrCheckbox = gui.checkBox(box, self, "useAttr", "Use column names",
                                             callback=self.updateInfoItems)
        self.geneAttrCheckbox.toggled[bool].connect(self.geneAttrComboBox.setDisabled)

        gui.auto_commit(self.controlArea, self, "auto_commit", "Commit")

        gui.rubber(self.controlArea)

        gui.lineEdit(self.mainArea, self, "search_string", "Filter",
                     callbackOnType=True, callback=self.searchUpdate)

        self.treeWidget = QTreeView(self.mainArea)

        self.treeWidget.setAlternatingRowColors(True)
        self.treeWidget.setSortingEnabled(True)
        self.treeWidget.setSelectionMode(QTreeView.ExtendedSelection)
        self.treeWidget.setUniformRowHeights(True)
        self.treeWidget.setRootIsDecorated(False)

        self.treeWidget.setItemDelegateForColumn(HEADER_SCHEMA['NCBI ID'], gui.LinkStyledItemDelegate(self.treeWidget))
        self.treeWidget.setItemDelegateForColumn(HEADER_SCHEMA['Ensembl ID'], gui.LinkStyledItemDelegate(self.treeWidget))

        self.treeWidget.viewport().setMouseTracking(True)
        self.mainArea.layout().addWidget(self.treeWidget)

        box = gui.widgetBox(self.mainArea, "", orientation="horizontal")
        gui.button(box, self, "Select Filtered", callback=self.selectFiltered)
        gui.button(box, self, "Clear Selection",
                   callback=self.treeWidget.clearSelection)

        self.geneinfo = []
        self.cells = []
        self.row2geneinfo = {}
        self.data = None

        # : (# input genes, # matches genes)
        self.matchedInfo = 0, 0

        self.setBlocking(True)
        self.executor = ThreadExecutor(self)

        self.progressBarInit()

        task = Task(
            function=partial(
                taxonomy.ensure_downloaded,
                callback=methodinvoke(self, "advance", ())
            )
        )

        task.resultReady.connect(self.initialize)
        task.exceptionReady.connect(self._onInitializeError)

        self.initfuture = self.executor.submit(task)

    def sizeHint(self):
        return QSize(1024, 720)

    @Slot()
    def advance(self):
        assert self.thread() is QThread.currentThread()
        self.progressBarSet(self.progressBarValue + 1,
                            processEvents=None)

    def _get_available_organisms(self):
        available_organism = sorted([(tax_id, taxonomy.name(tax_id)) for tax_id in taxonomy.common_taxids()],
                                    key=lambda x: x[1])

        self.organisms = [tax_id[0] for tax_id in available_organism]

        self.organismComboBox.addItems([tax_id[1] for tax_id in available_organism])

    def initialize(self):
        if self.__initialized:
            # Already initialized
            return
        self.__initialized = True

        self._get_available_organisms()
        self.organism_index = self.organisms.index(taxonomy.DEFAULT_ORGANISM)
        self.taxid = self.organisms[self.organism_index]

        self.infoLabel.setText("No data on input\n")
        self.initfuture = None

        self.setBlocking(False)
        self.progressBarFinished(processEvents=None)

    def _onInitializeError(self, exc):
        sys.excepthook(type(exc), exc, None)
        self.error(0, "Could not download the necessary files.")

    def _onSelectedOrganismChanged(self):
        assert 0 <= self.organism_index <= len(self.organisms)
        self.taxid = self.organisms[self.organism_index]

        if self.data is not None:
            self.updateInfoItems()

    @Inputs.data
    def setData(self, data=None):
        if not self.__initialized:
            self.initfuture.result()
            self.initialize()

        if self.itemsfuture is not None:
            raise Exception("Already processing")

        self.data = data

        if data is not None:
            self.geneAttrComboBox.clear()
            self.attributes = [attr for attr in data.domain.variables + data.domain.metas
                               if isinstance(attr, (Orange.data.StringVariable, Orange.data.DiscreteVariable))]

            for var in self.attributes:
                self.geneAttrComboBox.addItem(*gui.attributeItem(var))

            self.taxid = str(self.data.attributes.get(TAX_ID, ''))
            self.useAttr = self.data.attributes.get(GENE_AS_ATTRIBUTE_NAME, self.useAttr)

            self.gene_attr = min(self.gene_attr, len(self.attributes) - 1)

            if self.taxid in self.organisms:
                self.organism_index = self.organisms.index(self.taxid)

            self.updateInfoItems()
        else:
            self.clear()

    def updateInfoItems(self):
        self.warning(0)
        if self.data is None:
            return

        if self.useAttr:
            genes = [attr.name for attr in self.data.domain.attributes]
        elif self.attributes:
            attr = self.attributes[self.gene_attr]
            genes = [str(ex[attr]) for ex in self.data
                     if not math.isnan(ex[attr])]
        else:
            genes = []
        if not genes:
            self.warning(0, "Could not extract genes from input dataset.")

        self.warning(1)
        org = self.organisms[min(self.organism_index, len(self.organisms) - 1)]
        source_name, info_getter = ("NCBI Info", ncbi_info)

        self.error(0)

        self.progressBarInit()
        self.setBlocking(True)
        self.setEnabled(False)
        self.infoLabel.setText("Retrieving info records.\n")

        self.genes = genes

        task = Task(
            function=partial(
                info_getter, org, genes,
                advance=methodinvoke(self, "advance", ()))
        )
        self.itemsfuture = self.executor.submit(task)
        task.finished.connect(self._onItemsCompleted)

    def _onItemsCompleted(self):
        self.setBlocking(False)
        self.progressBarFinished()
        self.setEnabled(True)

        try:
            self.map_input_to_ensembl, geneinfo = self.itemsfuture.result()
        finally:
            self.itemsfuture = None

        self.geneinfo = geneinfo
        self.cells = cells = []
        self.row2geneinfo = {}

        for i, (input_name, gi) in enumerate(geneinfo):
            if gi:
                row = []
                for item in gi:
                    row.append(item)

                # parse synonyms
                row[HEADER_SCHEMA['Synonyms']] = ','.join(row[HEADER_SCHEMA['Synonyms']])
                cells.append(row)
                self.row2geneinfo[len(cells) - 1] = i

        model = TreeModel(cells, list(HEADER_SCHEMA.keys()), None)

        proxyModel = QSortFilterProxyModel(self)
        proxyModel.setSourceModel(model)
        self.treeWidget.setModel(proxyModel)
        self.treeWidget.selectionModel().selectionChanged.connect(self.commit)

        for i in range(len(HEADER_SCHEMA)):
            self.treeWidget.resizeColumnToContents(i)
            self.treeWidget.setColumnWidth(
                i, min(self.treeWidget.columnWidth(i), 200)
            )

        self.infoLabel.setText("%i genes\n%i matched NCBI's IDs" %
                               (len(self.genes), len(cells)))
        self.matchedInfo = len(self.genes), len(cells)

        if self.useAttr:
            new_data = self.data.from_table(self.data.domain, self.data)

            for gene_var in new_data.domain.attributes:
                 gene_var.attributes['Ensembl ID'] = str(self.map_input_to_ensembl[gene_var.name])

            self.Outputs.data.send(new_data)

        elif self.attributes:
            ensembl_ids = []
            for gene_name in self.data.get_column_view(self.attributes[self.gene_attr])[0]:
                if gene_name and gene_name in self.map_input_to_ensembl:
                    ensembl_ids.append(self.map_input_to_ensembl[gene_name])
                else:
                    ensembl_ids.append('')

            data_with_ensembl = append_columns(self.data,
                                               metas=[(Orange.data.StringVariable('Ensembl ID'), ensembl_ids)])
            self.Outputs.data.send(data_with_ensembl)

    def clear(self):
        self.infoLabel.setText("No data on input\n")
        self.treeWidget.setModel(
            TreeModel([], ["NCBI ID", "Symbol", "Locus Tag",
                           "Chromosome", "Description", "Synonyms",
                           "Nomenclature"], self.treeWidget))

        self.geneAttrComboBox.clear()
        self.Outputs.selected_genes.send(None)

    def commit(self):
        if self.data is None:
            self.Outputs.selected_genes.send(None)
            self.Outputs.data.send(None)
            return

        model = self.treeWidget.model()
        selection = self.treeWidget.selectionModel().selection()
        selection = model.mapSelectionToSource(selection)
        selectedRows = list(
            chain(*(range(r.top(), r.bottom() + 1) for r in selection))
        )
        model = model.sourceModel()

        selectedGeneids = [self.row2geneinfo[row] for row in selectedRows]
        selectedIds = [self.geneinfo[i][0] for i in selectedGeneids]
        selectedIds = set(selectedIds)
        gene2row = dict((self.geneinfo[self.row2geneinfo[row]][0], row)
                        for row in selectedRows)

        isselected = selectedIds.__contains__

        if selectedIds:

            if self.useAttr:
                attrs = [attr for attr in self.data.domain.attributes if isselected(attr.name)]
                domain = Orange.data.Domain(attrs, self.data.domain.class_vars, self.data.domain.metas)
                newdata = self.data.from_table(domain, self.data)

                self.Outputs.selected_genes.send(newdata)

            elif self.attributes:
                attr = self.attributes[self.gene_attr]
                gene_col = [attr.str_val(v) for v in self.data.get_column_view(attr)[0]]
                gene_col = [(i, name) for i, name in enumerate(gene_col) if isselected(name)]
                indices = [i for i, _ in gene_col]

                # SELECTED GENES OUTPUT
                selected_genes_metas = [Orange.data.StringVariable(name) for name in gene.GENE_INFO_HEADER_LABELS]
                selected_genes_domain = Orange.data.Domain(
                    self.data.domain.attributes, self.data.domain.class_vars,
                    self.data.domain.metas + tuple(selected_genes_metas))

                selected_genes_data = self.data.from_table(selected_genes_domain, self.data)[indices]

                model_rows = [gene2row[gene_name] for _, gene_name in gene_col]
                for col, meta in zip(range(model.columnCount()), selected_genes_metas):
                    col_data = [str(model.index(row, col).data(Qt.DisplayRole)) for row in model_rows]
                    col_data = np.array(col_data, dtype=object, ndmin=2).T
                    selected_genes_data[:, meta] = col_data

                if not len(selected_genes_data):
                    selected_genes_data = None

                self.Outputs.selected_genes.send(selected_genes_data)
        else:
            self.Outputs.selected_genes.send(None)

    def rowFiltered(self, row):
        searchStrings = self.search_string.lower().split()
        row = " ".join(self.cells[row]).lower()
        return not all([s in row for s in searchStrings])

    def searchUpdate(self):
        if not self.data:
            return
        searchStrings = self.search_string.lower().split()
        index = self.treeWidget.model().sourceModel().index
        mapFromSource = self.treeWidget.model().mapFromSource
        for i, row in enumerate(self.cells):
            row = " ".join(row).lower()
            self.treeWidget.setRowHidden(
                mapFromSource(index(i, 0)).row(),
                QModelIndex(),
                not all([s in row for s in searchStrings]))

    def selectFiltered(self):
        if not self.data:
            return
        itemSelection = QItemSelection()

        index = self.treeWidget.model().sourceModel().index
        mapFromSource = self.treeWidget.model().mapFromSource
        for i, row in enumerate(self.cells):
            if not self.rowFiltered(i):
                itemSelection.select(mapFromSource(index(i, 0)),
                                     mapFromSource(index(i, 0)))
        self.treeWidget.selectionModel().select(
            itemSelection,
            QItemSelectionModel.Select | QItemSelectionModel.Rows)

    def onAltSourceChange(self):
        self.updateInfoItems()

    def onDeleteWidget(self):
        # try to cancel pending tasks
        if self.initfuture:
            self.initfuture.cancel()
        if self.itemsfuture:
            self.itemsfuture.cancel()

        self.executor.shutdown(wait=False)
        super().onDeleteWidget()


def reportItemView(view):
    model = view.model()
    return reportItemModel(view, model)


def reportItemModel(view, model, index=QModelIndex()):
    if not index.isValid() or model.hasChildren(index):
        columnCount, rowCount = model.columnCount(index), model.rowCount(index)
        if not index.isValid():
            text = ('<table>\n<tr>' +
                    ''.join('<th>%s</th>' %
                            model.headerData(i, Qt.Horizontal, Qt.DisplayRole)
                            for i in range(columnCount)) +
                    '</tr>\n')
        else:
            pass
        text += ''.join('<tr>' +
                        ''.join('<td>' + reportItemModel(view, model, model.index(row, column, index)) +
                                '</td>' for column in range(columnCount)) +
                        '</tr>\n'
                        for row in range(rowCount)
                        if not view.isRowHidden(row, index))
        text += '</table>'
        return text
    else:
        variant = model.data(index, Qt.DisplayRole)
        return str(variant)


def test_main(argv=sys.argv):
    from AnyQt.QtWidgets import QApplication
    app = QApplication(argv)

    if len(argv) > 1:
        filename = argv[1]
    else:
        filename = "brown-selected"

    data = Orange.data.Table(filename)
    w = OWGeneInfo()
    w.setData(data)
    w.show()
    w.raise_()
    r = app.exec_()
    w.saveSettings()
    return r


if __name__ == "__main__":
    sys.exit(test_main())
