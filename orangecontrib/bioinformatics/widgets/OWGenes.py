""" Genes """
import sys
from typing import List
from functools import lru_cache

import numpy as np

from AnyQt.QtGui import QFont, QColor
from AnyQt.QtCore import Qt, QSize, QTimer, QVariant, QModelIndex, QAbstractTableModel
from AnyQt.QtWidgets import QStyle, QSplitter, QTableView, QHeaderView, QAbstractItemView

from Orange.data import Table, Domain, StringVariable, DiscreteVariable
from Orange.data import filter as table_filter
from Orange.widgets.gui import (
    LinkRole,
    LinkStyledItemDelegate,
    vBox,
    rubber,
    checkBox,
    comboBox,
    lineEdit,
    widgetBox,
    auto_commit,
    widgetLabel,
)
from Orange.widgets.utils import itemmodels
from Orange.widgets.widget import OWWidget
from Orange.widgets.settings import Setting, ContextSetting, DomainContextHandler
from Orange.widgets.utils.signals import Input, Output
from Orange.widgets.utils.concurrent import TaskState, ConcurrentWidgetMixin

from orangecontrib.bioinformatics.ncbi import taxonomy
from orangecontrib.bioinformatics.ncbi.gene import ENTREZ_ID, Gene, GeneMatcher
from orangecontrib.bioinformatics.widgets.utils.data import (
    TAX_ID,
    GENE_ID_COLUMN,
    GENE_ID_ATTRIBUTE,
    GENE_AS_ATTRIBUTE_NAME,
)

NCBI_DETAIL_LINK = 'http://www.ncbi.nlm.nih.gov/sites/entrez?Db=gene&Cmd=ShowDetailView&TermToSearch={}'
HEADER = [
    ['Input ID', 'Entrez ID', 'Name', 'Description', 'Synonyms', 'Other IDs'],
    ['input_identifier', 'gene_id', 'symbol', 'description', 'synonyms', 'db_refs'],
]


def run_gene_matcher(gene_matcher: GeneMatcher, state: TaskState):
    current_iter = 0
    max_iter = len(gene_matcher.genes)

    def callback():
        nonlocal current_iter
        current_iter += 1
        state.set_progress_value(100 * (current_iter / max_iter))

    state.set_status("Working ...")
    gene_matcher._progress_callback = callback
    gene_matcher.match_genes()


class GeneInfoModel(itemmodels.PyTableModel):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

        self.header_labels, self.gene_attributes = HEADER
        self.setHorizontalHeaderLabels(self.header_labels)

        try:
            # note: make sure gene_id is set in owgenes_header
            self.entrez_column_index = self.gene_attributes.index('gene_id')
        except ValueError:
            raise ValueError("Make sure 'gene_id' is set in header")

        self.genes = None
        self.table = None

        self.font = QFont()
        self.font.setUnderline(True)
        self.color = QColor(Qt.blue)

        @lru_cache(maxsize=10000)
        def _row_instance(row, column):
            return self[int(row)][int(column)]

        self._row_instance = _row_instance

    def initialize(self, list_of_genes):
        self.genes = list_of_genes
        self.__table_from_genes([gene for gene in list_of_genes if gene.gene_id])
        self.update_model()

    def rowCount(self, parent=QModelIndex()):
        return 0 if parent.isValid() else len(self._table)

    def columnCount(self, parent=QModelIndex()):
        return 0 if (parent.isValid() or self._table.size == 0) else self._table.shape[1]

    def clear(self):
        self.beginResetModel()
        self._table = np.array([[]])
        self.resetSorting()
        self._roleData.clear()
        self.endResetModel()

    def data(
        self,
        index,
        role,
        _str=str,
        _Qt_DisplayRole=Qt.DisplayRole,  # noqa: N803
        _Qt_EditRole=Qt.EditRole,
        _Qt_FontRole=Qt.FontRole,
        _Qt_ForegroundRole=Qt.ForegroundRole,
        _LinkRolee=LinkRole,
        _recognizedRoles=frozenset([Qt.DisplayRole, Qt.EditRole, Qt.FontRole, Qt.ForegroundRole, LinkRole]),
    ):

        if role not in _recognizedRoles:
            return None

        row, col = index.row(), index.column()
        if not 0 <= row <= self.rowCount():
            return None
        row = self.mapToSourceRows(row)

        try:
            # value = self[row][col]
            value = self._row_instance(row, col)
        except IndexError:
            return

        if role == Qt.DisplayRole:
            return QVariant(str(value))
        elif role == Qt.ToolTipRole:
            return QVariant(str(value))

        if col == self.entrez_column_index:
            if role == _Qt_ForegroundRole:
                return self.color
            elif role == _Qt_FontRole:
                return self.font
            elif role == _LinkRolee:
                return NCBI_DETAIL_LINK.format(value)

    def __table_from_genes(self, list_of_genes: List[Gene]) -> None:
        def to_list(gene: Gene) -> List[str]:
            _, header_tags = HEADER

            def parse_attribute(tag):
                gene_attr = getattr(gene, '{}'.format(tag))

                if isinstance(gene_attr, dict):
                    # note: db_refs are stored as dicts
                    gene_attr = (
                        ', '.join('{}: {}'.format(key, val) for (key, val) in gene_attr.items()) if gene_attr else ' '
                    )
                elif isinstance(gene_attr, list):
                    # note: synonyms are stored as lists
                    gene_attr = ', '.join(gene_attr) if gene_attr else ' '

                return gene_attr

            return [parse_attribute(tag) for tag in header_tags]

        self.table = np.asarray([to_list(gene) for gene in list_of_genes])

    def get_filtered_genes(self):
        return list(self._table[:, self.entrez_column_index]) if self._table.size else []

    def filter_table(self, filter_pattern: str):
        selection = np.full(self.table.shape, True)
        for search_word in filter_pattern.split():
            match_result = np.core.defchararray.find(np.char.lower(self.table), search_word.lower()) >= 0
            selection = selection & match_result
        return selection

    def update_model(self, filter_pattern=''):
        # clear cache if model changes
        self._row_instance.cache_clear()

        if self.table.size:
            self.wrap(self.table[self.filter_table(filter_pattern).any(axis=1), :])


class UnknownGeneInfoModel(itemmodels.PyListModel):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.header_labels = ['IDs from the input data without corresponding Entrez ID']
        self.genes = []

    def initialize(self, list_of_genes):
        self.genes = list_of_genes
        self.wrap([', '.join([gene.input_identifier for gene in list_of_genes if not gene.gene_id])])

    def data(self, index, role=Qt.DisplayRole):
        row = index.row()
        if role in [self.list_item_role, Qt.EditRole] and self._is_index_valid(index):
            return self[row]
        elif role == Qt.TextAlignmentRole:
            return Qt.AlignLeft | Qt.AlignTop
        elif self._is_index_valid(row):
            return self._other_data[row].get(role, None)

    def headerData(self, section, orientation, role=Qt.DisplayRole):

        if role == Qt.DisplayRole and orientation == Qt.Horizontal:
            return self.header_labels[section]
        return QAbstractTableModel.headerData(self, section, orientation, role)


class OWGenes(OWWidget, ConcurrentWidgetMixin):
    name = "Genes"
    description = "Tool for working with genes"
    icon = "../widgets/icons/OWGeneInfo.svg"
    priority = 40
    want_main_area = True

    selected_organism: int = Setting(11)
    search_pattern: str = Setting('')
    exclude_unmatched = Setting(True)
    replace_id_with_symbol = Setting(True)
    auto_commit = Setting(True)

    settingsHandler = DomainContextHandler()
    selected_gene_col = ContextSetting(None)
    use_attr_names = ContextSetting(True)

    replaces = ['orangecontrib.bioinformatics.widgets.OWGeneNameMatcher.OWGeneNameMatcher']

    class Inputs:
        data_table = Input("Data", Table)

    class Outputs:
        data_table = Output("Data", Table)
        gene_matcher_results = Output("Genes", Table)

    class Information(OWWidget.Information):
        pass

    def sizeHint(self):
        return QSize(1280, 960)

    def __init__(self):
        OWWidget.__init__(self)
        ConcurrentWidgetMixin.__init__(self)

        # ATTRIBUTES #
        self.target_database = ENTREZ_ID

        # input data
        self.input_data = None
        self.input_genes = None
        self.tax_id = None
        self.column_candidates = []

        # input options
        self.organisms = []

        # gene matcher
        self.gene_matcher = None

        # progress bar
        self.progress_bar = None

        self._timer = QTimer()
        self._timer.timeout.connect(self._apply_filter)
        self._timer.setSingleShot(True)

        # GUI SECTION #

        # Control area
        self.info_box = widgetLabel(widgetBox(self.controlArea, "Info", addSpace=True), 'No data on input.\n')

        organism_box = vBox(self.controlArea, 'Organism')
        self.organism_select_combobox = comboBox(
            organism_box, self, 'selected_organism', callback=self.on_input_option_change
        )

        self.get_available_organisms()
        self.organism_select_combobox.setCurrentIndex(self.selected_organism)

        box = widgetBox(self.controlArea, 'Gene IDs in the input data')
        self.gene_columns_model = itemmodels.DomainModel(valid_types=(StringVariable, DiscreteVariable))
        self.gene_column_combobox = comboBox(
            box,
            self,
            'selected_gene_col',
            label='Stored in data column',
            model=self.gene_columns_model,
            sendSelectedValue=True,
            callback=self.on_input_option_change,
        )

        self.attr_names_checkbox = checkBox(
            box,
            self,
            'use_attr_names',
            'Stored as feature (column) names',
            disables=[(-1, self.gene_column_combobox)],
            callback=self.on_input_option_change,
        )

        self.gene_column_combobox.setDisabled(bool(self.use_attr_names))

        output_box = vBox(self.controlArea, 'Output')

        # separator(output_box)
        # output_box.layout().addWidget(horizontal_line())
        # separator(output_box)
        self.exclude_radio = checkBox(
            output_box, self, 'exclude_unmatched', 'Exclude unmatched genes', callback=self.commit
        )

        self.replace_radio = checkBox(
            output_box, self, 'replace_id_with_symbol', 'Replace feature IDs with gene names', callback=self.commit
        )

        auto_commit(self.controlArea, self, "auto_commit", "&Commit", box=False)

        rubber(self.controlArea)

        # Main area
        self.filter = lineEdit(
            self.mainArea, self, 'search_pattern', 'Filter:', callbackOnType=True, callback=self.handle_filter_callback
        )
        # rubber(self.radio_group)
        self.mainArea.layout().addWidget(self.filter)

        # set splitter
        self.splitter = QSplitter()
        self.splitter.setOrientation(Qt.Vertical)

        self.table_model = GeneInfoModel()
        self.table_view = QTableView()
        self.table_view.setAlternatingRowColors(True)
        self.table_view.viewport().setMouseTracking(True)
        self.table_view.setSortingEnabled(True)
        self.table_view.setShowGrid(False)
        self.table_view.verticalHeader().hide()
        # self.table_view.horizontalHeader().setSectionResizeMode(QHeaderView.Stretch)

        self.unknown_model = UnknownGeneInfoModel()

        self.unknown_view = QTableView()
        self.unknown_view.setModel(self.unknown_model)
        self.unknown_view.verticalHeader().hide()
        self.unknown_view.setShowGrid(False)
        self.unknown_view.setSelectionMode(QAbstractItemView.NoSelection)
        self.unknown_view.horizontalHeader().setSectionResizeMode(QHeaderView.Stretch)

        self.splitter.addWidget(self.table_view)
        self.splitter.addWidget(self.unknown_view)

        self.splitter.setStretchFactor(0, 90)
        self.splitter.setStretchFactor(1, 10)

        self.mainArea.layout().addWidget(self.splitter)

    def handle_filter_callback(self):
        self._timer.stop()
        self._timer.start(500)

    def _apply_filter(self):
        # filter only if input data is present and model is populated
        if self.table_model.table is not None:
            self.table_model.update_model(filter_pattern=str(self.search_pattern))
            self.commit()

    def __reset_widget_state(self):
        self.table_view.clearSpans()
        self.table_view.setModel(None)
        self.table_model.clear()
        self.unknown_model.clear()
        self._update_info_box()

    def _update_info_box(self):

        if self.input_genes and self.gene_matcher:
            num_genes = len(self.gene_matcher.genes)
            known_genes = len(self.gene_matcher.get_known_genes())

            info_text = (
                '{} genes in input data\n'
                '{} genes match Entrez database\n'
                '{} genes with match conflicts\n'.format(num_genes, known_genes, num_genes - known_genes)
            )

        else:
            info_text = 'No data on input.'

        self.info_box.setText(info_text)

    def on_done(self, _):
        # update info box
        self._update_info_box()

        # set output options
        self.toggle_radio_options()

        # set known genes
        self.table_model.initialize(self.gene_matcher.genes)
        self.table_view.setModel(self.table_model)
        self.table_view.selectionModel().selectionChanged.connect(self.commit)
        self.table_view.setSelectionBehavior(QAbstractItemView.SelectRows)

        self.table_view.setItemDelegateForColumn(
            self.table_model.entrez_column_index, LinkStyledItemDelegate(self.table_view)
        )
        v_header = self.table_view.verticalHeader()
        option = self.table_view.viewOptions()
        size = self.table_view.style().sizeFromContents(QStyle.CT_ItemViewItem, option, QSize(20, 20), self.table_view)

        v_header.setDefaultSectionSize(size.height() + 2)
        v_header.setMinimumSectionSize(5)
        self.table_view.horizontalHeader().setStretchLastSection(True)

        # set unknown genes
        self.unknown_model.initialize(self.gene_matcher.genes)
        self.unknown_view.verticalHeader().setStretchLastSection(True)

        self._apply_filter()

    def get_available_organisms(self):
        available_organism = sorted(
            ((tax_id, taxonomy.name(tax_id)) for tax_id in taxonomy.common_taxids()), key=lambda x: x[1]
        )

        self.organisms = [tax_id[0] for tax_id in available_organism]
        self.organism_select_combobox.addItems([tax_id[1] for tax_id in available_organism])

    def gene_names_from_table(self):
        """ Extract and return gene names from `Orange.data.Table`.
        """
        self.input_genes = []
        if self.input_data:
            if self.use_attr_names:
                self.input_genes = [str(attr.name).strip() for attr in self.input_data.domain.attributes]
            else:
                if self.selected_gene_col is None:
                    self.selected_gene_col = self.gene_column_identifier()

                self.input_genes = [
                    str(e[self.selected_gene_col]) for e in self.input_data if not np.isnan(e[self.selected_gene_col])
                ]

    def _update_gene_matcher(self):
        self.gene_names_from_table()

        self.gene_matcher = GeneMatcher(self.get_selected_organism(), auto_start=False)
        self.gene_matcher.genes = self.input_genes
        # self.gene_matcher.organism = self.get_selected_organism()

    def get_selected_organism(self):
        return self.organisms[self.selected_organism]

    def _run(self):
        if self.gene_matcher is not None:
            self.start(run_gene_matcher, self.gene_matcher)

    def on_input_option_change(self):
        self.__reset_widget_state()
        self._update_gene_matcher()
        self._run()

    def gene_column_identifier(self):
        """
        Get most suitable column that stores genes. If there are
        several suitable columns, select the one with most unique
        values. Take the best one.
        """

        # candidates -> (variable, num of unique values)
        candidates = (
            (col, np.unique(self.input_data.get_column_view(col)[0]).size)
            for col in self.gene_columns_model
            if isinstance(col, DiscreteVariable) or isinstance(col, StringVariable)
        )

        best_candidate, _ = sorted(candidates, key=lambda x: x[1])[-1]
        return best_candidate

    def find_genes_location(self):
        """ Try locate the genes in the input data when we first load the data.

            Proposed rules:
                - when no suitable feature names are present, check the columns.
                - find the most suitable column, that is, the one with most unique values.

        """
        domain = self.input_data.domain
        if not domain.attributes:
            if self.selected_gene_col is None:
                self.selected_gene_col = self.gene_column_identifier()
                self.use_attr_names = False

    @Inputs.data_table
    def handle_input(self, data):
        self.closeContext()
        self.input_data = None
        self.input_genes = None
        self.__reset_widget_state()
        self.gene_columns_model.set_domain(None)
        self.selected_gene_col = None

        if data:
            self.input_data = data
            self.gene_columns_model.set_domain(self.input_data.domain)

            # check if input table has tax_id, human is used if tax_id is not found
            self.tax_id = str(self.input_data.attributes.get(TAX_ID, '9606'))
            # check for gene location. Default is that genes are attributes in the input table.
            self.use_attr_names = self.input_data.attributes.get(GENE_AS_ATTRIBUTE_NAME, self.use_attr_names)

            if self.tax_id in self.organisms and not self.selected_organism:
                self.selected_organism = self.organisms.index(self.tax_id)

            self.openContext(self.input_data.domain)
            self.find_genes_location()
            self.on_input_option_change()

    def commit(self):
        selection = self.table_view.selectionModel().selectedRows(self.table_model.entrez_column_index)

        selected_genes = [row.data() for row in selection]
        if not len(selected_genes):
            selected_genes = self.table_model.get_filtered_genes()

        gene_ids = self.get_target_ids()
        known_genes = [gid for gid in gene_ids if gid != '?']

        table = None
        gm_table = None
        if known_genes:
            # Genes are in rows (we have a column with genes).
            if not self.use_attr_names:

                if self.target_database in self.input_data.domain:
                    gene_var = self.input_data.domain[self.target_database]
                    metas = self.input_data.domain.metas
                else:
                    gene_var = StringVariable(self.target_database)
                    metas = self.input_data.domain.metas + (gene_var,)

                domain = Domain(self.input_data.domain.attributes, self.input_data.domain.class_vars, metas)

                table = self.input_data.transform(domain)
                col, _ = table.get_column_view(gene_var)
                col[:] = gene_ids

                # filter selected rows
                selected_genes_set = set(selected_genes)
                selected_rows = [
                    row_index for row_index, row in enumerate(table) if str(row[gene_var]) in selected_genes_set
                ]

                # handle table attributes
                table.attributes[TAX_ID] = self.get_selected_organism()
                table.attributes[GENE_AS_ATTRIBUTE_NAME] = False
                table.attributes[GENE_ID_COLUMN] = self.target_database
                table = table[selected_rows] if selected_rows else table

                if self.exclude_unmatched:
                    # create filter from selected column for genes
                    only_known = table_filter.FilterStringList(gene_var, known_genes)
                    # apply filter to the data
                    table = table_filter.Values([only_known])(table)

                self.Outputs.data_table.send(table)

            # genes are are in columns (genes are features).
            else:
                domain = self.input_data.domain.copy()
                table = self.input_data.transform(domain)

                for gene in self.gene_matcher.genes:
                    if gene.input_identifier in table.domain:

                        table.domain[gene.input_identifier].attributes[self.target_database] = (
                            str(gene.gene_id) if gene.gene_id else '?'
                        )

                        if self.replace_id_with_symbol:
                            try:
                                table.domain[gene.input_identifier].name = str(gene.symbol)
                            except AttributeError:
                                # TODO: missing gene symbol, need to handle this?
                                pass

                # filter selected columns
                selected_genes_set = set(selected_genes)
                selected = [
                    column
                    for column in table.domain.attributes
                    if self.target_database in column.attributes
                    and str(column.attributes[self.target_database]) in selected_genes_set
                ]

                output_attrs = table.domain.attributes

                if selected:
                    output_attrs = selected

                if self.exclude_unmatched:
                    known_genes_set = set(known_genes)
                    output_attrs = [
                        col for col in output_attrs if col.attributes[self.target_database] in known_genes_set
                    ]

                domain = Domain(output_attrs, table.domain.class_vars, table.domain.metas)

                table = table.from_table(domain, table)

                # handle table attributes
                table.attributes[TAX_ID] = self.get_selected_organism()
                table.attributes[GENE_AS_ATTRIBUTE_NAME] = True
                table.attributes[GENE_ID_ATTRIBUTE] = self.target_database

            gm_table = self.gene_matcher.to_data_table(selected_genes=selected_genes if selected_genes else None)

        self.Outputs.data_table.send(table)
        self.Outputs.gene_matcher_results.send(gm_table)

    def toggle_radio_options(self):
        self.replace_radio.setEnabled(bool(self.use_attr_names))

        if self.gene_matcher.genes:
            # enable checkbox if unknown genes are detected
            self.exclude_radio.setEnabled(len(self.gene_matcher.genes) != len(self.gene_matcher.get_known_genes()))
            self.exclude_unmatched = len(self.gene_matcher.genes) != len(self.gene_matcher.get_known_genes())

    def get_target_ids(self):
        return [str(gene.gene_id) if gene.gene_id else '?' for gene in self.gene_matcher.genes]


if __name__ == "__main__":

    def main_test():
        from AnyQt.QtWidgets import QApplication
        import sys

        app = QApplication([])
        w = OWGenes()
        if len(sys.argv) > 1:
            data = Table(sys.argv[1])
            w.handle_input(data)
        w.show()
        w.raise_()
        r = app.exec_()
        w.saveSettings()
        return r

    sys.exit(main_test())
