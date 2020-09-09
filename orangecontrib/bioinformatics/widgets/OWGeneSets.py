""" GeneSets """
from enum import IntEnum
from types import SimpleNamespace
from typing import Set, List, Tuple, Optional
from urllib.parse import urlparse

from AnyQt.QtGui import QColor, QStandardItem, QStandardItemModel
from AnyQt.QtCore import Qt, QSize
from AnyQt.QtWidgets import QTreeView, QHBoxLayout, QHeaderView

from Orange.data import Table, Domain
from Orange.data import filter as table_filter
from Orange.widgets.gui import LinkRole, LinkStyledItemDelegate, spin, vBox, lineEdit, widgetBox, auto_commit
from Orange.widgets.widget import Msg, OWWidget
from Orange.widgets.settings import Setting, SettingProvider
from Orange.widgets.utils.signals import Input, Output
from Orange.widgets.utils.concurrent import TaskState, ConcurrentWidgetMixin

from orangecontrib.bioinformatics.geneset import GeneSets
from orangecontrib.bioinformatics.widgets.utils.gui import FilterProxyModel, NumericalColumnDelegate
from orangecontrib.bioinformatics.widgets.components import GeneSetSelection
from orangecontrib.bioinformatics.widgets.utils.data import TableAnnotation, check_table_annotation


class Results(SimpleNamespace):
    items: List[QStandardItem] = []


def run(gene_sets: GeneSets, selected_gene_sets: List[Tuple[str, ...]], genes, state: TaskState) -> Results:
    results = Results()
    items = []
    step, steps = 0, len(gene_sets)

    if not genes:
        return results

    state.set_status('Calculating...')

    for gene_set in sorted(gene_sets):

        step += 1
        if step % (steps / 10) == 0:
            state.set_progress_value(100 * step / steps)

        if gene_set.hierarchy not in selected_gene_sets:
            continue

        if state.is_interruption_requested():
            return results

        matched_set = gene_set.genes & genes
        if len(matched_set) > 0:
            category_column = QStandardItem()
            term_column = QStandardItem()
            count_column = QStandardItem()
            genes_column = QStandardItem()

            category_column.setData(", ".join(gene_set.hierarchy), Qt.DisplayRole)
            term_column.setData(gene_set.name, Qt.DisplayRole)
            term_column.setData(gene_set.name, Qt.ToolTipRole)

            # there was some cases when link string was not empty string but not valid (e.g. "_")
            if gene_set.link and urlparse(gene_set.link).scheme:
                term_column.setData(gene_set.link, LinkRole)
                term_column.setForeground(QColor(Qt.blue))

            count_column.setData(matched_set, Qt.UserRole)
            count_column.setData(len(matched_set), Qt.DisplayRole)

            genes_column.setData(len(gene_set.genes), Qt.DisplayRole)
            genes_column.setData(set(gene_set.genes), Qt.UserRole)  # store genes to get then on output on selection

            items.append([count_column, genes_column, category_column, term_column])

    results.items = items
    return results


class Header(IntEnum):
    count = 0
    genes = 1
    category = 2
    term = 3

    @staticmethod
    def labels():
        return ['Count', 'Genes In Set', 'Category', 'Term']


class OWGeneSets(OWWidget, ConcurrentWidgetMixin):
    name = 'Gene Sets'
    description = ""
    icon = 'icons/OWGeneSets.svg'
    priority = 80
    want_main_area = True

    organism = Setting(None, schema_only=True)
    stored_gene_sets_selection = Setting([], schema_only=True)
    selected_rows = Setting([], schema_only=True)

    min_count: int
    min_count = Setting(5)

    use_min_count: bool
    use_min_count = Setting(True)

    auto_commit: bool
    auto_commit = Setting(False)

    search_pattern: str
    search_pattern = Setting('')

    # component settings
    gs_selection_component: SettingProvider = SettingProvider(GeneSetSelection)

    class Inputs:
        data = Input('Data', Table)
        custom_gene_sets = Input('Custom Gene Sets', Table)

    class Outputs:
        matched_genes = Output('Matched Genes', Table)

    class Warning(OWWidget.Warning):
        all_sets_filtered = Msg('All sets were filtered out.')

    class Error(OWWidget.Error):
        organism_mismatch = Msg('Organism in input data and custom gene sets does not match')
        cant_reach_host = Msg('Host orange.biolab.si is unreachable.')
        cant_load_organisms = Msg('No available organisms, please check your connection.')
        custom_gene_sets_table_format = Msg('Custom gene sets data must have genes represented as rows.')

    def __init__(self):
        super().__init__()
        # OWWidget.__init__(self)
        ConcurrentWidgetMixin.__init__(self)

        # Control area
        box = vBox(self.controlArea, True, margin=0)
        self.gs_selection_component: GeneSetSelection = GeneSetSelection(self, box)
        self.gs_selection_component.selection_changed.connect(self._on_selection_changed)

        # Main area
        self.filter_proxy_model = FilterProxyModel()
        self.filter_proxy_model.setFilterKeyColumn(Header.term)

        self.tree_view = QTreeView()
        self.tree_view.setAlternatingRowColors(True)
        self.tree_view.setSortingEnabled(True)
        self.tree_view.sortByColumn(Header.count, Qt.DescendingOrder)

        self.tree_view.setSelectionMode(QTreeView.ExtendedSelection)
        self.tree_view.setEditTriggers(QTreeView.NoEditTriggers)
        self.tree_view.viewport().setMouseTracking(True)
        self.tree_view.setItemDelegateForColumn(Header.term, LinkStyledItemDelegate(self.tree_view))
        self.tree_view.setItemDelegateForColumn(Header.genes, NumericalColumnDelegate(self))
        self.tree_view.setItemDelegateForColumn(Header.count, NumericalColumnDelegate(self))
        self.tree_view.setModel(self.filter_proxy_model)

        h_layout = QHBoxLayout()
        h_layout.setSpacing(100)
        h_widget = widgetBox(self.mainArea, orientation=h_layout)

        spin(
            h_widget,
            self,
            'min_count',
            0,
            1000,
            label='Count',
            tooltip='Minimum genes count',
            checked='use_min_count',
            callback=self.filter_view,
            callbackOnReturn=True,
            checkCallback=self.filter_view,
        )

        self.line_edit_filter = lineEdit(h_widget, self, 'search_pattern')
        self.line_edit_filter.setPlaceholderText('Filter gene sets ...')
        self.line_edit_filter.textChanged.connect(self.filter_view)

        self.mainArea.layout().addWidget(self.tree_view)
        self.tree_view.header().setSectionResizeMode(QHeaderView.ResizeToContents)

        self.commit_button = auto_commit(self.controlArea, self, 'auto_commit', '&Commit', box=False)

        self.input_data: Optional[Table] = None
        self.num_of_selected_genes: int = 0

    @property
    def tax_id(self) -> Optional[str]:
        if self.input_data:
            return self.input_data.attributes[TableAnnotation.tax_id]

    @property
    def gene_as_attr_name(self) -> Optional[bool]:
        if self.input_data:
            return self.input_data.attributes[TableAnnotation.gene_as_attr_name]

    @property
    def gene_location(self) -> Optional[str]:
        if not self.input_data:
            return

        if self.gene_as_attr_name:
            return self.input_data.attributes[TableAnnotation.gene_id_attribute]
        else:
            return self.input_data.attributes[TableAnnotation.gene_id_column]

    @property
    def input_genes(self) -> Set[str]:
        if not self.input_data:
            return set()

        if self.gene_as_attr_name:
            return {
                str(variable.attributes.get(self.gene_location, '?')) for variable in self.input_data.domain.attributes
            }
        else:
            return {str(g) for g in self.input_data.get_column_view(self.gene_location)[0]}

    def on_partial_result(self, _):
        pass

    def on_done(self, result: Results):
        model = QStandardItemModel()
        for item in result.items:
            model.appendRow(item)

        model.setSortRole(Qt.UserRole)
        model.setHorizontalHeaderLabels(Header.labels())

        self.filter_proxy_model.setSourceModel(model)
        self.tree_view.selectionModel().selectionChanged.connect(self.commit)
        self.filter_view()
        self.update_info_box()

    def on_exception(self, ex):
        # TODO: handle possible exceptions
        raise ex

    def onDeleteWidget(self):
        self.shutdown()
        super().onDeleteWidget()

    def _on_selection_changed(self):
        self.start(run, self.gs_selection_component.gene_sets, self.gs_selection_component.selection, self.input_genes)

    @Inputs.data
    @check_table_annotation
    def set_data(self, input_data: Table):
        self.Outputs.matched_genes.send(None)
        self.input_data = None
        self.num_of_selected_genes = 0

        if input_data:
            self.input_data = input_data
            self.gs_selection_component.initialize(self.tax_id)

        self.update_info_box()

    @Inputs.custom_gene_sets
    def handle_custom_gene_sets_input(self, custom_data):
        self.Outputs.matched_genes.send(None)

        if custom_data:
            self.gs_selection_component.initialize_custom_gene_sets(custom_data)
        else:
            self.gs_selection_component.initialize_custom_gene_sets(None)

        self.update_info_box()

    def commit(self):
        selection_model = self.tree_view.selectionModel()
        self.num_of_selected_genes = 0

        if selection_model:
            selection = selection_model.selectedRows(Header.count)
            self.selected_rows = [self.filter_proxy_model.mapToSource(sel).row() for sel in selection]

            if selection and self.input_genes:
                genes = [model_index.data(Qt.UserRole) for model_index in selection]
                output_genes = list(set.union(*genes))
                self.num_of_selected_genes = len(output_genes)

                if self.gene_as_attr_name:
                    selected = [
                        column
                        for column in self.input_data.domain.attributes
                        if self.gene_location in column.attributes
                        and str(column.attributes[self.gene_location]) in output_genes
                    ]

                    domain = Domain(selected, self.input_data.domain.class_vars, self.input_data.domain.metas)
                    new_data = self.input_data.from_table(domain, self.input_data)
                    self.Outputs.matched_genes.send(new_data)
                else:
                    # create filter from selected column for genes
                    only_known = table_filter.FilterStringList(self.gene_location, output_genes)
                    # apply filter to the data
                    data_table = table_filter.Values([only_known])(self.input_data)
                    self.Outputs.matched_genes.send(data_table)

        self.update_info_box()

    def update_info_box(self):
        input_string = ''
        input_number = ''

        if self.input_genes:
            input_string += '{} unique gene names on input.\n'.format(len(self.input_genes))
            input_number += str(len(self.input_genes))
            self.info.set_output_summary(
                str(self.num_of_selected_genes), '{} genes on output.\n'.format(self.num_of_selected_genes)
            )
        else:
            self.info.set_output_summary(self.info.NoOutput)

        if self.gs_selection_component.data:
            num_of_genes = self.gs_selection_component.num_of_genes
            num_of_sets = self.gs_selection_component.num_of_custom_sets
            input_number += f"{'' if input_number else '0'}|{num_of_genes}"
            input_string += '{} marker genes in {} sets\n'.format(num_of_genes, num_of_sets)

        if not input_number:
            self.info.set_input_summary(self.info.NoInput)
        else:
            self.info.set_input_summary(input_number, input_string)

    def create_filters(self):
        search_term: List[str] = self.search_pattern.lower().strip().split()

        filters = [
            FilterProxyModel.Filter(
                Header.term, Qt.DisplayRole, lambda value: all(fs in value.lower() for fs in search_term)
            )
        ]

        if self.use_min_count:
            filters.append(
                FilterProxyModel.Filter(Header.count, Qt.DisplayRole, lambda value: value >= self.min_count)
            )

        return filters

    def filter_view(self):
        filter_proxy: FilterProxyModel = self.filter_proxy_model
        model: QStandardItemModel = filter_proxy.sourceModel()

        if isinstance(model, QStandardItemModel):

            # apply filtering rules
            filter_proxy.set_filters(self.create_filters())

            if model.rowCount() and not filter_proxy.rowCount():
                self.Warning.all_sets_filtered()
            else:
                self.Warning.clear()

    def sizeHint(self):
        return QSize(800, 600)


if __name__ == "__main__":
    from Orange.widgets.utils.widgetpreview import WidgetPreview

    widget = WidgetPreview(OWGeneSets)
    widget.run()
