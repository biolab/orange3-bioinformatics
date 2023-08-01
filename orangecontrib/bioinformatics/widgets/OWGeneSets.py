""" GeneSets """
from typing import Set, List, Tuple, Optional
from urllib.parse import urlparse

from AnyQt.QtCore import Qt
from AnyQt.QtWidgets import QTableView, QHBoxLayout

from Orange.data import Table, Domain
from Orange.data import filter as table_filter
from Orange.widgets.gui import (
    LinkRole,
    spin,
    vBox,
    deferred,
    lineEdit,
    widgetBox,
    auto_commit,
)
from Orange.widgets.widget import Msg, OWWidget
from Orange.widgets.settings import Setting, SettingProvider
from Orange.widgets.utils.signals import Input, Output
from Orange.widgets.utils.concurrent import TaskState, ConcurrentWidgetMixin
from Orange.widgets.utils.itemmodels import PyTableModel

from orangecontrib.bioinformatics.geneset import GeneSet, GeneSets
from orangecontrib.bioinformatics.widgets.utils.gui import FilterProxyModel
from orangecontrib.bioinformatics.widgets.components import GeneSetSelection
from orangecontrib.bioinformatics.widgets.utils.data import (
    TableAnnotation,
    check_table_annotation,
)
from orangecontrib.bioinformatics.widgets.utils.itemmodels import TableModel
from orangecontrib.bioinformatics.widgets.utils.itemdelegates import (
    LinkStyledItemDelegate,
)


def run(
    gene_sets: GeneSets,
    selected_gene_sets: List[Tuple[str, ...]],
    genes,
    state: TaskState,
) -> Tuple[List[int], List[GeneSet]]:
    step, steps = 0, len(gene_sets)

    filtered_gene_sets = []
    mapped_genes = []

    for gene_set in sorted(gene_sets):
        step += 1
        if step % (steps / 10) == 0:
            state.set_progress_value(100 * step / steps)

        if gene_set.hierarchy not in selected_gene_sets:
            continue

        if state.is_interruption_requested():
            return mapped_genes, filtered_gene_sets

        mapped_genes.append(gene_set.genes & genes)
        filtered_gene_sets.append(gene_set)

    return filtered_gene_sets, mapped_genes


class GeneSetsModel(TableModel):
    (
        mapped_genes,
        genes_in_set,
        collection,
        gs_name,
    ) = range(4)

    def __init__(self, gs: List[GeneSet], matching_genes: List[str]):
        def gene_set_url(row):
            return (
                gs[row].link if gs[row].link and urlparse(gs[row].link).scheme else None
            )

        columns = [
            TableModel.Column(
                {Qt.DisplayRole: 'Mapped Genes'},
                {
                    Qt.DisplayRole: lambda row: len(matching_genes[row]),
                    Qt.UserRole: lambda row: matching_genes[row],
                },
            ),
            TableModel.Column(
                {Qt.DisplayRole: 'Genes In Set'},
                {
                    Qt.DisplayRole: lambda row: len(gs[row].genes),
                    Qt.UserRole: lambda row: gs[row],
                },
            ),
            TableModel.Column(
                {Qt.DisplayRole: 'Collection'},
                {
                    Qt.DisplayRole: lambda row: ', '.join(gs[row].hierarchy),
                },
            ),
            TableModel.Column(
                {Qt.DisplayRole: 'Name'},
                {
                    Qt.DisplayRole: lambda row: gs[row].name,
                    LinkRole: gene_set_url,
                    Qt.ToolTipRole: gene_set_url,
                },
            ),
        ]

        super().__init__(len(gs), columns)

    def sortColumnData(self, column):
        return super().sortColumnData(column)


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
    use_min_count = Setting(False)

    auto_commit: bool
    auto_commit = Setting(True)

    search_pattern: str
    search_pattern = Setting('')

    # component settings
    gs_selection_component: SettingProvider = SettingProvider(GeneSetSelection)

    class Inputs:
        data = Input('Data', Table)
        custom_gene_sets = Input('Custom Gene Sets', Table)

    class Outputs:
        gene_sets = Output('Gene Sets', GeneSets)
        reduced_data = Output('Reduced Data', Table)

    class Warning(OWWidget.Warning):
        all_sets_filtered = Msg('All sets were filtered out.')

    class Error(OWWidget.Error):
        organism_mismatch = Msg(
            'Organism in input data and custom gene sets does not match'
        )
        cant_reach_host = Msg('Host orange.biolab.si is unreachable.')
        cant_load_organisms = Msg(
            'No available organisms, please check your connection.'
        )
        custom_gene_sets_table_format = Msg(
            'Custom gene sets data must have genes represented as rows.'
        )

    def __init__(self):
        super().__init__()
        # OWWidget.__init__(self)
        ConcurrentWidgetMixin.__init__(self)

        self.input_data: Optional[Table] = None

        # Control area
        box = vBox(self.controlArea, True, margin=0)
        self.gs_selection_component: GeneSetSelection = GeneSetSelection(self, box)
        self.gs_selection_component.selection_changed.connect(
            self._on_selection_changed
        )
        self.gs_selection_component.initialize()

        # Main area
        self.filter_proxy_model = FilterProxyModel()

        self.view = QTableView()
        self.view.setSelectionBehavior(QTableView.SelectRows)
        self.view.horizontalHeader().setStretchLastSection(True)
        self.view.setAlternatingRowColors(True)
        self.view.setSortingEnabled(True)
        self.view.verticalHeader().hide()
        self.view.viewport().setMouseTracking(True)
        self.view.setItemDelegateForColumn(
            GeneSetsModel.gs_name, LinkStyledItemDelegate(self.view)
        )
        self.view.setModel(self.filter_proxy_model)

        h_layout = QHBoxLayout()
        h_layout.setSpacing(100)
        h_widget = widgetBox(self.mainArea, orientation=h_layout)

        spin(
            h_widget,
            self,
            'min_count',
            0,
            1000,
            label='Min mapped genes',
            tooltip='Minimum amount of input genes that map to gene set.',
            checked='use_min_count',
            callback=self.filter_view,
            callbackOnReturn=True,
            checkCallback=self.filter_view,
        )

        self.line_edit_filter = lineEdit(h_widget, self, 'search_pattern')
        self.line_edit_filter.setPlaceholderText('Filter gene sets ...')
        self.line_edit_filter.textChanged.connect(self.filter_view)

        self.mainArea.layout().addWidget(self.view)
        self.commit_button = auto_commit(
            self.controlArea, self, 'auto_commit', '&Commit', box=False
        )

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
                str(variable.attributes.get(self.gene_location, '?'))
                for variable in self.input_data.domain.attributes
            }
        else:
            return {str(g) for g in self.input_data.get_column(self.gene_location)}

    def on_partial_result(self, _):
        pass

    def on_done(self, result: Tuple):
        self.filter_proxy_model.setSourceModel(GeneSetsModel(*result))

        if self.input_data:
            self.view.sortByColumn(GeneSetsModel.mapped_genes, Qt.DescendingOrder)
            self.view.setColumnHidden(GeneSetsModel.mapped_genes, False)
        else:
            self.view.sortByColumn(GeneSetsModel.genes_in_set, Qt.DescendingOrder)
            self.view.setColumnHidden(GeneSetsModel.mapped_genes, True)

        self.view.resizeColumnsToContents()
        # this is slow
        # self.view.resizeRowsToContents()
        self.view.selectionModel().selectionChanged.connect(self.commit.deferred)
        self.filter_view()

    def on_exception(self, ex):
        # TODO: handle possible exceptions
        raise ex

    def onDeleteWidget(self):
        self.shutdown()
        super().onDeleteWidget()

    def _on_selection_changed(self):
        if hasattr(self, 'input_genes'):
            self.start(
                run,
                self.gs_selection_component.gene_sets,
                self.gs_selection_component.selection,
                self.input_genes,
            )

    @check_table_annotation
    @Inputs.data
    def set_data(self, input_data: Table):
        self.Outputs.gene_sets.send(None)
        self.Outputs.reduced_data.send(None)
        self.input_data = None

        if input_data:
            self.input_data = input_data
            tax_id = self.input_data.attributes[TableAnnotation.tax_id]
            self.gs_selection_component.set_selected_organism_by_tax_id(tax_id)

    @Inputs.custom_gene_sets
    def handle_custom_gene_sets_input(self, custom_data):
        self.Outputs.gene_sets.send(None)
        self.Outputs.reduced_data.send(None)

        if custom_data:
            self.gs_selection_component.initialize_custom_gene_sets(custom_data)
        else:
            self.gs_selection_component.initialize_custom_gene_sets(None)

    @deferred
    def commit(self):
        selection_model = self.view.selectionModel()
        gene_sets = None
        reduced_data = None

        if selection_model:
            selection_gene_set = selection_model.selectedRows(
                GeneSetsModel.genes_in_set
            )
            gene_sets = GeneSets(
                [model_index.data(Qt.UserRole) for model_index in selection_gene_set]
            )

            selection_mapped_genes = selection_model.selectedRows(
                GeneSetsModel.mapped_genes
            )
            self.selected_rows = [
                self.filter_proxy_model.mapToSource(sel).row()
                for sel in selection_mapped_genes
            ]

            if selection_mapped_genes and self.input_genes:
                genes = [
                    model_index.data(Qt.UserRole)
                    for model_index in selection_mapped_genes
                ]
                output_genes = list(set.union(*genes))

                if self.gene_as_attr_name:
                    selected = [
                        column
                        for column in self.input_data.domain.attributes
                        if self.gene_location in column.attributes
                        and str(column.attributes[self.gene_location]) in output_genes
                    ]
                    domain = Domain(
                        selected,
                        self.input_data.domain.class_vars,
                        self.input_data.domain.metas,
                    )
                    reduced_data = self.input_data.from_table(domain, self.input_data)
                else:
                    # create filter from selected column for genes
                    only_known = table_filter.FilterStringList(
                        self.gene_location, output_genes
                    )
                    # apply filter to the data
                    reduced_data = table_filter.Values([only_known])(self.input_data)

        self.Outputs.gene_sets.send(gene_sets)
        self.Outputs.reduced_data.send(reduced_data)

    def create_filters(self):
        search_term: List[str] = self.search_pattern.lower().strip().split()
        filters = []

        if search_term:
            filters.append(
                FilterProxyModel.Filter(
                    GeneSetsModel.gs_name,
                    Qt.DisplayRole,
                    lambda value: all(fs in value.lower() for fs in search_term),
                )
            )

        if self.use_min_count:
            filters.append(
                FilterProxyModel.Filter(
                    GeneSetsModel.mapped_genes,
                    Qt.DisplayRole,
                    lambda value: value >= self.min_count,
                )
            )

        return filters

    def filter_view(self):
        filter_proxy: FilterProxyModel = self.filter_proxy_model
        model: PyTableModel = filter_proxy.sourceModel()

        # apply filtering rules
        filter_proxy.set_filters(self.create_filters())

        if model.rowCount() and not filter_proxy.rowCount():
            self.Warning.all_sets_filtered()
        else:
            self.Warning.clear()

        # this is slow
        # self.view.resizeRowsToContents()

    # def sizeHint(self):
    #     return QSize(800, 600)


if __name__ == "__main__":
    from Orange.widgets.utils.widgetpreview import WidgetPreview

    widget = WidgetPreview(OWGeneSets)
    widget.run()
