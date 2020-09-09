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
from Orange.widgets.gui import (
    LinkRole,
    LinkStyledItemDelegate,
    spin,
    vBox,
    lineEdit,
    widgetBox,
    doubleSpin,
    auto_commit,
    radioButtonsInBox,
)
from Orange.widgets.widget import Msg, OWWidget
from Orange.widgets.settings import Setting, SettingProvider
from Orange.widgets.utils.signals import Input, Output
from Orange.widgets.utils.concurrent import TaskState, ConcurrentWidgetMixin

from orangecontrib.bioinformatics.geneset import GeneSets
from orangecontrib.bioinformatics.ncbi.gene import GeneInfo
from orangecontrib.bioinformatics.utils.statistics import FDR
from orangecontrib.bioinformatics.widgets.utils.gui import FilterProxyModel, NumericalColumnDelegate
from orangecontrib.bioinformatics.widgets.components import GeneSetSelection
from orangecontrib.bioinformatics.widgets.utils.data import TableAnnotation, check_table_annotation


class Results(SimpleNamespace):
    items: List[QStandardItem] = []


def run(
    gene_sets: GeneSets, selected_gene_sets: List[Tuple[str, ...]], genes, state: TaskState, reference_genes=None
) -> Results:
    results = Results()
    items = []
    step, steps = 0, len(gene_sets)

    def set_progress():
        nonlocal step
        step += 1
        state.set_progress_value(100 * (step / steps))

    if not genes:
        return results

    state.set_status('Calculating...')

    for gene_set in sorted(gene_sets):
        set_progress()

        if gene_set.hierarchy not in selected_gene_sets:
            continue

        if state.is_interruption_requested():
            return results

        reference_genes = [] if reference_genes is None else reference_genes
        enrichemnt_result = gene_set.set_enrichment(reference_genes, genes.intersection(reference_genes))

        if len(enrichemnt_result.query) > 0:
            category_column = QStandardItem()
            term_column = QStandardItem()
            count_column = QStandardItem()
            genes_column = QStandardItem()
            ref_column = QStandardItem()
            pval_column = QStandardItem()
            fdr_column = QStandardItem()
            enrichment_column = QStandardItem()

            category_column.setData(", ".join(gene_set.hierarchy), Qt.DisplayRole)
            term_column.setData(gene_set.name, Qt.DisplayRole)
            term_column.setData(gene_set.name, Qt.ToolTipRole)
            # there was some cases when link string was not empty string but not valid (e.g. "_")
            if gene_set.link and urlparse(gene_set.link).scheme:
                term_column.setData(gene_set.link, LinkRole)
                term_column.setForeground(QColor(Qt.blue))

            count_column.setData(len(enrichemnt_result.query), Qt.DisplayRole)
            count_column.setData(set(enrichemnt_result.query), Qt.UserRole)

            genes_column.setData(len(gene_set.genes), Qt.DisplayRole)
            genes_column.setData(set(gene_set.genes), Qt.UserRole)  # store genes to get then on output on selection

            ref_column.setData(len(enrichemnt_result.reference), Qt.DisplayRole)

            pval_column.setData(enrichemnt_result.p_value, Qt.DisplayRole)
            pval_column.setData(enrichemnt_result.p_value, Qt.ToolTipRole)

            enrichment_column.setData(enrichemnt_result.enrichment_score, Qt.DisplayRole)
            enrichment_column.setData(enrichemnt_result.enrichment_score, Qt.ToolTipRole)

            items.append(
                [
                    count_column,
                    ref_column,
                    pval_column,
                    fdr_column,
                    enrichment_column,
                    genes_column,
                    category_column,
                    term_column,
                ]
            )

    results.items = items
    return results


class Header(IntEnum):
    count = 0
    reference = 1
    p_val = 2
    fdr = 3
    enrichment = 4
    genes = 5
    category = 6
    term = 7

    @staticmethod
    def labels():
        return ['Count', 'Reference', 'p-Value', 'FDR', 'Enrichment', 'Genes In Set', 'Category', 'Term']


class OWGeneSets(OWWidget, ConcurrentWidgetMixin):
    name = "Gene Set Enrichment"
    description = ""
    icon = "icons/OWGeneSets.svg"
    priority = 90
    want_main_area = True

    organism = Setting(None, schema_only=True)
    stored_gene_sets_selection = Setting([], schema_only=True)
    selected_rows = Setting([], schema_only=True)

    min_count: int
    min_count = Setting(5, schema_only=True)

    use_min_count: bool
    use_min_count = Setting(True, schema_only=True)

    max_p_value: int
    max_p_value = Setting(0.0001, schema_only=True)

    use_p_value: bool
    use_p_value = Setting(False, schema_only=True)

    max_fdr: int
    max_fdr = Setting(0.01, schema_only=True)

    use_max_fdr: bool
    use_max_fdr = Setting(True, schema_only=True)

    search_pattern: str
    search_pattern = Setting('')

    use_reference_data = Setting(False, schema_only=True)

    auto_commit: bool
    auto_commit = Setting(False)

    # component settings
    gs_selection_component: SettingProvider = SettingProvider(GeneSetSelection)

    class Inputs:
        data = Input('Data', Table)
        custom_gene_sets = Input('Custom Gene Sets', Table)
        reference = Input('Reference Genes', Table)

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

        self.reference_radio_box = radioButtonsInBox(
            self.controlArea,
            self,
            'use_reference_data',
            ['Entire genome', 'Reference gene set (input)'],
            tooltips=['Use entire genome (for gene set enrichment)', 'Use reference set of genes'],
            box='Reference',
            callback=self._on_selection_changed,
        )
        self.reference_radio_box.setEnabled(False)
        self.reference_genes: Optional[List[str]] = None

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
        self.tree_view.setItemDelegateForColumn(Header.p_val, NumericalColumnDelegate(self, precision=2, notation='e'))
        self.tree_view.setItemDelegateForColumn(Header.fdr, NumericalColumnDelegate(self, precision=2, notation='e'))
        self.tree_view.setItemDelegateForColumn(Header.enrichment, NumericalColumnDelegate(self, precision=1))
        self.tree_view.setItemDelegateForColumn(Header.reference, NumericalColumnDelegate(self))
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

        doubleSpin(
            h_widget,
            self,
            'max_p_value',
            0.0,
            1.0,
            0.0001,
            label='p-value',
            tooltip='Maximum p-value of the enrichment score',
            checked='use_p_value',
            callback=self.filter_view,
            callbackOnReturn=True,
            checkCallback=self.filter_view,
        )

        doubleSpin(
            h_widget,
            self,
            'max_fdr',
            0.0,
            1.0,
            0.0001,
            label='FDR',
            tooltip='Maximum false discovery rate',
            checked='use_max_fdr',
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
        self.gene_info: Optional[GeneInfo] = None

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

    def get_data_from_source_model(self, proxy_row_index, column):
        proxy_index = self.filter_proxy_model.index(proxy_row_index, column)
        source_index = self.filter_proxy_model.mapToSource(proxy_index)
        return source_index.data(role=Qt.DisplayRole)

    def _update_fdr(self):
        proxy = self.filter_proxy_model
        model = self.filter_proxy_model.sourceModel()

        if model is not None:
            assert isinstance(model, QStandardItemModel)
            p_values = [(i, self.get_data_from_source_model(i, Header.p_val)) for i in range(proxy.rowCount())]
            fdr_values = FDR([p_val for _, p_val in p_values])

            for i, fdr_val in zip([i for i, _ in p_values], fdr_values):
                proxy_index = proxy.index(i, Header.fdr)
                source_index = self.filter_proxy_model.mapToSource(proxy_index)
                source_item = model.item(source_index.row(), Header.fdr)
                source_item.setData(fdr_val, role=Qt.DisplayRole)
                source_item.setData(fdr_val, role=Qt.ToolTipRole)

    def on_partial_result(self, _):
        pass

    def on_done(self, result: Results):
        model = QStandardItemModel()
        for item in result.items:
            model.appendRow(item)

        model.setSortRole(Qt.UserRole)
        model.setHorizontalHeaderLabels(Header.labels())

        self.filter_proxy_model.setSourceModel(model)
        self.filter_proxy_model.reset_filters()
        self._update_fdr()
        self.filter_view()
        self.update_info_box()
        self.tree_view.selectionModel().selectionChanged.connect(self.commit)

    def on_exception(self, ex):
        # TODO: handle possible exceptions
        raise ex

    def onDeleteWidget(self):
        self.shutdown()
        super().onDeleteWidget()

    def _on_selection_changed(self):
        if self.use_reference_data and self.reference_genes:
            ref_genes = self.reference_genes
        else:
            ref_genes = self.gene_info.keys() if self.gene_info else []

        self.start(
            run,
            self.gs_selection_component.gene_sets,
            self.gs_selection_component.selection,
            self.input_genes,
            reference_genes=ref_genes,
        )

    @Inputs.data
    @check_table_annotation
    def set_data(self, input_data: Table):
        self.Outputs.matched_genes.send(None)
        self.input_data = None
        self.num_of_selected_genes = 0

        if input_data:
            self.input_data = input_data
            self.gene_info = GeneInfo(self.tax_id)
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

    @Inputs.reference
    @check_table_annotation
    def set_reference_data(self, input_data: Table):
        self.reference_genes = []
        self.use_reference_data = bool(input_data)
        self.reference_radio_box.setEnabled(bool(input_data))

        if input_data:
            attrs = input_data.attributes
            if attrs[TableAnnotation.gene_as_attr_name]:
                for variable in input_data.domain.attributes:
                    self.reference_genes.append(
                        str(variable.attributes.get(attrs[TableAnnotation.gene_id_attribute], '?'))
                    )
            else:
                genes, _ = self.reference_data.get_column_view(attrs[TableAnnotation.gene_id_column])
                self.reference_genes = [str(g) for g in genes]

        self._on_selection_changed()

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

        if self.use_p_value:
            filters.append(
                FilterProxyModel.Filter(Header.p_val, Qt.DisplayRole, lambda value: value < self.max_p_value)
            )

        if self.use_max_fdr:
            filters.append(FilterProxyModel.Filter(Header.fdr, Qt.DisplayRole, lambda value: value < self.max_fdr))

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
