""" GeneSets """
import threading
import concurrent.futures

from functools import partial

from AnyQt.QtWidgets import (
    QTreeView, QHeaderView, QHBoxLayout
)
from AnyQt.QtCore import (
    Qt, Slot, QThread
)
from AnyQt.QtGui import (
    QColor, QStandardItemModel, QStandardItem,
)

from Orange.widgets.gui import (
    vBox, lineEdit, LinkRole, LinkStyledItemDelegate, doubleSpin,
    auto_commit, widgetLabel, spin,  widgetBox, radioButtonsInBox
)
from Orange.data import Table
from Orange.widgets.settings import Setting
from Orange.widgets.utils.signals import Output, Input

from orangecontrib.bioinformatics.widgets.utils.data import (
    TAX_ID, GENE_AS_ATTRIBUTE_NAME, GENE_ID_COLUMN, GENE_ID_ATTRIBUTE
)

from orangecontrib.bioinformatics.widgets.utils.gui import GeneSetsSelection, NumericalColumnDelegate, FilterProxyModel
from orangecontrib.bioinformatics.widgets import OWGeneSets as gene_sets
from orangecontrib.bioinformatics.utils.statistics import FDR


class OWGeneSetEnrichment(gene_sets.OWGeneSets):
    name = "Gene Set Enrichment"
    description = ""
    icon = "icons/OWGeneSets.svg"
    priority = 9

    max_p_value = Setting(0.0001)
    use_p_value = Setting(False)
    max_fdr = Setting(0.01)
    use_max_fdr = Setting(True)
    use_reference_data = Setting(True, schema_only=True)

    COUNT, REFERENCE, P_VAL, FDR, ENRICHMENT, GENES, CATEGORY, TERM = range(8)
    DATA_HEADER_LABELS = ["Count", 'Reference', 'p-Value', 'FDR', 'Enrichment', 'Genes In Set', 'Category', 'Term']

    class Inputs(gene_sets.OWGeneSets.Inputs):
        reference = Input("Reference Genes", Table)

    class Outputs:
        matched_genes = Output("Matched Genes", Table)

    def __init__(self):
        # reference data attributes
        self.reference_data = None
        self.reference_genes = None
        self.reference_tax_id = None
        self.reference_attr_names = None
        self.reference_gene_id_attribute = None
        self.reference_gene_id_column = None
        self.reference_radio_box = None

        super().__init__()


    @Inputs.reference
    def handle_reference_genes(self, data):
        """
        Set the (optional) input dataset with reference gene names.
        """

        if data:
            self.reference_data = data
            self.reference_tax_id = str(self.reference_data.attributes.get(TAX_ID, None))
            self.reference_attr_names = self.reference_data.attributes.get(GENE_AS_ATTRIBUTE_NAME, None)
            self.reference_gene_id_attribute = self.reference_data.attributes.get(GENE_ID_ATTRIBUTE, None)
            self.reference_gene_id_column = self.reference_data.attributes.get(GENE_ID_COLUMN, None)

            if not (self.reference_attr_names is not None
                    and ((self.reference_gene_id_attribute is None) ^ (self.reference_gene_id_column is None))):

                if self.reference_tax_id is None:
                    self.Error.missing_annotation()
                    return

                self.Error.missing_gene_id()
                return

            elif self.reference_tax_id is None:
                self.Error.missing_tax_id()
                return

        self.__get_reference_genes()
        self.reference_radio_box.setEnabled(bool(self.reference_data))
        self.invalidate()

    def __get_source_data(self, proxy_row_index, column):
        proxy_index = self.filter_proxy_model.index(proxy_row_index, column)
        source_index = self.filter_proxy_model.mapToSource(proxy_index)
        return source_index.data(role=Qt.DisplayRole)

    def _update_fdr(self):
        # Update the FDR in place due to a changed selected categories set and
        # results for all of these categories are already available.
        proxy = self.filter_proxy_model
        model = self.filter_proxy_model.sourceModel()

        if model is not None:
            assert isinstance(model, QStandardItemModel)

            p_values = [(i, self.__get_source_data(i, self.P_VAL)) for i in range(proxy.rowCount())]
            fdr_values = FDR([p_val for _, p_val in p_values])

            for i, fdr_val in zip([i for i, _ in p_values], fdr_values):
                proxy_index = proxy.index(i, self.FDR)
                source_index = self.filter_proxy_model.mapToSource(proxy_index)
                source_item = model.item(source_index.row(), self.FDR)
                source_item.setData(fdr_val, role=Qt.DisplayRole)
                source_item.setData(fdr_val, role=Qt.ToolTipRole)

    def __get_reference_genes(self):
        self.reference_genes = []

        if self.reference_attr_names:
            for variable in self.reference_data.domain.attributes:
                self.reference_genes.append(str(variable.attributes.get(self.reference_gene_id_attribute, '?')))
        else:
            genes, _ = self.reference_data.get_column_view(self.reference_gene_id_column)
            self.reference_genes = [str(g) for g in genes]

    def create_filters(self):
        search_term = self.search_pattern.lower().strip().split()

        # apply filtering rules
        filters = [
            FilterProxyModel.Filter(
                self.TERM, Qt.DisplayRole,
                lambda value: all(fs in value.lower() for fs in search_term))
        ]

        if self.use_min_count:
           filters.append(
               FilterProxyModel.Filter(
                   self.COUNT, Qt.DisplayRole,
                   lambda value: value >= self.min_count,
               )
           )

        if self.use_p_value:
            filters.append(
                FilterProxyModel.Filter(
                    self.P_VAL, Qt.DisplayRole,
                    lambda value: value < self.max_p_value
                )
            )

        if self.use_max_fdr:
            filters.append(
                FilterProxyModel.Filter(
                    self.FDR, Qt.DisplayRole,
                    lambda value: value < self.max_fdr
                )
            )

        return filters

    def create_partial(self):
        reference_genes = self.reference_genes if (self.use_reference_data and self.reference_data)\
                                               else self.gs_widget.gs_object.genes()

        return partial(self.set_items,
                       self.gs_widget.gs_object,
                       self.stored_gene_sets_selection,
                       set(self.input_genes),
                       self.callback,
                       reference_genes=reference_genes)

    @staticmethod
    def set_items(gene_sets, sets_to_display, genes, callback, reference_genes=None):
        model_items = []
        if not genes:
            return

        for gene_set in sorted(gene_sets):
            if gene_set.hierarchy not in sets_to_display:
                continue

            reference_genes = [] if reference_genes is None else reference_genes
            enrichemnt_result = gene_set.set_enrichment(reference_genes, genes.intersection(reference_genes))
            callback()

            if len(enrichemnt_result.query) > 0:
                category_column = QStandardItem()
                name_column = QStandardItem()
                count_column = QStandardItem()
                genes_column = QStandardItem()
                ref_column = QStandardItem()
                pval_column = QStandardItem()
                fdr_column = QStandardItem()
                enrichment_column = QStandardItem()

                category_column.setData(", ".join(gene_set.hierarchy), Qt.DisplayRole)
                name_column.setData(gene_set.name, Qt.DisplayRole)
                name_column.setData(gene_set.name, Qt.ToolTipRole)
                name_column.setData(gene_set.link, LinkRole)
                name_column.setForeground(QColor(Qt.blue))

                count_column.setData(len(enrichemnt_result.query), Qt.DisplayRole)
                count_column.setData(set(enrichemnt_result.query), Qt.UserRole)

                genes_column.setData(len(gene_set.genes), Qt.DisplayRole)
                genes_column.setData(set(gene_set.genes), Qt.UserRole)  # store genes to get then on output on selection

                ref_column.setData(len(enrichemnt_result.reference), Qt.DisplayRole)

                pval_column.setData(enrichemnt_result.p_value, Qt.DisplayRole)
                pval_column.setData(enrichemnt_result.p_value, Qt.ToolTipRole)

                enrichment_column.setData(enrichemnt_result.enrichment_score, Qt.DisplayRole)
                enrichment_column.setData(enrichemnt_result.enrichment_score, Qt.ToolTipRole)

                model_items.append([count_column, ref_column, pval_column, fdr_column, enrichment_column,
                                    genes_column, category_column, name_column])
        return model_items

    # We must extend this, because we need to update FDR values after workers finish enrichment
    @Slot(concurrent.futures.Future)
    def _init_gene_sets_finished(self, f):
        assert self.thread() is QThread.currentThread()
        assert threading.current_thread() == threading.main_thread()
        assert self._task is not None
        assert self._task.future is f
        assert f.done()

        self._task = None
        self.progress_bar.finish()
        self.setStatusMessage('')

        try:
            results = f.result()  # type: list
            [self.data_model.appendRow(model_item) for model_item in results]
            self.filter_proxy_model.setSourceModel(self.data_model)
            self.data_view.selectionModel().selectionChanged.connect(self.commit)
            self._update_fdr()
            self.filter_data_view()
            self.set_selection()
            self.update_info_box()

        except Exception as ex:
            print(ex)

    def assign_delegates(self):
        self.data_view.setItemDelegateForColumn(
            self.GENES, NumericalColumnDelegate(self)
        )

        self.data_view.setItemDelegateForColumn(
            self.COUNT, NumericalColumnDelegate(self)
        )

        self.data_view.setItemDelegateForColumn(
            self.REFERENCE, NumericalColumnDelegate(self)
        )

        self.data_view.setItemDelegateForColumn(
            self.P_VAL, NumericalColumnDelegate(self, precision=2, notation='e')
        )

        self.data_view.setItemDelegateForColumn(
            self.FDR, NumericalColumnDelegate(self, precision=2, notation='e')
        )

        self.data_view.setItemDelegateForColumn(
            self.ENRICHMENT, NumericalColumnDelegate(self, precision=1)
        )

    def setup_control_area(self):
        # Control area
        self.input_info = widgetLabel(
            widgetBox(self.controlArea, "Info", addSpace=True), 'No data on input.\n'
        )
        self.custom_gs_col_box = box = vBox(self.controlArea, 'Custom Gene Set Term Column')
        box.hide()

        self.reference_radio_box = radioButtonsInBox(
            self.controlArea, self, "use_reference_data", ["Entire genome", "Reference gene set (input)"],
            tooltips=["Use entire genome (for gene set enrichment)", "Use reference set of genes"],
            box="Reference", callback=self.invalidate)

        self.reference_radio_box.setEnabled(False)

        gene_sets_box = widgetBox(self.controlArea, "Gene Sets")
        self.gs_widget = GeneSetsSelection(gene_sets_box, self, 'stored_gene_sets_selection')
        self.gs_widget.hierarchy_tree_widget.itemClicked.connect(self.update_tree_view)

        self.commit_button = auto_commit(self.controlArea, self, "auto_commit", "&Commit", box=False)

    def setup_filter_area(self):
        h_layout = QHBoxLayout()
        h_layout.setSpacing(100)
        h_widget = widgetBox(self.mainArea, orientation=h_layout)

        spin(h_widget, self, 'min_count', 0, 100,
             label='Count',
             tooltip='Minimum genes count',
             checked='use_min_count',
             callback=self.filter_data_view,
             callbackOnReturn=True,
             checkCallback=self.filter_data_view)

        doubleSpin(h_widget, self, 'max_p_value', 0.0, 1.0, 0.0001,
                   label='p-value',
                   tooltip='Maximum p-value of the enrichment score',
                   checked='use_p_value',
                   callback=self.filter_data_view,
                   callbackOnReturn=True,
                   checkCallback=self.filter_data_view
                   )

        doubleSpin(h_widget, self, 'max_fdr', 0.0, 1.0, 0.0001,
                   label='FDR',
                   tooltip='Maximum false discovery rate',
                   checked='use_max_fdr',
                   callback=self.filter_data_view,
                   callbackOnReturn=True,
                   checkCallback=self.filter_data_view
                   )

        self.line_edit_filter = lineEdit(h_widget, self, 'search_pattern')
        self.line_edit_filter.setPlaceholderText('Filter gene sets ...')
        self.line_edit_filter.textChanged.connect(self.filter_data_view)

    def setup_gui(self):
        # control area
        self.setup_control_area()

        # main area
        self.data_view = QTreeView()
        self.setup_filter_model()
        self.setup_filter_area()
        self.data_view.setAlternatingRowColors(True)
        self.data_view.sortByColumn(self.COUNT, Qt.DescendingOrder)
        self.data_view.setSortingEnabled(True)
        self.data_view.setSelectionMode(QTreeView.ExtendedSelection)
        self.data_view.setEditTriggers(QTreeView.NoEditTriggers)
        self.data_view.viewport().setMouseTracking(False)
        self.data_view.setItemDelegateForColumn(self.TERM, LinkStyledItemDelegate(self.data_view))

        self.mainArea.layout().addWidget(self.data_view)

        self.data_view.header().setSectionResizeMode(QHeaderView.ResizeToContents)
        self.assign_delegates()


if __name__ == "__main__":
    from AnyQt.QtWidgets import QApplication
    app = QApplication([])
    ow = OWGeneSetEnrichment()
    ow.show()
    app.exec_()
