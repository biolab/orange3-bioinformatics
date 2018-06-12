""" Gene Sets Enrichment"""
from AnyQt.QtWidgets import (
    QTreeView, QTreeWidget
)
from AnyQt.QtCore import (
    Qt
)
from AnyQt.QtGui import (
    QColor, QStandardItemModel, QStandardItem
)

from Orange.data import Table
from Orange.widgets.utils.signals import Input
from Orange.widgets.settings import Setting
from Orange.widgets.widget import OWWidget, Msg

from Orange.widgets.gui import (
    vBox, lineEdit, LinkRole,
    auto_commit, widgetLabel, spin, comboBox, widgetBox, QHBoxLayout, doubleSpin, radioButtonsInBox
)

from orangecontrib.bioinformatics.widgets.utils.gui import NumericalColumnDelegate, FilterProxyModel
from orangecontrib.bioinformatics.widgets import OWGeneSets
from orangecontrib.bioinformatics.utils.statistics import FDR


from orangecontrib.bioinformatics.widgets.utils.data import (
    TAX_ID, GENE_AS_ATTRIBUTE_NAME, GENE_ID_COLUMN, GENE_ID_ATTRIBUTE
)


class OWGeneSetsEnrichment(OWGeneSets.OWGeneSets):
    name = "Set Enrichment"
    description = "This is widget description"
    icon = "icons/OWSetEnrichment.svg"
    priority = 10

    COUNT, REFERENCE, P_VAL, FDR, ENRICHMENT, GENES, CATEGORY, TERM = range(8)
    DATA_HEADER_LABELS = ["Count", 'Reference', 'p-Value', 'FDR', 'Enrichment', 'Genes In Set', 'Category', 'Term']

    # min_count = Setting(5)
    # use_min_count = Setting(True)

    max_p_value = Setting(0.0001)
    use_p_value = Setting(False)

    max_fdr = Setting(0.01)
    use_max_fdr = Setting(True)

    use_reference_data = Setting(True)

    class Inputs(OWGeneSets.OWGeneSets.Inputs):
        reference = Input("Reference genes", Table)

    class Warning(OWWidget.Warning):
        all_sets_filtered = Msg('All sets were filtered out.')

    def __init__(self):
        self.line_edit_filter = None
        self.reference_radio_box = None
        self.reference_data = None
        self.reference_genes = None

        self.reference_tax_id = None
        self.reference_attr_names = None
        self.reference_gene_id_attribute = None
        self.reference_gene_id_column = None

        super().__init__()

    def __get_reference_genes(self):
        self.reference_genes = []

        if self.reference_attr_names:
            for variable in self.reference_data.domain.attributes:
                self.reference_genes.append(str(variable.attributes.get(self.reference_gene_id_attribute, '?')))
        else:
            genes, _ = self.reference_data.get_column_view(self.reference_gene_id_column)
            self.reference_genes = [str(g) for g in genes]

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

    def update_tree_view(self):
        if self.use_reference_data and self.reference_data:
            self.init_gene_sets(reference_genes=self.reference_genes)
        else:
            self.init_gene_sets()

    def invalidate(self):
        super().invalidate()

    def _init_gene_sets_finished(self, *args, **kwargs):
        super()._init_gene_sets_finished(*args, **kwargs)
        self.filter_data_view()
        self._update_fdr()

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
            self.P_VAL, NumericalColumnDelegate(self, precision=4)
        )

        self.data_view.setItemDelegateForColumn(
            self.FDR, NumericalColumnDelegate(self, precision=4)
        )

        self.data_view.setItemDelegateForColumn(
            self.ENRICHMENT, NumericalColumnDelegate(self, precision=1)
        )

    def filter_data_view(self):
        filter_proxy = self.filter_proxy_model  # type: FilterProxyModel
        model = filter_proxy.sourceModel()      # type: QStandardItemModel

        assert isinstance(model, QStandardItemModel)

        # update fdr values
        self._update_fdr()

        search_term = self.search_pattern.lower().strip().split()

        # apply filtering rules
        filters = [
            FilterProxyModel.Filter(
                self.TERM, Qt.DisplayRole,
                lambda value: all(fs in value.lower() for fs in search_term))
        ]

        # if self.use_min_count:
        #    filters.append(
        #        FilterProxyModel.Filter(
        #            self.COUNT, Qt.DisplayRole,
        #            lambda value: value >= self.min_count,
        #        )
        #    )

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

        filter_proxy.set_filters(filters)

        if model.rowCount() and not filter_proxy.rowCount():
            self.Warning.all_sets_filtered()
        else:
            self.Warning.clear()

    def setup_filter_model(self):
        self.filter_proxy_model = FilterProxyModel()
        self.filter_proxy_model.setFilterKeyColumn(self.TERM)

    def setup_filter_area(self):
        h_layout = QHBoxLayout()
        h_layout.setSpacing(100)
        h_widget = widgetBox(self.mainArea, orientation=h_layout)

        # spin(h_widget, self, 'min_count', 0, 100,
        #    label='Count',
        #    tooltip='Minimum genes count',
        #     checked='use_min_count',
        #     callback=self.filter_data_view,
        #     callbackOnReturn=True,
        #     checkCallback=self.filter_data_view)

        doubleSpin(h_widget, self, 'max_p_value', 0.0, 1.0, 0.0001,
                   label='p-value',
                   tooltip='Maximum p-value',
                   checked='use_p_value',
                   callback=self.filter_data_view,
                   callbackOnReturn=True,
                   checkCallback=self.filter_data_view
                   )

        doubleSpin(h_widget, self, 'max_fdr', 0.0, 1.0, 0.0001,
                   label='FDR',
                   tooltip='Maximum False discovery rate',
                   checked='use_max_fdr',
                   callback=self.filter_data_view,
                   callbackOnReturn=True,
                   checkCallback=self.filter_data_view
                   )

        self.line_edit_filter = lineEdit(h_widget, self, 'search_pattern')
        self.line_edit_filter.setPlaceholderText('Filter gene sets ...')
        self.line_edit_filter.textChanged.connect(self.filter_data_view)

    def setup_control_area(self):
        info_box = vBox(self.controlArea, 'Info')
        self.input_info = widgetLabel(info_box)

        box = vBox(self.controlArea, "Minimum Count")
        self.spin_widget = spin(box, self, 'matched_treshold', 0, 1000,
                                callback=self.invalidate)

        box = vBox(self.controlArea, "Custom Gene Sets")
        self.gs_label_combobox = comboBox(box, self, "gene_set_label", sendSelectedValue=True,
                                          model=self.feature_model, callback=self.invalidate)
        self.gs_label_combobox.setDisabled(True)

        self.reference_radio_box = radioButtonsInBox(
            self.controlArea, self, "use_reference_data", ["Entire genome", "Reference gene set (input)"],
            tooltips=["Use entire genome (for gene set enrichment)", "Use reference set of genes"],
            box="Reference", callback=self.invalidate)

        self.reference_radio_box.setEnabled(False)

        hierarchy_box = widgetBox(self.controlArea, "Gene Set Databases")
        self.hierarchy_widget = QTreeWidget(self)
        self.hierarchy_widget.setEditTriggers(QTreeView.NoEditTriggers)
        self.hierarchy_widget.setHeaderLabels(self.HIERARCHY_HEADER_LABELS)
        self.hierarchy_widget.itemClicked.connect(self.update_tree_view)
        hierarchy_box.layout().addWidget(self.hierarchy_widget)

        self.commit_button = auto_commit(self.controlArea, self, "auto_commit", "&Commit", box=False)

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

    @staticmethod
    def set_items(gene_sets, sets_to_display, genes, matched_treshold, ref, callback):
        model_items = []
        if not genes:
            return

        for gene_set in gene_sets:
            if gene_set.hierarchy not in sets_to_display:
                continue
            enrichemnt_result = gene_set.set_enrichment(ref, genes.intersection(ref))
            callback()

            if len(enrichemnt_result.query) >= matched_treshold:
                category_column = QStandardItem()
                name_column = QStandardItem()
                count_column = QStandardItem()
                genes_column = QStandardItem()
                ref_column = QStandardItem()
                pval_column = QStandardItem()
                fdr_column = QStandardItem()
                enrichemnt_column = QStandardItem()

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

                enrichemnt_column.setData(enrichemnt_result.enrichment_score, Qt.DisplayRole)
                enrichemnt_column.setData(enrichemnt_result.enrichment_score, Qt.ToolTipRole)

                model_items.append([count_column, ref_column, pval_column, fdr_column, enrichemnt_column,
                                    genes_column, category_column, name_column])
        return model_items


if __name__ == "__main__":
    from AnyQt.QtWidgets import QApplication
    app = QApplication([])
    ow = OWGeneSetsEnrichment()
    ow.show()
    app.exec_()
