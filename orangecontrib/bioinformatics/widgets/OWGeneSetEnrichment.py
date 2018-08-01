""" GeneSets """
import threading
import concurrent.futures

from functools import partial
from collections import defaultdict

from AnyQt.QtWidgets import (
    QTreeView, QHeaderView, QHBoxLayout
)
from AnyQt.QtCore import (
    Qt, QSize, QThreadPool, Slot, QThread
)
from AnyQt.QtGui import (
    QColor, QStandardItemModel, QStandardItem,
)

from Orange.widgets.gui import (
    vBox, lineEdit, ProgressBar, LinkRole, LinkStyledItemDelegate, doubleSpin,
    auto_commit, widgetLabel, spin, comboBox, widgetBox, radioButtonsInBox
)
from Orange.data import Domain, Table, DiscreteVariable, StringVariable
from Orange.widgets.widget import OWWidget, Msg
from Orange.widgets.settings import Setting, ContextSetting
from Orange.widgets.utils.signals import Output, Input
from Orange.widgets.utils.itemmodels import DomainModel
from Orange.widgets.utils.concurrent import ThreadExecutor, FutureWatcher, methodinvoke

from orangecontrib.bioinformatics.widgets.utils.data import (
    TAX_ID, GENE_AS_ATTRIBUTE_NAME, GENE_ID_COLUMN, GENE_ID_ATTRIBUTE,
    ERROR_ON_MISSING_ANNOTATION, ERROR_ON_MISSING_GENE_ID, ERROR_ON_MISSING_TAX_ID
)

from orangecontrib.bioinformatics.widgets.utils.gui import GeneSetsSelection
from orangecontrib.bioinformatics.widgets.utils.gui import NumericalColumnDelegate, FilterProxyModel
from orangecontrib.bioinformatics.widgets.utils.settings import OrganismContextHandler
from orangecontrib.bioinformatics.utils.statistics import FDR
from orangecontrib.bioinformatics import geneset


class Task:
    future = None
    watcher = None
    cancelled = False

    def cancel(self):
        self.cancelled = True
        self.future.cancel()
        concurrent.futures.wait([self.future])


class OWGeneSets(OWWidget):
    name = "Gene Set Enrichment"
    description = ""
    icon = "icons/OWGeneSets.svg"
    priority = 9
    want_main_area = True
    settingsHandler = OrganismContextHandler()

    # settings
    auto_commit = Setting(True)
    stored_gene_sets_selection = ContextSetting([])
    organism = ContextSetting(None)

    min_count = Setting(5)
    use_min_count = Setting(True)

    max_p_value = Setting(0.0001)
    use_p_value = Setting(False)

    max_fdr = Setting(0.01)
    use_max_fdr = Setting(True)

    use_reference_data = Setting(True)

    custom_gene_set_indicator = Setting(None)

    COUNT, REFERENCE, P_VAL, FDR, ENRICHMENT, GENES, CATEGORY, TERM = range(8)
    DATA_HEADER_LABELS = ["Count", 'Reference', 'p-Value', 'FDR', 'Enrichment', 'Genes In Set', 'Category', 'Term']

    class Inputs:
        genes = Input("Data", Table)
        custom_sets = Input('Custom Gene Sets', Table)
        reference = Input("Reference Genes", Table)

    class Outputs:
        matched_genes = Output("Matched Genes", Table)

    class Information(OWWidget.Information):
        pass

    class Warning(OWWidget.Warning):
        all_sets_filtered = Msg('All sets were filtered out.')

    class Error(OWWidget.Error):
        organism_mismatch = Msg('Organism in input data and custom gene sets does not match')
        missing_annotation = Msg(ERROR_ON_MISSING_ANNOTATION)
        missing_gene_id = Msg(ERROR_ON_MISSING_GENE_ID)
        missing_tax_id = Msg(ERROR_ON_MISSING_TAX_ID)
        cant_reach_host = Msg("Host orange.biolab.si is unreachable.")
        cant_load_organisms = Msg("No available organisms, please check your connection.")

    def __init__(self):
        super().__init__()

        # commit
        self.commit_button = None

        # progress bar
        self.progress_bar = None
        self.progress_bar_iterations = None

        # data
        self.input_data = None
        self.input_genes = []
        self.tax_id = None
        self.use_attr_names = None
        self.gene_id_attribute = None
        self.gene_id_column = None

        # custom gene sets
        self.custom_data = None
        self.feature_model = DomainModel(valid_types=(DiscreteVariable, StringVariable))
        self.gs_label_combobox = None
        self.custom_tax_id = None
        self.custom_use_attr_names = None
        self.custom_gene_id_attribute = None
        self.custom_gene_id_column = None
        self.num_of_custom_sets = None

        # Gene Sets widget
        self.gs_widget = None

        # reference genes
        self.reference_radio_box = None
        self.reference_data = None
        self.reference_genes = None

        self.reference_tax_id = None
        self.reference_attr_names = None
        self.reference_gene_id_attribute = None
        self.reference_gene_id_column = None

        # info box
        self.input_info = None
        self.num_of_sel_genes = 0

        # filter
        self.line_edit_filter = None
        self.search_pattern = ''
        self.organism_select_combobox = None

        # data model view
        self.data_view = None
        self.data_model = None

        # gene matcher NCBI
        self.gene_matcher = None

        # filter proxy model
        self.filter_proxy_model = None

        # hierarchy widget
        self.hierarchy_widget = None
        self.hierarchy_state = None

        # spinbox
        self.spin_widget = None

        # threads
        self.threadpool = QThreadPool(self)
        self.workers = None

        self._task = None  # type: Optional[Task]
        self._executor = ThreadExecutor()

        # gui
        self.setup_gui()

    def __reset_widget_state(self):
        self.update_info_box()
        # clear data view
        self.init_item_model()
        # reset filters
        self.setup_filter_model()

    def cancel(self):
        """
        Cancel the current task (if any).
        """
        if self._task is not None:
            self._task.cancel()
            assert self._task.future.done()
            # disconnect the `_task_finished` slot
            self._task.watcher.done.disconnect(self._init_gene_sets_finished)
            self._task = None

    @Slot()
    def progress_advance(self):
        # GUI should be updated in main thread. That's why we are calling advance method here
        if self.progress_bar:
            self.progress_bar.advance()

    def __get_input_genes(self):
        self.input_genes = []

        if self.use_attr_names:
            for variable in self.input_data.domain.attributes:
                self.input_genes.append(str(variable.attributes.get(self.gene_id_attribute, '?')))
        else:
            genes, _ = self.input_data.get_column_view(self.gene_id_column)
            self.input_genes = [str(g) for g in genes]

    def handle_custom_gene_sets(self, select_customs_flag=False):
        if self.custom_gene_set_indicator:
            if self.custom_data is not None and self.custom_gene_id_column is not None:

                if self.__check_organism_mismatch():
                    self.gs_label_combobox.setDisabled(True)
                    self.Error.organism_mismatch()
                    self.gs_widget.update_gs_hierarchy()
                    return

                if isinstance(self.custom_gene_set_indicator, DiscreteVariable):
                    gene_sets_names = self.custom_gene_set_indicator.values
                else:
                    gene_sets_names, _ = self.custom_data.get_column_view(self.custom_gene_set_indicator)

                self.num_of_custom_sets = len(set(gene_sets_names))
                gene_names, _ = self.custom_data.get_column_view(self.custom_gene_id_column)
                hierarchy_title = (self.custom_data.name if self.custom_data.name else 'Custom sets', )
                try:
                    self.gs_widget.add_custom_sets(
                        gene_sets_names, gene_names,
                        hierarchy_title=hierarchy_title, select_customs_flag=select_customs_flag)
                except geneset.GeneSetException:
                    pass
                self.gs_label_combobox.setDisabled(False)
            else:
                self.gs_widget.update_gs_hierarchy()

        self.update_info_box()

    def update_tree_view(self):
        if self.use_reference_data and self.reference_data:
            self.init_gene_sets(reference_genes=self.reference_genes)
        else:
            self.init_gene_sets()

    def invalidate(self):
        # clear
        self.__reset_widget_state()
        self.update_info_box()

        if self.input_data is not None:
            # setup
            self.__get_input_genes()
            self.update_tree_view()

    def __check_organism_mismatch(self):
        """ Check if organisms from different inputs match.

        :return: True if there is a mismatch
        """
        if self.tax_id is not None and self.custom_tax_id is not None:
            return self.tax_id != self.custom_tax_id
        return False

    def __get_reference_genes(self):
        self.reference_genes = []

        if self.reference_attr_names:
            for variable in self.reference_data.domain.attributes:
                self.reference_genes.append(str(variable.attributes.get(self.reference_gene_id_attribute, '?')))
        else:
            genes, _ = self.reference_data.get_column_view(self.reference_gene_id_column)
            self.reference_genes = [str(g) for g in genes]

    def _handle_future_model(self):
        if self.custom_gene_set_indicator in self.feature_model:
            index = self.feature_model.indexOf(self.custom_gene_set_indicator)
            self.custom_gene_set_indicator = self.feature_model[index]
        else:
            if self.feature_model:
                self.custom_gene_set_indicator = self.feature_model[0]
            else:
                self.custom_gene_set_indicator = None

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

    @Inputs.custom_sets
    def handle_custom_input(self, data):
        self.Error.clear()
        self.__reset_widget_state()
        self.custom_data = None
        self.custom_tax_id = None
        self.custom_use_attr_names = None
        self.custom_gene_id_attribute = None
        self.custom_gene_id_column = None
        self.gs_label_combobox.setDisabled(True)
        self.feature_model.set_domain(None)

        if data:
            self.custom_data = data
            self.feature_model.set_domain(self.custom_data.domain)
            self.custom_tax_id = str(self.custom_data.attributes.get(TAX_ID, None))
            self.custom_use_attr_names = self.custom_data.attributes.get(GENE_AS_ATTRIBUTE_NAME, None)
            self.custom_gene_id_attribute = self.custom_data.attributes.get(GENE_ID_ATTRIBUTE, None)
            self.custom_gene_id_column = self.custom_data.attributes.get(GENE_ID_COLUMN, None)
            self._handle_future_model()

        self.gs_label_combobox.setDisabled(True)
        self.gs_widget.clear_custom_sets()
        self.handle_custom_gene_sets(select_customs_flag=True)
        self.invalidate()

    @Inputs.genes
    def handle_genes_input(self, data):
        self.closeContext()
        self.Error.clear()
        self.__reset_widget_state()
        # clear output
        self.Outputs.matched_genes.send(None)
        # clear input values
        self.input_genes = []
        self.input_data = None
        self.stored_gene_sets_selection = tuple()
        self.tax_id = None
        self.use_attr_names = None
        self.gene_id_attribute = None
        self.gs_label_combobox.setDisabled(True)
        self.gs_widget.clear()
        self.gs_widget.clear_gene_sets()
        self.update_info_box()

        if data:
            self.input_data = data
            self.tax_id = str(self.input_data.attributes.get(TAX_ID, None))
            self.use_attr_names = self.input_data.attributes.get(GENE_AS_ATTRIBUTE_NAME, None)
            self.gene_id_attribute = self.input_data.attributes.get(GENE_ID_ATTRIBUTE, None)
            self.gene_id_column = self.input_data.attributes.get(GENE_ID_COLUMN, None)
            self.update_info_box()

            if not(self.use_attr_names is not None
                   and ((self.gene_id_attribute is None) ^ (self.gene_id_column is None))):

                if self.tax_id is None:
                    self.Error.missing_annotation()
                    return

                self.Error.missing_gene_id()
                return

            elif self.tax_id is None:
                self.Error.missing_tax_id()
                return

            if self.__check_organism_mismatch():
                self.Error.organism_mismatch()
                return

            self.openContext(self.tax_id)

            self.gs_widget.load_gene_sets(self.tax_id)

            # if input data change, we need to set feature model again
            if self.custom_data:
                self.gs_widget.clear_custom_sets()
                self._handle_future_model()
                self.handle_custom_gene_sets()

            self.invalidate()

    def update_info_box(self):
        info_string = ''
        if self.input_genes:
            info_string += '{} unique gene names on input.\n'.format(len(self.input_genes))
            info_string += '{} genes on output.\n'.format(self.num_of_sel_genes)
        else:
            if self.input_data:
                if not any([self.gene_id_column, self.gene_id_attribute]):
                    info_string += 'Input data with incorrect meta data.\nUse Gene Name Matcher widget.'
            else:
                info_string += 'No data on input.\n'

        if self.custom_data:
            info_string += '{} marker genes in {} sets\n'.format(self.custom_data.X.shape[0], self.num_of_custom_sets)


        self.input_info.setText(info_string)

    def init_gene_sets(self, reference_genes=None):
        if self._task is not None:
            self.cancel()
        assert self._task is None

        self._task = Task()
        progress_advance = methodinvoke(self, "progress_advance")

        def callback():
            if self._task.cancelled:
                raise KeyboardInterrupt()
            if self.progress_bar:
                progress_advance()

        if reference_genes is None:
            reference_genes = self.gs_widget.gs_object.genes()

        self.init_item_model()

        sets_to_display = self.gs_widget.get_hierarchies(only_selected=True)
        # save setting on selected hierarchies
        self.stored_gene_sets_selection = sets_to_display
        # save context
        self.closeContext()

        f = partial(self.set_items,
                    self.gs_widget.gs_object,
                    sets_to_display,
                    set(self.input_genes),
                    reference_genes,
                    self.min_count if self.use_min_count else 1,
                    callback=callback)

        progress_iterations = sum([len(g_set) for hier, g_set
                                   in self.gs_widget.gs_object.map_hierarchy_to_sets().items()
                                   if hier in sets_to_display])

        self.progress_bar = ProgressBar(self, iterations=progress_iterations)

        self._task.future = self._executor.submit(f)

        self._task.watcher = FutureWatcher(self._task.future)
        self._task.watcher.done.connect(self._init_gene_sets_finished)

        self.openContext(self.tax_id)

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
            self.update_info_box()
        except Exception as ex:
            print(ex)

    def filter_data_view(self):
        filter_proxy = self.filter_proxy_model  # type: FilterProxyModel
        model = filter_proxy.sourceModel()      # type: QStandardItemModel

        if isinstance(model, QStandardItemModel):

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

    def commit(self):
        selection_model = self.data_view.selectionModel()

        if selection_model:
            # genes_from_set = selection_model.selectedRows(GENES)
            matched_genes = selection_model.selectedRows(self.COUNT)

            if matched_genes and self.input_genes:
                genes = [model_index.data(Qt.UserRole) for model_index in matched_genes]
                output_genes = [gene_name for gene_name in list(set.union(*genes))]
                self.num_of_sel_genes = len(output_genes)
                self.update_info_box()

                if self.use_attr_names:
                    selected = [column for column in self.input_data.domain.attributes
                                if self.gene_id_attribute in column.attributes and
                                str(column.attributes[self.gene_id_attribute]) in output_genes]

                    domain = Domain(selected, self.input_data.domain.class_vars, self.input_data.domain.metas)
                    new_data = self.input_data.from_table(domain, self.input_data)
                    self.Outputs.matched_genes.send(new_data)

                else:
                    selected_rows = []
                    for row_index, row in enumerate(self.input_data):
                        gene_in_row = str(row[self.gene_id_column])
                        if gene_in_row in self.input_genes and gene_in_row in output_genes:
                            selected_rows.append(row_index)

                    if selected_rows:
                        selected = self.input_data[selected_rows]
                    else:
                        selected = None

                    self.Outputs.matched_genes.send(selected)

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

    def setup_filter_model(self):
        self.filter_proxy_model = FilterProxyModel()
        self.filter_proxy_model.setFilterKeyColumn(self.TERM)
        self.data_view.setModel(self.filter_proxy_model)

    def setup_filter_area(self):
        h_layout = QHBoxLayout()
        h_layout.setSpacing(100)
        h_widget = widgetBox(self.mainArea, orientation=h_layout)

        spin(h_widget, self, 'min_count', 0, 100,
             label='Count',
             tooltip='Minimum genes count',
             checked='use_min_count',
             callback=self.invalidate,
             callbackOnReturn=True,
             checkCallback=self.invalidate)

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

    def on_gene_set_indicator_changed(self):
        self._handle_future_model()
        self.gs_widget.clear_custom_sets()
        self.handle_custom_gene_sets()
        self.invalidate()

    def setup_control_area(self):
        # Control area
        self.input_info = widgetLabel(
            widgetBox(self.controlArea, "Info", addSpace=True), 'No data on input.\n'
        )

        box = vBox(self.controlArea, 'Custom Gene Set Term Column')
        self.gs_label_combobox = comboBox(box, self, "custom_gene_set_indicator", sendSelectedValue=True,
                                          model=self.feature_model, callback=self.on_gene_set_indicator_changed)
        self.gs_label_combobox.setDisabled(True)

        self.reference_radio_box = radioButtonsInBox(
            self.controlArea, self, "use_reference_data", ["Entire genome", "Reference gene set (input)"],
            tooltips=["Use entire genome (for gene set enrichment)", "Use reference set of genes"],
            box="Reference", callback=self.invalidate)

        self.reference_radio_box.setEnabled(False)

        gene_sets_box = widgetBox(self.controlArea, "Gene Sets")
        self.gs_widget = GeneSetsSelection(gene_sets_box, self, 'stored_gene_sets_selection')
        self.gs_widget.hierarchy_tree_widget.itemClicked.connect(self.update_tree_view)

        self.commit_button = auto_commit(self.controlArea, self, "auto_commit", "&Commit", box=False)

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

    @staticmethod
    def set_items(gene_sets, sets_to_display, genes, ref, count_treshold, callback):
        model_items = []
        if not genes:
            return

        for gene_set in gene_sets:
            if gene_set.hierarchy not in sets_to_display:
                continue
            enrichemnt_result = gene_set.set_enrichment(ref, genes.intersection(ref))
            callback()

            if len(enrichemnt_result.query) >= count_treshold:
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

    def init_item_model(self):
        if self.data_model:
            self.data_model.clear()
            self.setup_filter_model()
        else:
            self.data_model = QStandardItemModel()

        self.data_model.setSortRole(Qt.UserRole)
        self.data_model.setHorizontalHeaderLabels(self.DATA_HEADER_LABELS)

    def sizeHint(self):
        return QSize(1280, 960)


if __name__ == "__main__":
    from AnyQt.QtWidgets import QApplication
    app = QApplication([])
    ow = OWGeneSets()
    ow.show()
    app.exec_()
