""" GeneSets """
import threading
import concurrent.futures

from functools import partial
from typing import Optional

from AnyQt.QtWidgets import (
    QTreeView, QHeaderView, QHBoxLayout,
)
from AnyQt.QtCore import (
    Qt, QSize, QThreadPool, Slot, QThread, QItemSelection, QItemSelectionRange, QItemSelectionModel
)
from AnyQt.QtGui import (
    QColor, QStandardItemModel, QStandardItem,
)

from Orange.widgets.gui import (
    vBox, lineEdit, ProgressBar, LinkRole, LinkStyledItemDelegate,
    auto_commit, widgetLabel, spin, comboBox, widgetBox
)
from Orange.data import Domain, Table, DiscreteVariable, StringVariable, filter as table_filter
from Orange.widgets.widget import OWWidget, Msg
from Orange.widgets.settings import Setting
from Orange.widgets.utils.signals import Output, Input
from Orange.widgets.utils.itemmodels import DomainModel
from Orange.widgets.utils.concurrent import ThreadExecutor, FutureWatcher, methodinvoke

from orangecontrib.bioinformatics.widgets.utils.data import (
    TAX_ID, GENE_AS_ATTRIBUTE_NAME, GENE_ID_COLUMN, GENE_ID_ATTRIBUTE,
    ERROR_ON_MISSING_ANNOTATION, ERROR_ON_MISSING_GENE_ID, ERROR_ON_MISSING_TAX_ID
)

from orangecontrib.bioinformatics.widgets.utils.gui import GeneSetsSelection
from orangecontrib.bioinformatics.widgets.utils.gui import NumericalColumnDelegate, FilterProxyModel
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
    name = "Gene Sets"
    description = ""
    icon = "icons/OWGeneSets.svg"
    priority = 9
    want_main_area = True

    COUNT, GENES, CATEGORY, TERM = range(4)
    DATA_HEADER_LABELS = ["Count", 'Genes In Set', 'Category', 'Term']

    organism = Setting(None, schema_only=True)
    stored_gene_sets_selection = Setting([], schema_only=True)
    selected_rows = Setting([], schema_only=True)
    custom_gene_set_indicator = Setting(None, schema_only=True)

    min_count = Setting(5)
    use_min_count = Setting(True)
    auto_commit = Setting(True)

    class Inputs:
        genes = Input("Data", Table)
        custom_sets = Input('Custom Gene Sets', Table)

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
        self.custom_gs_col_box = None
        self.gs_label_combobox = None
        self.custom_tax_id = None
        self.custom_use_attr_names = None
        self.custom_gene_id_attribute = None
        self.custom_gene_id_column = None
        self.num_of_custom_sets = None

        # Gene Sets widget
        self.gs_widget = None

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
                    # self.gs_label_combobox.setDisabled(True)
                    self.Error.organism_mismatch()
                    self.gs_widget.update_gs_hierarchy()
                    return

                if isinstance(self.custom_gene_set_indicator, DiscreteVariable):
                    labels = self.custom_gene_set_indicator.values
                    gene_sets_names = [labels[int(idx)] for idx
                                       in self.custom_data.get_column_view(self.custom_gene_set_indicator)[0]]
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
                # self.gs_label_combobox.setDisabled(False)
            else:
                self.gs_widget.update_gs_hierarchy()

        self.update_info_box()

    def update_tree_view(self):
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

    @Inputs.custom_sets
    def handle_custom_input(self, data):
        self.Error.clear()
        self.__reset_widget_state()
        self.custom_data = None
        self.custom_tax_id = None
        self.custom_use_attr_names = None
        self.custom_gene_id_attribute = None
        self.custom_gene_id_column = None
        self.feature_model.set_domain(None)

        if data:
            self.custom_data = data
            self.feature_model.set_domain(self.custom_data.domain)
            self.custom_tax_id = str(self.custom_data.attributes.get(TAX_ID, None))
            self.custom_use_attr_names = self.custom_data.attributes.get(GENE_AS_ATTRIBUTE_NAME, None)
            self.custom_gene_id_attribute = self.custom_data.attributes.get(GENE_ID_ATTRIBUTE, None)
            self.custom_gene_id_column = self.custom_data.attributes.get(GENE_ID_COLUMN, None)

            if self.gs_label_combobox is None:
                self.gs_label_combobox = comboBox(self.custom_gs_col_box, self, "custom_gene_set_indicator",
                                                  sendSelectedValue=True, model=self.feature_model,
                                                  callback=self.on_gene_set_indicator_changed)
            self.custom_gs_col_box.show()

            if self.custom_gene_set_indicator in self.feature_model:
                index = self.feature_model.indexOf(self.custom_gene_set_indicator)
                self.custom_gene_set_indicator = self.feature_model[index]
            else:
                self.custom_gene_set_indicator = self.feature_model[0]
        else:
            self.custom_gs_col_box.hide()

        self.gs_widget.clear_custom_sets()
        self.handle_custom_gene_sets(select_customs_flag=self.custom_gene_set_indicator is not None)
        self.invalidate()

    @Inputs.genes
    def handle_genes_input(self, data):
        self.Error.clear()
        self.__reset_widget_state()
        # clear output
        self.Outputs.matched_genes.send(None)
        # clear input values
        self.input_genes = []
        self.input_data = None
        self.tax_id = None
        self.use_attr_names = None
        self.gene_id_attribute = None
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

            self.gs_widget.load_gene_sets(self.tax_id)

            # if input data change, we need to refresh custom sets
            if self.custom_data:
                self.gs_widget.clear_custom_sets()
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

    def create_partial(self):
        return partial(self.set_items,
                       self.gs_widget.gs_object,
                       self.stored_gene_sets_selection,
                       set(self.input_genes),
                       self.callback)

    def callback(self):
        if self._task.cancelled:
            raise KeyboardInterrupt()
        if self.progress_bar:
            methodinvoke(self, "progress_advance")()

    def init_gene_sets(self):
        if self._task is not None:
            self.cancel()
        assert self._task is None

        self._task = Task()
        self.init_item_model()

        # save setting on selected hierarchies
        self.stored_gene_sets_selection = self.gs_widget.get_hierarchies(only_selected=True)

        f = self.create_partial()

        progress_iterations = sum([len(g_set) for hier, g_set
                                   in self.gs_widget.gs_object.map_hierarchy_to_sets().items()
                                   if hier in self.stored_gene_sets_selection])

        self.progress_bar = ProgressBar(self, iterations=progress_iterations)

        self._task.future = self._executor.submit(f)

        self._task.watcher = FutureWatcher(self._task.future)
        self._task.watcher.done.connect(self._init_gene_sets_finished)

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
            self.filter_data_view()
            self.set_selection()
            self.update_info_box()
        except Exception as ex:
            print(ex)

    def create_filters(self):
        search_term = self.search_pattern.lower().strip().split()

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

        return filters

    def filter_data_view(self):
        filter_proxy = self.filter_proxy_model  # type: FilterProxyModel
        model = filter_proxy.sourceModel()      # type: QStandardItemModel

        if isinstance(model, QStandardItemModel):

            # apply filtering rules
            filter_proxy.set_filters(self.create_filters())

            if model.rowCount() and not filter_proxy.rowCount():
                self.Warning.all_sets_filtered()
            else:
                self.Warning.clear()

    def set_selection(self):
        if len(self.selected_rows):
            view = self.data_view
            model = self.data_model

            row_model_indexes = [model.indexFromItem(model.item(i)) for i in self.selected_rows]
            proxy_rows = [self.filter_proxy_model.mapFromSource(i).row() for i in row_model_indexes]

            if model.rowCount() <= self.selected_rows[-1]:
                return

            header_count = view.header().count() - 1
            selection = QItemSelection()

            for row_index in proxy_rows:
                selection.append(
                    QItemSelectionRange(
                        self.filter_proxy_model.index(row_index, 0),
                        self.filter_proxy_model.index(row_index, header_count)
                    )
                )

            view.selectionModel().select(selection, QItemSelectionModel.ClearAndSelect)

    def commit(self):
        selection_model = self.data_view.selectionModel()

        if selection_model:
            selection = selection_model.selectedRows(self.COUNT)
            self.selected_rows = [self.filter_proxy_model.mapToSource(sel).row() for sel in selection]

            if selection and self.input_genes:
                genes = [model_index.data(Qt.UserRole) for model_index in selection]
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
                    # create filter from selected column for genes
                    only_known = table_filter.FilterStringList(self.gene_id_column, output_genes)
                    # apply filter to the data
                    data_table = table_filter.Values([only_known])(self.input_data)

                    self.Outputs.matched_genes.send(data_table)

    def assign_delegates(self):
        self.data_view.setItemDelegateForColumn(
            self.GENES, NumericalColumnDelegate(self)
        )

        self.data_view.setItemDelegateForColumn(
            self.COUNT, NumericalColumnDelegate(self)
        )

    def setup_filter_model(self):
        self.filter_proxy_model = FilterProxyModel()
        self.filter_proxy_model.setFilterKeyColumn(self.TERM)
        self.data_view.setModel(self.filter_proxy_model)

    def setup_filter_area(self):
        h_layout = QHBoxLayout()
        h_layout.setSpacing(100)
        h_widget = widgetBox(self.mainArea, orientation=h_layout)

        spin(h_widget, self, 'min_count', 0, 1000,
             label='Count',
             tooltip='Minimum genes count',
             checked='use_min_count',
             callback=self.filter_data_view,
             callbackOnReturn=True,
             checkCallback=self.filter_data_view)

        self.line_edit_filter = lineEdit(h_widget, self, 'search_pattern')
        self.line_edit_filter.setPlaceholderText('Filter gene sets ...')
        self.line_edit_filter.textChanged.connect(self.filter_data_view)

    def on_gene_set_indicator_changed(self):
        # self._handle_future_model()
        self.gs_widget.clear_custom_sets()
        self.handle_custom_gene_sets()
        self.invalidate()

    def setup_control_area(self):
        # Control area
        self.input_info = widgetLabel(
            widgetBox(self.controlArea, "Info", addSpace=True), 'No data on input.\n'
        )
        self.custom_gs_col_box = box = vBox(self.controlArea, 'Custom Gene Set Term Column')
        box.hide()

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
    def set_items(gene_sets, sets_to_display, genes, callback):
        model_items = []
        if not genes:
            return

        for gene_set in sorted(gene_sets):
            if gene_set.hierarchy not in sets_to_display:
                continue

            callback()

            matched_set = gene_set.genes & genes
            if len(matched_set) > 0:
                category_column = QStandardItem()
                term_column = QStandardItem()
                count_column = QStandardItem()
                genes_column = QStandardItem()

                category_column.setData(", ".join(gene_set.hierarchy), Qt.DisplayRole)
                term_column.setData(gene_set.name, Qt.DisplayRole)
                term_column.setData(gene_set.name, Qt.ToolTipRole)
                term_column.setData(gene_set.link, LinkRole)
                term_column.setForeground(QColor(Qt.blue))

                count_column.setData(matched_set, Qt.UserRole)
                count_column.setData(len(matched_set), Qt.DisplayRole)

                genes_column.setData(len(gene_set.genes), Qt.DisplayRole)
                genes_column.setData(set(gene_set.genes), Qt.UserRole)  # store genes to get then on output on selection

                model_items.append([count_column, genes_column, term_column, category_column])

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
