""" GeneSets """
import threading
import concurrent.futures


from functools import partial
from collections import defaultdict

from AnyQt.QtWidgets import (
    QTreeView, QTreeWidget, QTreeWidgetItem, QTreeWidgetItemIterator, QHeaderView
)
from AnyQt.QtCore import (
    Qt, QSize, QThreadPool, QSortFilterProxyModel, Slot, QThread
)
from AnyQt.QtGui import (
    QColor, QStandardItemModel, QStandardItem,
)

from Orange.widgets.gui import (
    vBox, lineEdit, ProgressBar, LinkRole, LinkStyledItemDelegate,
    auto_commit, widgetLabel, spin, comboBox, widgetBox
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

from orangecontrib.bioinformatics.widgets.utils.gui import NumericalColumnDelegate
from orangecontrib.bioinformatics.widgets.utils.concurrent import Worker
from orangecontrib.bioinformatics.widgets.utils.settings import OrganismContextHandler
from orangecontrib.bioinformatics.utils import serverfiles
from orangecontrib.bioinformatics import geneset


def hierarchy_tree(gene_sets):
    def tree():
        return defaultdict(tree)

    collection = tree()

    def collect(col, set_hierarchy):
        if set_hierarchy:
            collect(col[set_hierarchy[0]], set_hierarchy[1:])

    for hierarchy in gene_sets:
        collect(collection, hierarchy)

    return collection


def download_gene_sets(tax_id, gene_sets):

    # get only those sets that are not already downloaded
    for hierarchy in [hierarchy for hierarchy in gene_sets]:
        serverfiles.localpath_download(geneset.DOMAIN, geneset.filename(hierarchy, tax_id))
    return gene_sets


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
    settingsHandler = OrganismContextHandler()

    # settings
    auto_commit = Setting(True)
    stored_selections = ContextSetting([])
    organism = ContextSetting(None)
    matched_treshold = Setting(5)

    GENES, COUNT, CATEGORY, TERM,  = range(4)
    DATA_HEADER_LABELS = ['Genes In Set', 'Count', 'Category', 'Term']
    HIERARCHY_HEADER_LABELS = ['Category']

    class Inputs:
        genes = Input("Genes", Table)
        custom_sets = Input('Custom Gene Sets', Table)

    class Outputs:
        matched_genes = Output("Matched Genes", Table)

    class Information(OWWidget.Information):
        pass

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

        # gene sets object
        self.gene_sets_obj = geneset.GeneSets()

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
        self.gene_set_label = None
        self.gs_label_combobox = None
        self.custom_tax_id = None
        self.custom_use_attr_names = None
        self.custom_gene_id_attribute = None
        self.custom_gene_id_column = None

        # info box
        self.input_info = None
        self.num_of_sel_genes = 0

        # filter
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
        # reset hierarchy widget state
        self.hierarchy_widget.clear()
        # clear data view
        self.init_item_model()

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

    def __construct_custom_gene_sets(self):
        custom_set_hier = ('Custom sets',)

        # delete any custom sets if they exists
        self.gene_sets_obj.delete_sets_by_hierarchy(custom_set_hier)

        if self.custom_data and self.custom_gene_id_column:

            gene_sets_names, _ = self.custom_data.get_column_view(self.gene_set_label)
            gene_names, _ = self.custom_data.get_column_view(self.custom_gene_id_column)

            temp_dict = defaultdict(list)
            for set_name, gene_name in zip(gene_sets_names, gene_names):
                temp_dict[set_name].append(gene_name)

            g_sets = []
            for key, value in temp_dict.items():
                g_sets.append(geneset.GeneSet(gs_id=key,
                                              hierarchy=custom_set_hier,
                                              organism=self.custom_tax_id,
                                              name=key,
                                              genes=set(value)))

            self.gene_sets_obj.update(g_sets)

    def __update_hierarchy(self):
        self.set_hierarchy_model(self.hierarchy_widget, hierarchy_tree(self.gene_sets_obj.hierarchies()))
        self.set_selected_hierarchies()

    def update_tree_view(self):
        self.init_gene_sets()

    def invalidate(self):
        # clear
        self.__reset_widget_state()
        self.update_info_box()

        if self.input_data is not None:
            # setup
            self.__construct_custom_gene_sets()
            self.__get_input_genes()
            self.__update_hierarchy()
            self.update_tree_view()

    def __check_organism_mismatch(self):
        """ Check if organisms from different inputs match.

        :return: True if there is a mismatch
        """
        if self.tax_id is not None and self.custom_tax_id is not None:
            return self.tax_id != self.custom_tax_id
        return False

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
            self.custom_tax_id = str(self.custom_data.attributes.get(TAX_ID, None))
            self.custom_use_attr_names = self.custom_data.attributes.get(GENE_AS_ATTRIBUTE_NAME, None)
            self.custom_gene_id_attribute = self.custom_data.attributes.get(GENE_ID_ATTRIBUTE, None)
            self.custom_gene_id_column = self.custom_data.attributes.get(GENE_ID_COLUMN, None)

            if not (self.custom_use_attr_names is not None
                    and ((self.custom_gene_id_attribute is None) ^ (self.custom_gene_id_column is None))):

                if self.custom_tax_id is None:
                    self.Error.missing_annotation()
                    return

                self.Error.missing_gene_id()
                return

            elif self.custom_tax_id is None:
                self.Error.missing_tax_id()
                return

            if self.__check_organism_mismatch():
                self.Error.organism_mismatch()
                return

            self.gs_label_combobox.setDisabled(False)
            self.feature_model.set_domain(self.custom_data.domain)

            if self.feature_model:
                self.gene_set_label = self.feature_model[0]

        self.invalidate()

    @Inputs.genes
    def handle_genes_input(self, data):
        self.closeContext()
        self.Error.clear()
        self.__reset_widget_state()
        # clear output
        self.Outputs.matched_genes.send(None)
        # clear input genes
        self.input_genes = []
        self.gs_label_combobox.setDisabled(True)
        self.update_info_box()

        if data:
            self.input_data = data
            self.tax_id = str(self.input_data.attributes.get(TAX_ID, None))
            self.use_attr_names = self.input_data.attributes.get(GENE_AS_ATTRIBUTE_NAME, None)
            self.gene_id_attribute = self.input_data.attributes.get(GENE_ID_ATTRIBUTE, None)
            self.gene_id_column = self.input_data.attributes.get(GENE_ID_COLUMN, None)

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

            # if input data change, we need to set feature model again
            if self.custom_data:
                self.gs_label_combobox.setDisabled(False)
                self.feature_model.set_domain(self.custom_data.domain)

                if self.feature_model:
                    self.gene_set_label = self.feature_model[0]

            self.download_gene_sets()

    def update_info_box(self):
        info_string = ''
        if self.input_genes:
            info_string += '{} unique gene names on input.\n'.format(len(self.input_genes))
            info_string += '{} genes on output.\n'.format(self.num_of_sel_genes)
        else:
            info_string += 'No genes on input.\n'

        self.input_info.setText(info_string)

    def on_gene_sets_download(self, result):
        # make sure this happens in the main thread.
        # Qt insists that widgets be created within the GUI(main) thread.
        assert threading.current_thread() == threading.main_thread()
        self.setStatusMessage('')

        if result:
            for res in result:
                g_sets = geneset.load_gene_sets(res, self.tax_id)
                self.gene_sets_obj.update([g_set for g_set in g_sets])

        # add custom sets if there are any
        self.invalidate()
        self.update_info_box()

    def download_gene_sets(self):
        self.Error.clear()
        # reset hierarchy widget state
        self.hierarchy_widget.clear()
        # clear data view
        self.init_item_model()

        # get all gene sets for selected organism
        gene_sets = geneset.list_all(organism=self.tax_id)
        # status message
        self.setStatusMessage('downloading sets')

        worker = Worker(download_gene_sets, self.tax_id, gene_sets)
        worker.signals.result.connect(self.on_gene_sets_download)

        # move download process to worker thread
        self.threadpool.start(worker)

    def set_hierarchy_model(self, tree_widget, sets):

        def beautify_displayed_text(text):
            if '_' in text:
                return text.replace('_', ' ').title()
            else:
                return text

        # TODO: maybe optimize this code?
        for key, value in sets.items():
            item = QTreeWidgetItem(tree_widget, [beautify_displayed_text(key)])
            item.setFlags(item.flags() & (Qt.ItemIsUserCheckable | ~Qt.ItemIsSelectable | Qt.ItemIsEnabled))
            item.setExpanded(True)
            item.hierarchy = key

            if value:
                item.setFlags(item.flags() | Qt.ItemIsTristate)
                self.set_hierarchy_model(item, value)
            else:
                if item.parent():
                    item.hierarchy = (item.parent().hierarchy, key)

            if not item.childCount() and not item.parent():
                item.hierarchy = (key, )

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
            reference_genes = self.gene_sets_obj.genes()

        self.init_item_model()

        sets_to_display = self.get_hierarchies(only_selected=True)
        # save setting on selected hierarchies
        self.stored_selections = sets_to_display
        # save context
        self.closeContext()

        f = partial(self.set_items,
                    self.gene_sets_obj,
                    sets_to_display,
                    set(self.input_genes),
                    self.matched_treshold,
                    reference_genes,
                    callback=callback)

        progress_iterations = sum([len(g_set) for hier, g_set in self.gene_sets_obj.map_hierarchy_to_sets().items()
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
        except Exception as ex:
            print(ex)

    def set_selected_hierarchies(self):
        iterator = QTreeWidgetItemIterator(self.hierarchy_widget, QTreeWidgetItemIterator.All)

        while iterator.value():
            # note: if hierarchy value is not a tuple, then this is just top level qTreeWidgetItem that
            #       holds subcategories. We don't want to display all sets from category
            if type(iterator.value().hierarchy) is not str:
                if iterator.value().hierarchy in self.stored_selections:
                    iterator.value().setCheckState(0, Qt.Checked)
                else:
                    iterator.value().setCheckState(0, Qt.Unchecked)

            iterator += 1

        # if no items are checked, we check first one at random
        if len(self.get_hierarchies(only_selected=True)) == 0:
            iterator = QTreeWidgetItemIterator(self.hierarchy_widget, QTreeWidgetItemIterator.NotChecked)

            while iterator.value():
                if type(iterator.value().hierarchy) is not str:
                    iterator.value().setCheckState(0, Qt.Checked)
                    return

                iterator += 1

    def get_hierarchies(self, **kwargs):
        """ return selected hierarchy
        """
        only_selected = kwargs.get('only_selected', None)

        sets_to_display = list()

        if only_selected:
            iterator = QTreeWidgetItemIterator(self.hierarchy_widget, QTreeWidgetItemIterator.Checked)
        else:
            iterator = QTreeWidgetItemIterator(self.hierarchy_widget)

        while iterator.value():
            # note: if hierarchy value is not a tuple, then this is just top level qTreeWidgetItem that
            #       holds subcategories. We don't want to display all sets from category
            if type(iterator.value().hierarchy) is not str:

                if not only_selected:
                    sets_to_display.append(iterator.value().hierarchy)
                else:
                    if not iterator.value().isDisabled():
                        sets_to_display.append(iterator.value().hierarchy)

            iterator += 1

        return sets_to_display

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

    def setup_filter_model(self):
        self.filter_proxy_model = QSortFilterProxyModel()
        self.filter_proxy_model.setFilterKeyColumn(self.TERM)

    def setup_filter_area(self):
        line_edit_filter = lineEdit(self.mainArea, self, 'search_pattern', 'Filter gene sets:')
        line_edit_filter.setPlaceholderText('search pattern ...')
        line_edit_filter.textChanged.connect(self.filter_proxy_model.setFilterRegExp)

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

        hierarchy_box = widgetBox(self.controlArea, "Gene Set Databases")
        self.hierarchy_widget = QTreeWidget(self)
        self.hierarchy_widget.setEditTriggers(QTreeView.NoEditTriggers)
        self.hierarchy_widget.setHeaderLabels(self.HIERARCHY_HEADER_LABELS)
        self.hierarchy_widget.itemClicked.connect(self.update_tree_view)
        hierarchy_box.layout().addWidget(self.hierarchy_widget)

        self.commit_button = auto_commit(self.controlArea, self, "auto_commit", "&Commit", box=False)

    def setup_gui(self):
        # control area
        self.setup_control_area()

        # main area
        self.setup_filter_model()
        self.setup_filter_area()

        self.data_view = QTreeView()
        self.data_view.setModel(self.filter_proxy_model)
        self.data_view.setAlternatingRowColors(True)
        self.data_view.sortByColumn(self.COUNT, Qt.DescendingOrder)
        self.data_view.setSortingEnabled(True)
        self.data_view.setSelectionMode(QTreeView.ExtendedSelection)
        self.data_view.setEditTriggers(QTreeView.NoEditTriggers)
        self.data_view.viewport().setMouseTracking(False)
        self.data_view.setItemDelegateForColumn(self.TERM, LinkStyledItemDelegate(self.data_view))
        self.data_view.selectionModel().selectionChanged.connect(self.commit)

        self.mainArea.layout().addWidget(self.data_view)

        self.data_view.header().setSectionResizeMode(QHeaderView.ResizeToContents)
        self.assign_delegates()

    @staticmethod
    def set_items(gene_sets, sets_to_display, genes, matched_treshold, ref, callback):
        model_items = []
        if not genes:
            return

        for gene_set in gene_sets:
            if gene_set.hierarchy not in sets_to_display:
                continue

            category_column = QStandardItem()
            name_column = QStandardItem()
            matched_column = QStandardItem()
            genes_column = QStandardItem()

            category_column.setData(", ".join(gene_set.hierarchy), Qt.DisplayRole)
            name_column.setData(gene_set.name, Qt.DisplayRole)
            name_column.setData(gene_set.name, Qt.ToolTipRole)
            name_column.setData(gene_set.link, LinkRole)
            name_column.setForeground(QColor(Qt.blue))

            matched_set = gene_set.genes & genes
            matched_column.setData(matched_set, Qt.UserRole)
            matched_column.setData(len(matched_set), Qt.DisplayRole)

            genes_column.setData(len(gene_set.genes), Qt.DisplayRole)
            genes_column.setData(gene_set.genes, Qt.UserRole)  # store genes to get then on output on selection

            callback()

            if len(matched_set) >= matched_treshold:
                model_items.append([genes_column, matched_column, name_column, category_column])

        return model_items

    def init_item_model(self):
        if self.data_model:
            self.data_model.clear()
            self.filter_proxy_model.setSourceModel(None)
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
