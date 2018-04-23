""" GeneSets """
import threading
import numpy as np

from requests.exceptions import ConnectionError
from collections import defaultdict

from AnyQt.QtWidgets import (
    QTreeView, QTreeWidget, QTreeWidgetItem, QTreeWidgetItemIterator
)
from AnyQt.QtCore import (
    Qt, QSize, QThreadPool, QSortFilterProxyModel
)
from AnyQt.QtGui import (
    QColor, QStandardItemModel, QStandardItem
)

from Orange.widgets.gui import (
    vBox, comboBox, lineEdit, ProgressBar, widgetBox, LinkRole, LinkStyledItemDelegate,
    auto_commit, widgetLabel, checkBox, attributeItem
)

from Orange.widgets.widget import OWWidget, Msg
from Orange.widgets.utils import itemmodels
from Orange.widgets.settings import Setting, ContextSetting
from Orange.widgets.utils.signals import Output, Input

from Orange.data import DiscreteVariable, StringVariable, Domain, Table

from orangecontrib.bioinformatics.widgets.utils.data import TAX_ID, GENE_AS_ATTRIBUTE_NAME
from orangecontrib.bioinformatics.widgets.utils.concurrent import Worker
from orangecontrib.bioinformatics.widgets.utils.settings import OrganismContextHandler
from orangecontrib.bioinformatics.utils import serverfiles
from orangecontrib.bioinformatics.ncbi import gene, taxonomy
from orangecontrib.bioinformatics import geneset


CATEGORY, GENES, MATCHED, TERM = range(4)
DATA_HEADER_LABELS = ["Category", "Genes", "Matched", "Term"]
HIERARCHY_HEADER_LABELS = ["Category"]


def hierarchy_tree(tax_id, gene_sets):
    def tree():
        return defaultdict(tree)

    collection = tree()

    def collect(col, set_hierarchy):
        if set_hierarchy:
            collect(col[set_hierarchy[0]], set_hierarchy[1:])

    for hierarchy, t_id in gene_sets:
        collect(collection[t_id], hierarchy)

    return tax_id, collection[tax_id]


def get_collections(gene_sets, mapped_genes, partial_result, progress_callback):

    for gene_set in gene_sets:
        category_column = QStandardItem()
        name_column = QStandardItem()
        matched_column = QStandardItem()
        genes_column = QStandardItem()

        category_column.setData(", ".join(gene_set.hierarchy), Qt.DisplayRole)
        name_column.setData(gene_set.name, Qt.DisplayRole)
        name_column.setData(gene_set.name, Qt.ToolTipRole)
        name_column.setData(gene_set.link, LinkRole)
        name_column.setForeground(QColor(Qt.blue))

        if mapped_genes:
            matched_set = gene_set.genes & {ncbi_id for input_name, ncbi_id in mapped_genes.items()}
            matched_column.setData(matched_set, Qt.UserRole)
            matched_column.setData(len(matched_set), Qt.DisplayRole)

        genes_column.setData(len(gene_set.genes), Qt.DisplayRole)
        genes_column.setData(gene_set.genes, Qt.UserRole)  # store genes to get then on output on selection

        progress_callback.emit()
        partial_result.emit([category_column, genes_column, matched_column, name_column])


def download_gene_sets(gene_sets, progress_callback):

    # get only those sets that are not already downloaded
    for hierarchy, tax_id in [(hierarchy, tax_id) for hierarchy, tax_id in gene_sets]:

        serverfiles.localpath_download(geneset.DOMAIN, geneset.filename(hierarchy, tax_id),
                                       callback=progress_callback.emit)

    return tax_id, gene_sets


class OWGeneSets(OWWidget):
    name = "Gene Sets"
    description = ""
    icon = "icons/OWGeneSets.svg"
    priority = 9
    want_main_area = True

    # settings
    selected_organism = Setting(0)
    auto_commit = Setting(True)
    auto_apply = Setting(True)
    gene_col_index = Setting(0)
    use_attr_names = Setting(False)

    settingsHandler = OrganismContextHandler()
    stored_selections = ContextSetting([])
    organism = ContextSetting(None)

    class Inputs:
        genes = Input("Genes", Table)

    class Outputs:
        matched_genes = Output("Matched Genes", Table)

    class Information(OWWidget.Information):
        pass

    class Error(OWWidget.Error):
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
        self.tax_id = None
        self.input_genes = None
        self.mapped_genes = None
        self.organisms = list()
        self.input_info = None

        self.column_candidates = []

        # filter
        self.lineEdit_filter = None
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

        # input options
        self.gene_columns = None
        self.gene_column_combobox = None
        self.attr_names_checkbox = None

        # threads
        self.threadpool = QThreadPool(self)
        self.workers = None

        # gui
        self.setup_gui()
        self._get_available_organisms()

        # self.handle_input(self.input_genes)
        # self.on_organism_change()

    def _progress_advance(self):
        # GUI should be updated in main thread. That's why we are calling advance method here
        if self.progress_bar:
            self.progress_bar.advance()

    def _get_available_organisms(self):
        available_organism = sorted([(tax_id, taxonomy.name(tax_id)) for tax_id in taxonomy.common_taxids()],
                                    key=lambda x: x[1])

        self.organisms = [tax_id[0] for tax_id in available_organism]

        self.organism_select_combobox.addItems([tax_id[1] for tax_id in available_organism])

    def _gene_names_from_table(self):
        """ Extract and return gene names from `Orange.data.Table`.
        """
        self.input_genes = []
        if self.input_data:
            if self.use_attr_names:
                self.input_genes = [str(attr.name).strip() for attr in self.input_data.domain.attributes]
            elif self.gene_columns:
                column = self.gene_columns[self.gene_col_index]
                self.input_genes = [str(e[column]) for e in self.input_data if not np.isnan(e[column])]

    def _update_gene_matcher(self):
        self._gene_names_from_table()
        if not self.gene_matcher:
            self.gene_matcher = gene.GeneMatcher(self.get_selected_organism())

        self.gene_matcher.genes = self.input_genes
        self.gene_matcher.organism = self.get_selected_organism()

    def get_selected_organism(self):
        return self.organisms[self.selected_organism]

    def on_input_option_change(self):
        self._update_gene_matcher()
        self.match_genes()

    @Inputs.genes
    def handle_input(self, data):
        self.closeContext()

        if data:
            self.input_data = data

            self.gene_column_combobox.clear()
            self.column_candidates = [attr for attr in data.domain.variables + data.domain.metas
                                      if isinstance(attr, (StringVariable, DiscreteVariable))]

            for var in self.column_candidates:
                self.gene_column_combobox.addItem(*attributeItem(var))

            self.tax_id = str(self.input_data.attributes.get(TAX_ID, ''))
            self.use_attr_names = self.input_data.attributes.get(GENE_AS_ATTRIBUTE_NAME, self.use_attr_names)

            self.gene_col_index = min(self.gene_col_index, len(self.column_candidates) - 1)

            if self.tax_id in self.organisms:
                self.selected_organism = self.organisms.index(self.tax_id)

            self.openContext(self.tax_id)

        self.on_input_option_change()

    def update_info_box(self):
        info_string = ''
        if self.input_genes:
            info_string += '{} unique gene names on input.\n'.format(len(self.input_genes))
            mapped = self.gene_matcher.get_known_genes()
            if mapped:
                ratio = (len(mapped) / len(self.input_genes)) * 100
                info_string += '{} ({:.2f}%) gene names matched.\n'.format(len(mapped), ratio)
        else:
            info_string += 'No genes on input.\n'

        self.input_info.setText(info_string)

    def match_genes(self):
        if self.gene_matcher:
            # init progress bar
            self.progress_bar = ProgressBar(self, iterations=len(self.gene_matcher.genes))
            # status message
            self.setStatusMessage('gene matcher running')

            worker = Worker(self.gene_matcher.run_matcher, progress_callback=True)
            worker.signals.progress.connect(self._progress_advance)
            worker.signals.finished.connect(self.handle_matcher_results)

            # move download process to worker thread
            self.threadpool.start(worker)

    def handle_matcher_results(self):
        assert threading.current_thread() == threading.main_thread()
        self.openContext(self.get_selected_organism())

        if self.progress_bar:
            self.progress_bar.finish()
            self.setStatusMessage('')

        if self.gene_matcher.map_input_to_ncbi():
            self.download_gene_sets()
            self.update_info_box()
        else:
            # reset gene sets
            self.init_item_model()
            self.hierarchy_widget.clear()
            self.update_info_box()

    def on_gene_sets_download(self, result):
        # make sure this happens in the main thread.
        # Qt insists that widgets be created within the GUI(main) thread.
        assert threading.current_thread() == threading.main_thread()
        self.progress_bar.finish()
        self.setStatusMessage('')

        tax_id, sets = result
        self.set_hierarchy_model(self.hierarchy_widget, *hierarchy_tree(tax_id, sets))
        self.set_selected_hierarchies()

        self.organism_select_combobox.setEnabled(True)  # re-enable combobox

        self.update_info_box()
        self.mapped_genes = self.gene_matcher.map_input_to_ncbi()
        self.workers = defaultdict(list)
        self.progress_bar_iterations = dict()

        for selected_hierarchy in self.get_hierarchies():
            gene_sets = geneset.load_gene_sets(selected_hierarchy)
            worker = Worker(get_collections, gene_sets, self.mapped_genes, progress_callback=True, partial_result=True)
            worker.signals.error.connect(self.handle_error)
            worker.signals.finished.connect(self.handle_worker_finished)
            worker.signals.progress.connect(self._progress_advance)
            worker.signals.partial_result.connect(self.populate_data_model)
            worker.setAutoDelete(False)

            self.workers[selected_hierarchy] = worker
            self.progress_bar_iterations[selected_hierarchy] = len(gene_sets)

        self.display_gene_sets()

    def handle_worker_finished(self):
        # We check if all workers have completed. If not, continue
        # dirty hax, is this ok?
        if self.progress_bar and self.progress_bar.widget.progressBarValue == 100:
            self.progress_bar.finish()
            self.setStatusMessage('')
            self.hierarchy_widget.setDisabled(False)

            # adjust column width
            for i in range(len(DATA_HEADER_LABELS) - 1):

                self.data_view.resizeColumnToContents(i)

            self.filter_proxy_model.setSourceModel(self.data_model)

    def populate_data_model(self, partial_result):
        assert threading.current_thread() == threading.main_thread()

        if partial_result:
            self.data_model.appendRow(partial_result)

    def set_hierarchy_model(self, model, tax_id, sets):

        def beautify_displayed_text(text):
            if '_' in text:
                return text.replace('_', ' ').title()
            else:
                return text

        # TODO: maybe optimize this code?
        for key, value in sets.items():
            item = QTreeWidgetItem(model, [beautify_displayed_text(key)])
            item.setFlags(item.flags() & (Qt.ItemIsUserCheckable | ~Qt.ItemIsSelectable | Qt.ItemIsEnabled))
            item.setExpanded(True)
            item.tax_id = tax_id
            item.hierarchy = key

            if value:
                item.setFlags(item.flags() | Qt.ItemIsTristate)
                self.set_hierarchy_model(item, tax_id, value)
            else:
                if item.parent():
                    item.hierarchy = ((item.parent().hierarchy, key), tax_id)

            if not item.childCount() and not item.parent():
                item.hierarchy = ((key,), tax_id)

    def download_gene_sets(self):
        tax_id = self.get_selected_organism()

        self.Error.clear()
        # do not allow user to change organism when download task is running
        self.organism_select_combobox.setEnabled(False)
        # reset hierarchy widget state
        self.hierarchy_widget.clear()
        # clear data view
        self.init_item_model()

        # get all gene sets for selected organism
        gene_sets = geneset.list_all(organism=tax_id)
        # init progress bar
        self.progress_bar = ProgressBar(self, iterations=len(gene_sets) * 100)
        # status message
        self.setStatusMessage('downloading sets')

        worker = Worker(download_gene_sets, gene_sets, progress_callback=True)
        worker.signals.progress.connect(self._progress_advance)
        worker.signals.result.connect(self.on_gene_sets_download)
        worker.signals.error.connect(self.handle_error)

        # move download process to worker thread
        self.threadpool.start(worker)

    def display_gene_sets(self):
        self.init_item_model()
        self.hierarchy_widget.setDisabled(True)

        only_selected_hier = self.get_hierarchies(only_selected=True)

        # init progress bar
        iterations = sum([self.progress_bar_iterations[hier] for hier in only_selected_hier])
        self.progress_bar = ProgressBar(self, iterations=iterations)
        self.setStatusMessage('displaying gene sets')

        if not only_selected_hier:
            self.progress_bar.finish()
            self.setStatusMessage('')
            self.hierarchy_widget.setDisabled(False)
            return

        # save setting on selected hierarchies
        self.stored_selections = only_selected_hier
        # save context
        self.closeContext()

        for selected_hierarchy in only_selected_hier:
            self.threadpool.start(self.workers[selected_hierarchy])

        self.openContext(self.get_selected_organism())

    def handle_error(self, ex):
        self.progress_bar.finish()
        self.setStatusMessage('')

        if isinstance(ex, ConnectionError):
            self.organism_select_combobox.setEnabled(True)  # re-enable combobox
            self.Error.cant_reach_host()

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
            matched_genes = selection_model.selectedRows(MATCHED)

            if matched_genes and self.input_genes:
                genes = [model_index.data(Qt.UserRole) for model_index in matched_genes]
                output_genes = [gene_name for gene_name in list(set.union(*genes))]
                input_to_ncbi = self.gene_matcher.map_input_to_ncbi()
                ncbi_to_input = {ncbi_id: input_name for input_name, ncbi_id
                                 in self.gene_matcher.map_input_to_ncbi().items()}

                if self.use_attr_names:
                    selected = [self.input_data.domain[ncbi_to_input[output_gene]] for output_gene in output_genes]
                    domain = Domain(selected, self.input_data.domain.class_vars, self.input_data.domain.metas)
                    new_data = self.input_data.from_table(domain, self.input_data)
                    self.Outputs.matched_genes.send(new_data)

                elif self.column_candidates:
                    column = self.column_candidates[self.gene_col_index]
                    selected_rows = []

                    for row_index, row in enumerate(self.input_data):
                        if str(row[column]) in input_to_ncbi.keys() and input_to_ncbi[str(row[column])] in output_genes:
                            selected_rows.append(row_index)
                    if selected_rows:
                        selected = self.input_data[selected_rows]
                    else:
                        selected = None
                    self.Outputs.matched_genes.send(selected)

    def setup_gui(self):
        # control area
        info_box = vBox(self.controlArea, 'Input info')
        self.input_info = widgetLabel(info_box)

        organism_box = vBox(self.controlArea, 'Organisms')
        self.organism_select_combobox = comboBox(organism_box, self,
                                                 'selected_organism',
                                                 callback=self.on_input_option_change)

        # Selection of genes attribute
        box = widgetBox(self.controlArea, 'Gene attribute')
        self.gene_columns = itemmodels.VariableListModel(parent=self)
        self.gene_column_combobox = comboBox(box, self, 'gene_col_index', callback=self.on_input_option_change)
        self.gene_column_combobox.setModel(self.gene_columns)

        self.attr_names_checkbox = checkBox(box, self, 'use_attr_names', 'Use attribute names',
                                            disables=[(-1, self.gene_column_combobox)],
                                            callback=self.on_input_option_change)

        self.gene_column_combobox.setDisabled(bool(self.use_attr_names))

        hierarchy_box = widgetBox(self.controlArea, "Entity Sets")
        self.hierarchy_widget = QTreeWidget(self)
        self.hierarchy_widget.setEditTriggers(QTreeView.NoEditTriggers)
        self.hierarchy_widget.setHeaderLabels(HIERARCHY_HEADER_LABELS)
        self.hierarchy_widget.itemClicked.connect(self.display_gene_sets)
        hierarchy_box.layout().addWidget(self.hierarchy_widget)

        self.commit_button = auto_commit(self.controlArea, self, "auto_commit", "&Commit", box=False)

        # rubber(self.controlArea)

        # main area
        self.filter_proxy_model = QSortFilterProxyModel(self.data_view)
        self.filter_proxy_model.setFilterKeyColumn(3)

        self.data_view = QTreeView()
        self.data_view.setModel(self.filter_proxy_model)
        self.data_view.setAlternatingRowColors(True)
        self.data_view.sortByColumn(2, Qt.DescendingOrder)
        self.data_view.setSortingEnabled(True)
        self.data_view.setSelectionMode(QTreeView.ExtendedSelection)
        self.data_view.setEditTriggers(QTreeView.NoEditTriggers)
        self.data_view.viewport().setMouseTracking(True)
        self.data_view.setItemDelegateForColumn(TERM, LinkStyledItemDelegate(self.data_view))

        self.data_view.selectionModel().selectionChanged.connect(self.commit)

        self.lineEdit_filter = lineEdit(self.mainArea, self, 'search_pattern', 'Filter gene sets:')
        self.lineEdit_filter.setPlaceholderText('search pattern ...')
        self.lineEdit_filter.textChanged.connect(self.filter_proxy_model.setFilterRegExp)

        self.mainArea.layout().addWidget(self.data_view)

    def init_item_model(self):
        if self.data_model:
            self.data_model.clear()
            self.filter_proxy_model.setSourceModel(None)
        else:
            self.data_model = QStandardItemModel()

        self.data_model.setSortRole(Qt.UserRole)
        self.data_model.setHorizontalHeaderLabels(DATA_HEADER_LABELS)

    def sizeHint(self):
        return QSize(1280, 960)


if __name__ == "__main__":
    from AnyQt.QtWidgets import QApplication
    app = QApplication([])
    ow = OWGeneSets()
    ow.show()
    app.exec_()
