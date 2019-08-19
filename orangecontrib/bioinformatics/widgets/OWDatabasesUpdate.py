""" Databases update widget """
import os
import sys
import json
import threading
from shutil import copyfile
from datetime import datetime as d_time
from functools import partial
from collections import OrderedDict, namedtuple

from requests.exceptions import Timeout, ConnectionError

from AnyQt.QtCore import Qt, QSize, Signal, QThreadPool
from AnyQt.QtWidgets import (
    QLabel,
    QDialog,
    QWidget,
    QCheckBox,
    QComboBox,
    QLineEdit,
    QFileDialog,
    QHBoxLayout,
    QPushButton,
    QToolButton,
    QTreeWidget,
    QVBoxLayout,
    QApplication,
    QTreeWidgetItem,
    QDialogButtonBox,
    QAbstractItemView,
    QStyledItemDelegate,
)

from Orange.widgets import gui
from Orange.widgets.widget import OWWidget

from orangecontrib.bioinformatics.utils import serverfiles
from orangecontrib.bioinformatics.geneset import DOMAIN as gene_sets_domain
from orangecontrib.bioinformatics.geneset import filename
from orangecontrib.bioinformatics.go.config import DOMAIN as gene_ontology_domain
from orangecontrib.bioinformatics.go.config import FILENAME_ANNOTATION
from orangecontrib.bioinformatics.ncbi.taxonomy import common_taxids, common_taxid_to_name, species_name_to_taxid
from orangecontrib.bioinformatics.widgets.utils.gui import TokenListCompleter
from orangecontrib.bioinformatics.widgets.utils.concurrent import Worker

# File states
AVAILABLE, CURRENT, OUTDATED, DEPRECATED, USER_FILE = range(5)
# File sources
SOURCE_SERVER = 'server_file'  # files on the serverfiles-bio repository
SOURCE_USER = 'user_file'  # user defined files
INFO_FILE_SCHEMA = {
    'domain': None,
    'filename': None,
    'source': None,
    'title': None,
    'tags': [],
    'size': None,
    'datetime': None
    # used only if files are compressed
    # 'uncompressed': None,
    # 'compression': None,
}


def file_size_bytes(file_path):
    """ returns file size in bytes """
    return os.stat(file_path).st_size


def create_info_file(file_path, **kwargs):
    info_dict = OrderedDict(INFO_FILE_SCHEMA)

    info_dict.update(**kwargs)
    info_dict['datetime'] = '{0:%Y-%m-%d %H:%M:%S.%f}'.format(d_time.today())
    info_dict['size'] = file_size_bytes(file_path)

    with open(file_path + '.info', 'wt') as f:
        json.dump(info_dict, f)


def create_folder(path):
    try:
        os.makedirs(path)
    except OSError:
        if os.path.exists(path):
            pass
        else:
            # There was an error on creation, so make sure we know about it
            raise


class UpdateOptionsItemDelegate(QStyledItemDelegate):
    """ An item delegate for the updates tree widget.

    note:
        Must be a child of a QTreeWidget.

    """

    def sizeHint(self, option, index):
        size = QStyledItemDelegate.sizeHint(self, option, index)
        parent = self.parent()
        item = parent.itemFromIndex(index)
        widget = parent.itemWidget(item, 0)
        if widget:
            size = QSize(size.width(), widget.sizeHint().height() / 2)
        return size


file_state = namedtuple('file_state', ['info_local', 'info_server', 'state'])

header_labels = ['', 'Title', 'Update', 'Updated', 'Size', 'Source']
header_index = namedtuple('header_index', ['Download'] + header_labels[1:])
header = header_index(*[index for index, _ in enumerate(header_labels)])


def UpdateItem_match(item, string):
    """
    Return `True` if the `UpdateItem` item contains a string in tags
    or in the title.

    """
    string = string.lower()
    return any(string.lower() in tag.lower() for tag in item.tags + [item.title])


def evaluate_all_info(local, server):
    """ Return FileState

    Args:
        local: info files from LocalFiles
        server: info files from ServerFiles

    """
    files = set(local.keys()).union(server.keys())

    for domain, file_name in sorted(files):
        yield FileState(domain, file_name, server.get((domain, file_name), None), local.get((domain, file_name), None))


def evaluate_files_state(progress_callback):
    progress_callback.emit()
    files = []

    # fetch remote info
    try:
        server_info = serverfiles.ServerFiles().allinfo()
    except (Timeout, ConnectionError) as e:
        raise e
    progress_callback.emit()

    # fetch local info
    local_info = serverfiles.allinfo()

    all_info = set(local_info.keys()).union(server_info.keys())

    for domain, file_name in sorted(all_info):
        files.append(
            FileState(
                domain,
                file_name,
                server_info.get((domain, file_name), None),
                local_info.get((domain, file_name), None),
            )
        )
    progress_callback.emit()
    return files


def download_server_file(fs, index, progress_callback):
    try:
        serverfiles.download(fs.domain, fs.filename, callback=progress_callback.emit)
    except Exception:
        # send FileState and index with Exception
        raise ValueError(fs, index)

    return fs, index


class OWDatabasesUpdate(OWWidget):

    name = "Databases Update"
    description = "Update local systems biology databases."
    icon = "../widgets/icons/OWDatabasesUpdate.svg"
    priority = 1

    inputs = []
    outputs = []

    want_main_area = False

    def __init__(self, parent=None, signalManager=None, name="Databases update"):
        OWWidget.__init__(self, parent, signalManager, name, wantMainArea=False)

        self.searchString = ""

        fbox = gui.widgetBox(self.controlArea, "Filter")
        self.completer = TokenListCompleter(self, caseSensitivity=Qt.CaseInsensitive)
        self.lineEditFilter = QLineEdit(textChanged=self.search_update)
        self.lineEditFilter.setCompleter(self.completer)

        fbox.layout().addWidget(self.lineEditFilter)

        box = gui.widgetBox(self.controlArea, "Files")
        self.filesView = QTreeWidget(self)
        self.filesView.setHeaderLabels(header_labels)
        self.filesView.setRootIsDecorated(False)
        self.filesView.setUniformRowHeights(True)
        self.filesView.setSelectionMode(QAbstractItemView.NoSelection)
        self.filesView.setSortingEnabled(True)
        self.filesView.sortItems(header.Title, Qt.AscendingOrder)
        self.filesView.setItemDelegateForColumn(0, UpdateOptionsItemDelegate(self.filesView))

        self.filesView.model().layoutChanged.connect(self.search_update)

        box.layout().addWidget(self.filesView)

        layout = QHBoxLayout()
        gui.widgetBox(self.controlArea, margin=0, orientation=layout)

        self.updateButton = gui.button(
            box, self, "Update all", callback=self.update_all, tooltip="Update all updatable files"
        )

        self.downloadButton = gui.button(
            box, self, "Download all", callback=self.download_filtered, tooltip="Download all filtered files shown"
        )

        self.cancelButton = gui.button(
            box, self, "Cancel", callback=self.cancel_active_threads, tooltip="Cancel scheduled downloads/updates."
        )

        self.addButton = gui.button(
            box, self, "Add ...", callback=self.__handle_dialog, tooltip="Add files for personal use."
        )

        layout.addWidget(self.updateButton)
        layout.addWidget(self.downloadButton)
        layout.addWidget(self.cancelButton)
        layout.addStretch()
        layout.addWidget(self.addButton)

        # Enable retryButton once connection is established
        # self.retryButton = gui.button(
        #     box, self, "Reconnect", callback=self.initialize_files_view
        # )
        # self.retryButton.hide()

        self.resize(800, 600)

        self.update_items = []
        self._dialog = None
        self.progress_bar = None

        # threads
        self.threadpool = QThreadPool(self)
        # self.threadpool.setMaxThreadCount(1)
        self.workers = []

        self.initialize_files_view()

    def __handle_dialog(self):
        if not self._dialog:
            self._dialog = FileUploadHelper(self)
        self._dialog.show()

    def __progress_advance(self):
        # GUI should be updated in main thread. That's why we are calling advance method here
        if self.progress_bar:
            self.progress_bar.advance()

    def handle_worker_exception(self, ex):
        self.progress_bar.finish()
        self.setStatusMessage('')

        if isinstance(ex, ConnectionError):
            # TODO: set warning messages
            pass

        print(ex)

    def initialize_files_view(self):
        # self.retryButton.hide()

        # clear view
        self.filesView.clear()
        # init progress bar
        self.progress_bar = gui.ProgressBar(self, iterations=3)
        # status message
        self.setStatusMessage('initializing')

        worker = Worker(evaluate_files_state, progress_callback=True)
        worker.signals.progress.connect(self.__progress_advance)
        worker.signals.result.connect(self.set_files_list)
        worker.signals.error.connect(self.handle_worker_exception)

        # move download process to worker thread
        self.threadpool.start(worker)
        self.setEnabled(False)

    def __create_action_button(self, fs, retry=None):
        if not fs.state not in [OUTDATED, USER_FILE] or not retry:
            self.filesView.setItemWidget(fs.tree_item, header.Update, None)

        button = QToolButton(None)
        if not retry:
            if fs.state == OUTDATED:
                button.setText('Update')
                button.clicked.connect(partial(self.submit_download_task, fs.domain, fs.filename, True))
            elif fs.state == USER_FILE:
                if not fs.info_server:
                    button.setText('Remove')
                    button.clicked.connect(partial(self.submit_remove_task, fs.domain, fs.filename))
                else:
                    button.setText('Use server version')
                    button.clicked.connect(partial(self.submit_download_task, fs.domain, fs.filename, True))
        else:
            button.setText('Retry')
            button.clicked.connect(partial(self.submit_download_task, fs.domain, fs.filename, True))

        button.setMaximumWidth(120)
        button.setMaximumHeight(20)
        button.setMinimumHeight(20)

        if sys.platform == "darwin":
            button.setAttribute(Qt.WA_MacSmallSize)

        self.filesView.setItemWidget(fs.tree_item, header.Update, button)

    def set_files_list(self, result):
        """ Set the files to show.
        """
        assert threading.current_thread() == threading.main_thread()
        self.progress_bar.finish()
        self.setStatusMessage('')
        self.setEnabled(True)

        self.update_items = result
        all_tags = set()

        for fs in self.update_items:
            fs.tree_item = FileStateItem(fs)
            fs.download_option = DownloadOption(state=fs.state)

            fs.download_option.download_clicked.connect(partial(self.submit_download_task, fs.domain, fs.filename))
            fs.download_option.remove_clicked.connect(partial(self.submit_remove_task, fs.domain, fs.filename))

        # add widget items to the QTreeWidget
        self.filesView.addTopLevelItems([fs.tree_item for fs in self.update_items])

        # add action widgets to tree items
        for fs in self.update_items:
            self.filesView.setItemWidget(fs.tree_item, header.Download, fs.download_option)
            if fs.state in [USER_FILE, OUTDATED]:
                self.__create_action_button(fs)

            all_tags.update(fs.tags)

        self.filesView.setColumnWidth(header.Download, self.filesView.sizeHintForColumn(header.Download))

        for column in range(1, len(header_labels)):
            self.filesView.resizeColumnToContents(column)

        hints = [hint for hint in sorted(all_tags) if not hint.startswith("#")]
        self.completer.setTokenList(hints)
        self.search_update()
        self.toggle_action_buttons()
        self.cancelButton.setEnabled(False)

    def toggle_action_buttons(self):
        selected_items = [fs for fs in self.update_items if not fs.tree_item.isHidden()]

        def button_check(sel_items, state, button):
            for item in sel_items:
                if item.state != state:
                    button.setEnabled(False)
                else:
                    button.setEnabled(True)
                    break

        button_check(selected_items, OUTDATED, self.updateButton)
        button_check(selected_items, AVAILABLE, self.downloadButton)

    def search_update(self, searchString=None):
        strings = str(self.lineEditFilter.text()).split()
        for fs in self.update_items:
            hide = not all(UpdateItem_match(fs, string) for string in strings)
            fs.tree_item.setHidden(hide)
        self.toggle_action_buttons()

    def update_all(self):
        for fs in self.update_items:
            if fs.state == OUTDATED and not fs.tree_item.isHidden():
                self.submit_download_task(fs.domain, fs.filename)

    def download_filtered(self):
        for fs in self.update_items:
            if not fs.tree_item.isHidden() and fs.state in [AVAILABLE, OUTDATED]:
                self.submit_download_task(fs.domain, fs.filename, start=False)

        self.run_download_tasks()

    def submit_download_task(self, domain, filename, start=True):
        """ Submit the (domain, filename) to be downloaded/updated.
        """
        # get selected tree item
        index = self.tree_item_index(domain, filename)
        fs = self.update_items[index]

        worker = Worker(download_server_file, fs, index, progress_callback=True)
        worker.signals.progress.connect(self.__progress_advance)
        worker.signals.result.connect(self.on_download_finished)
        worker.signals.error.connect(self.on_download_exception)

        self.workers.append(worker)

        if start:
            self.run_download_tasks()

    def run_download_tasks(self):
        self.cancelButton.setEnabled(True)
        # init progress bar

        self.progress_bar = gui.ProgressBar(self, iterations=len(self.workers) * 100)

        # status message
        self.setStatusMessage('downloading')

        # move workers to threadpool
        [self.threadpool.start(worker) for worker in self.workers]
        self.filesView.setDisabled(True)
        # reset list of workers
        self.workers = []

    def on_download_exception(self, ex):
        assert threading.current_thread() == threading.main_thread()
        self.progress_bar.finish()
        self.setStatusMessage('')
        print(ex)
        if isinstance(ex, ValueError):
            fs, index = ex.args

            # restore state and retry
            fs.refresh_state()
            fs.tree_item.update_data(fs)
            fs.download_option.state = fs.state
            self.__create_action_button(fs, retry=True)

    def on_download_finished(self, result):
        assert threading.current_thread() == threading.main_thread()

        # We check if all workers have completed. If not, continue
        if self.progress_bar.count == 100 or self.threadpool.activeThreadCount() == 0:
            self.filesView.setDisabled(False)
            self.progress_bar.finish()
            self.setStatusMessage('')

        fs, index = result
        # re-evaluate File State
        info = serverfiles.info(fs.domain, fs.filename)
        fs.refresh_state(info_local=info, info_server=info)
        # reinitialize treeWidgetItem
        fs.tree_item.update_data(fs)
        # reinitialize OptionWidget
        fs.download_option.state = fs.state
        self.filesView.setItemWidget(fs.tree_item, header.Update, None)

        self.toggle_action_buttons()
        for column in range(1, len(header_labels)):
            self.filesView.resizeColumnToContents(column)

    def submit_remove_task(self, domain, filename):
        serverfiles.LOCALFILES.remove(domain, filename)

        index = self.tree_item_index(domain, filename)
        fs = self.update_items[index]

        if fs.state == USER_FILE:
            self.filesView.takeTopLevelItem(self.filesView.indexOfTopLevelItem(fs.tree_item))
            self.update_items.remove(fs)
            # self.filesView.removeItemWidget(index)
        else:
            # refresh item state
            fs.info_local = None
            fs.refresh_state()
            # reinitialize treeWidgetItem
            fs.tree_item.update_data(fs)
            # reinitialize OptionWidget
            fs.download_option.state = fs.state

        self.toggle_action_buttons()

    def cancel_active_threads(self):
        """ Cancel all pending update/download tasks (that have not yet started).
        """
        if self.threadpool:
            self.threadpool.clear()

    def tree_item_index(self, domain, filename):
        for i, fs in enumerate(self.update_items):
            if fs.domain == domain and fs.filename == filename:
                return i
        raise ValueError("%r, %r not in update list" % (domain, filename))

    def onDeleteWidget(self):
        self.cancel_active_threads()
        OWWidget.onDeleteWidget(self)


class DownloadOption(QWidget):
    """ A Widget with download/update/remove options.
    """

    download_clicked = Signal()
    remove_clicked = Signal()

    def __init__(self, state=AVAILABLE, parent=None):
        QWidget.__init__(self, parent)
        layout = QHBoxLayout()
        layout.setSpacing(1)
        layout.setContentsMargins(1, 1, 1, 1)

        self.checkButton = QCheckBox()

        layout.addWidget(self.checkButton)
        self.setLayout(layout)

        self.setMinimumHeight(20)
        self.setMaximumHeight(20)

        self._state = state
        self._update()

    @property
    def state(self):
        return self._state

    @state.setter
    def state(self, state):
        self._state = state
        self._update()

    def _update(self):
        self.checkButton.setDisabled(False)

        if self.state == AVAILABLE:
            self.checkButton.setChecked(False)
        elif self.state == CURRENT:
            self.checkButton.setChecked(True)
        elif self.state == OUTDATED:
            self.checkButton.setChecked(True)
        elif self.state == DEPRECATED:
            self.checkButton.setChecked(True)
        elif self.state == USER_FILE:
            self.checkButton.setChecked(False)
            self.checkButton.setDisabled(True)
        else:
            raise ValueError("Invalid state %r" % self.state)

        try:
            self.checkButton.clicked.disconnect()  # Remove old signals if they exist
        except Exception:
            pass

        if not self.checkButton.isChecked():  # Switch signals if the file is present or not
            self.checkButton.clicked.connect(self.download_clicked)
        else:
            self.checkButton.clicked.connect(self.remove_clicked)


class FileState:
    def __init__(self, domain, file_name, info_server, info_local):
        self.domain = domain
        self.filename = file_name

        self.info_server = self.parse_info_file(info_server)
        self.info_local = self.parse_info_file(info_local)

        # self.source = None
        self.state = self.__item_state()

        self.tree_item = None
        self.download_option = None

    def refresh_state(self, info_server=None, info_local=None):
        if info_local:
            self.info_local = self.parse_info_file(info_local)
        if info_server:
            self.info_server = self.parse_info_file(info_server)

        self.state = self.__item_state()

    @property
    def tags(self):
        if self.state in [AVAILABLE]:
            return self.info_server.tags
        elif self.state in [OUTDATED, USER_FILE, CURRENT]:
            return self.info_local.tags
        else:
            return None

    @property
    def title(self):
        if self.state in [AVAILABLE]:
            return self.info_server.title
        elif self.state in [OUTDATED, USER_FILE, CURRENT]:
            return self.info_local.title
        else:
            return None

    @property
    def size(self):
        if self.state in [AVAILABLE]:
            return self.info_server.size
        elif self.state in [OUTDATED, USER_FILE, CURRENT]:
            return self.info_local.size
        else:
            return None

    @property
    def datetime(self):
        if self.state in [USER_FILE, CURRENT, OUTDATED]:
            return self.info_local.datetime
        else:
            return self.info_server.datetime

    @property
    def source(self):
        if self.state is USER_FILE:
            return self.info_local.source.capitalize().replace('_', ' ')
        else:
            return self.info_server.source.capitalize().replace('_', ' ')

    def __item_state(self):
        """ Return item state

        Note:
            available  -> available for download
            current    -> latest version downloaded
            outdated   -> needs update (newer version on serverfiles-bio repository)
            deprecated -> removed from serverfiles-bio repository
            user_file  -> not in serverfiles-bio repository (user defined)

        """

        if not self.info_server and self.info_local:
            # we check source of the file
            if self.info_local.source == SOURCE_USER:
                # this is user defined file
                return USER_FILE

            elif self.info_local.source == SOURCE_SERVER:
                # there is no record of this file on the server
                return DEPRECATED

        if not self.info_local and self.info_server:
            return AVAILABLE

        if self.info_server and self.info_local:
            if not self.info_local.source == SOURCE_USER:

                if self.info_local.datetime < self.info_server.datetime:
                    return OUTDATED
                else:
                    return CURRENT

            else:
                return USER_FILE

    @staticmethod
    def parse_info_file(info):
        """ Parse .info file from JSON like format to namedtuple
        """
        if info is not None:
            if not isinstance(info['datetime'], d_time):
                info['datetime'] = d_time.strptime(info['datetime'], "%Y-%m-%d %H:%M:%S.%f")
            return namedtuple('file_info', info.keys())(**info)


class FileStateItem(QTreeWidgetItem):

    STATE_STRINGS = {
        0: 'not downloaded',
        1: 'downloaded, current',
        2: 'downloaded, needs update',
        3: 'obsolete',
        4: 'custom file',
    }

    #: A role for the state item data.
    StateRole = next(gui.OrangeUserRole)

    # QTreeWidgetItem stores the DisplayRole and EditRole as the same role,
    # so we can't use EditRole to store the actual item data, instead we use
    # custom role.

    #: A custom edit role for the item's data
    #: (QTreeWidget treats Qt.EditRole as a alias for Qt.DisplayRole)
    EditRole2 = next(gui.OrangeUserRole)

    def __init__(self, fs):
        """ A QTreeWidgetItem for displaying a FileState.
        """
        QTreeWidgetItem.__init__(self, type=QTreeWidgetItem.UserType)
        self.update_data(fs)
        self._update_tool_tip(fs)

    def update_data(self, fs):
        self.setData(header.Download, FileStateItem.StateRole, fs.state)

        self.setData(header.Source, Qt.DisplayRole, fs.source)
        self.setData(header.Source, self.EditRole2, fs.source)

        self.setData(header.Title, Qt.DisplayRole, fs.title)
        self.setData(header.Title, self.EditRole2, fs.title)

        if not fs.state == AVAILABLE:
            self.setData(header.Updated, Qt.DisplayRole, fs.datetime.date().isoformat())
            self.setData(header.Updated, self.EditRole2, fs.datetime)
        else:
            self.setData(header.Updated, Qt.DisplayRole, '')
            self.setData(header.Updated, self.EditRole2, '')

        self.setData(header.Size, Qt.DisplayRole, serverfiles.sizeformat(fs.size))
        self.setData(header.Size, self.EditRole2, fs.size)

    def _update_tool_tip(self, fs):
        state_str = self.STATE_STRINGS[fs.state]
        if fs == DEPRECATED:
            diff_date = fs.info_server.datetime - fs.info_local.datetime
        else:
            diff_date = None

        tooltip = "State: {}\nTags: {}".format(state_str, ', '.join(tag for tag in fs.tags if not tag.startswith("#")))

        if fs.state in [CURRENT, OUTDATED, DEPRECATED]:
            tooltip += "\nFile: {}".format(serverfiles.localpath(fs.domain, fs.filename))

        if fs.state == OUTDATED and diff_date:
            tooltip += "\nServer version: {}\nStatus: old {} days".format(fs.datetime, diff_date.days)
        else:
            tooltip += "\nServer version: {}".format(fs.datetime)

        for i in range(1, len(header_labels) - 1):
            self.setToolTip(i, tooltip)

    def __lt__(self, other):
        widget = self.treeWidget()
        column = widget.sortColumn()
        if column == 0:
            role = FileStateItem.StateRole
        else:
            role = self.EditRole2

        left = self.data(column, role)
        right = other.data(column, role)
        try:
            return left < right
        except TypeError:
            pass
        # order lexically by str representation, but ensure `None`
        # always orders on one side
        left = (0, "") if left is None else (1, str(left))
        right = (0, "") if right is None else (1, str(right))
        return left < right


class FileUploadHelper(QDialog):

    # settings
    kegg_domain = 'KEGG'

    supported_domains = OrderedDict({'Gene Ontology': gene_ontology_domain, 'Gene Sets': gene_sets_domain})

    supported_organisms = [common_taxid_to_name(tax_id) for tax_id in common_taxids()]

    hierarchies = {
        'GO - Biological Process': ('GO', 'biological_process'),
        'GO - Molecular Function': ('GO', 'molecular_function'),
        'GO - Cellular Component': ('GO', 'cellular_component'),
        'KEGG - Pathways': ('KEGG', 'pathways'),
        'KEGG - Orthologs': ('KEGG', 'orthologs'),
    }

    def __init__(self, parent=None):
        super(FileUploadHelper, self).__init__(
            parent,
            Qt.Window
            | Qt.WindowTitleHint
            | Qt.CustomizeWindowHint
            | Qt.WindowCloseButtonHint
            | Qt.WindowMaximizeButtonHint,
        )
        self.setAttribute(Qt.WA_DeleteOnClose)
        self.setWindowTitle('Add new file')

        self.info_state = INFO_FILE_SCHEMA
        self.layout = QVBoxLayout(self)

        # domain selection combobox
        self.domain_selection = QComboBox()
        self.domain_selection.addItems(self.supported_domains.keys())
        self.domain_selection.currentIndexChanged.connect(self.__on_domain_selection)
        self.__create_selection_row('Domain: ', self.domain_selection)

        # domain selection combobox
        self.hierarchy_selection = QComboBox()
        self.hierarchy_selection.addItems(self.hierarchies.keys())
        self.layout.addWidget(self.hierarchy_selection, alignment=Qt.AlignVCenter)
        self.__on_domain_selection()

        # select organism
        self.organism_selection = QComboBox()
        self.organism_selection.addItems(self.supported_organisms)
        self.__create_selection_row('Organism: ', self.organism_selection)

        # title
        self.line_edit_title = QLineEdit()
        self.__create_selection_row('Title: ', self.line_edit_title)

        # tags
        self.line_edit_tags = QLineEdit()
        self.__create_selection_row('Tags (comma-separated): ', self.line_edit_tags)

        # file selector
        self.file_info = QLabel()
        self.file_select_btn = QPushButton('Select File', self)
        self.file_select_btn.clicked.connect(self.__handle_file_selector)
        self.__create_selection_row(' ', self.file_select_btn)

        # add file info section
        self.layout.addWidget(self.file_info, alignment=Qt.AlignCenter)

        self.layout.addStretch(1)

        # Ok and Cancel buttons
        self.buttons = QDialogButtonBox(QDialogButtonBox.Ok | QDialogButtonBox.Cancel, Qt.Horizontal, self)
        self.layout.addWidget(self.buttons, alignment=Qt.AlignJustify)

        self.buttons.accepted.connect(self.__accept)
        self.buttons.rejected.connect(self.__close)

        # path to a selected file
        self.file_path = None

    def __on_domain_selection(self):
        selected = self.__get_selected_domain() == gene_sets_domain
        self.hierarchy_selection.setVisible(selected)

    def __get_selected_domain(self):
        domain_label = list(self.supported_domains.keys())[self.domain_selection.currentIndex()]
        return self.supported_domains[domain_label]

    def __get_selected_hier(self):
        hier_label = list(self.hierarchies.keys())[self.hierarchy_selection.currentIndex()]
        return self.hierarchies[hier_label]

    def __create_selection_row(self, label, widget):
        self.layout.addWidget(QLabel(label), alignment=Qt.AlignLeft)
        self.layout.addWidget(widget, alignment=Qt.AlignVCenter)

    def __accept(self):
        if self.file_path:
            self.info_state = self.__parse_selection()
            self.__move_to_serverfiles_folder(self.file_path)

            self.parent().initialize_files_view()
            self.close()

    def __close(self):
        self.close()

    def closeEvent(self, event):
        # clean-up
        self.parent()._dialog = None

    def __filename(self, domain, organism):
        """ Create filename based od domain name and organism.
        """

        if domain in self.supported_domains.values() and domain == gene_ontology_domain and organism:
            return FILENAME_ANNOTATION.format(organism)

        elif domain in self.supported_domains.values() and domain == gene_sets_domain and organism:
            return filename((self.__get_selected_hier()), organism)

    def __parse_selection(self):
        try:
            domain = self.__get_selected_domain()
            organism = species_name_to_taxid(self.supported_organisms[self.organism_selection.currentIndex()])
        except KeyError as e:
            raise e

        return {
            'domain': domain,
            'organism': organism,
            'filename': self.__filename(domain, organism),
            'title': self.line_edit_title.text(),
            'tags': self.line_edit_tags.text().split(','),
            'source': SOURCE_USER,
        }

    def __move_to_serverfiles_folder(self, selected_file_path):
        domain_path = serverfiles.localpath(self.info_state['domain'])
        file_path = os.path.join(domain_path, self.info_state['filename'])
        create_folder(domain_path)

        try:
            copyfile(selected_file_path, file_path)
        except IOError as e:
            # TODO: handle error properly
            raise e

        # if copy successful create .info file
        create_info_file(file_path, **self.info_state)

    def __handle_file_selector(self):
        self.file_path = QFileDialog.getOpenFileName(self, 'Open File')[0]
        self.file_info.setText('Selected File: {}'.format(os.path.basename(self.file_path)))


if __name__ == "__main__":

    def main_test():
        app = QApplication(sys.argv)
        w = OWDatabasesUpdate()
        w.show()
        w.raise_()
        return app.exec_()

    sys.exit(main_test())
