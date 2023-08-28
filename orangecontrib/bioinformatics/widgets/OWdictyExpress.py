""" dictyExpress widget """
from typing import Optional

from requests.exceptions import ConnectionError

from AnyQt.QtGui import QFont
from AnyQt.QtCore import Qt, QSize
from AnyQt.QtWidgets import QLabel, QTreeWidget, QTreeWidgetItem

from Orange.data import Table, StringVariable
from Orange.widgets import gui, settings
from Orange.widgets.widget import Msg, OWWidget
from Orange.widgets.utils.signals import Output
from Orange.widgets.utils.concurrent import ConcurrentWidgetMixin

from orangecontrib.bioinformatics import resolwe
from orangecontrib.bioinformatics.resolwe import genapi, connect
from orangecontrib.bioinformatics.ncbi.gene import ENTREZ_ID, GeneMatcher
from orangecontrib.bioinformatics.resolwe.utils import etc_to_table
from orangecontrib.bioinformatics.widgets.utils.data import (
    TAX_ID,
    GENE_ID_COLUMN,
    GENE_ID_ATTRIBUTE,
    GENE_AS_ATTRIBUTE_NAME,
)
from orangecontrib.bioinformatics.widgets.components.resolwe import (
    SignIn,
    get_credential_manager,
)

Labels = [
    (" ", " "),
    ("var.project", "Project"),
    ("static.name", "Experiment"),
    ("static.cite", "Citation"),
    ("var.growth", "Growth"),
    ("var.treatment", "Treatment"),
    ("var.strain", "Strain"),
]


class OWdictyExpress(OWWidget, ConcurrentWidgetMixin):
    name = "dictyExpress"
    description = "Time-course gene expression data"
    icon = "../widgets/icons/OWdictyExpress.svg"
    want_main_area = True
    priority = 20

    class Inputs:
        pass

    class Outputs:
        etc_data = Output("Data", Table)

    class Error(OWWidget.Error):
        unreachable_host = Msg('Host not reachable')
        invalid_credentials = Msg('Invalid credentials')

    gene_as_attr_name = settings.Setting(0)

    selected_item = settings.Setting(None, schema_only=True)
    auto_commit = settings.Setting(False, schema_only=True)

    def __init__(self):
        super().__init__()
        ConcurrentWidgetMixin.__init__(self)

        self._res: Optional[genapi.GenAPI] = None
        self.organism = '44689'
        self.server = 'https://dictyexpress.research.bcm.edu'
        self.headerLabels = [x[1] for x in Labels]
        self.searchString = ""
        self.items = []
        self.genapi_pub_auth = {
            'url': genapi.DEFAULT_URL,
            'username': genapi.DEFAULT_EMAIL,
            'password': genapi.DEFAULT_PASSWD,
        }

        # Login Section
        box = gui.widgetBox(self.controlArea, 'Sign in')
        self.user_info = gui.label(box, self, '')
        self.server_info = gui.label(box, self, '')

        box = gui.widgetBox(box, orientation=Qt.Horizontal)
        self.sign_in_btn = gui.button(
            box, self, 'Sign In', callback=self.sign_in, autoDefault=False
        )
        self.sign_out_btn = gui.button(
            box, self, 'Sign Out', callback=self.sign_out, autoDefault=False
        )

        box = gui.widgetBox(self.controlArea, "Output")
        gui.radioButtonsInBox(
            box,
            self,
            "gene_as_attr_name",
            ["Genes in rows", "Genes in columns"],
            callback=self.invalidate,
        )

        self.clear_cache_btn = gui.button(
            self.controlArea,
            self,
            "Clear cache",
            autoDefault=False,
            callback=self.clear_cache,
        )

        gui.rubber(self.controlArea)

        self.commit_button = gui.auto_commit(
            self.controlArea, self, "auto_commit", "&Commit", box=False
        )

        # Experiment Section

        label = QLabel("Available projects:")
        my_font = QFont()
        my_font.setBold(True)
        label.setFont(my_font)
        self.mainArea.layout().addWidget(label)

        self.filter = gui.lineEdit(
            self.mainArea,
            self,
            "searchString",
            "Filter:",
            callbackOnType=True,
            callback=self.search_update,
        )

        self.experimentsWidget = QTreeWidget(
            alternatingRowColors=True,
            rootIsDecorated=False,
            uniformRowHeights=True,
            sortingEnabled=True,
        )

        self.experimentsWidget.setItemDelegateForColumn(
            0, gui.IndicatorItemDelegate(self, role=Qt.DisplayRole)
        )

        self.experimentsWidget.selectionModel().selectionChanged.connect(
            self.on_selection_changed
        )

        self.experimentsWidget.setHeaderLabels(self.headerLabels)
        self.mainArea.layout().addWidget(self.experimentsWidget)

        self.sign_in(silent=True)
        self.sizeHint()

    def sizeHint(self):
        return QSize(1400, 680)

    @property
    def res(self):
        return self._res

    @res.setter
    def res(self, value: genapi.GenAPI):
        if isinstance(value, genapi.GenAPI):
            self._res = value
            self.Error.clear()
            self.reset()
            self.load_experiments()
            self.update_user_status()
            self.Outputs.etc_data.send(None)

    def sign_in(self, silent=False):
        dialog = SignIn(self, server_type='genesis')

        if silent:
            dialog.sign_in()
            if dialog.resolwe_instance is not None:
                self.res = dialog.resolwe_instance
            else:
                self.res = connect(
                    **self.genapi_pub_auth, server_type=resolwe.GENESIS_PLATFORM
                )

        if not silent and dialog.exec_():
            self.res = dialog.resolwe_instance

    def sign_out(self):
        # Remove username and password
        cm = get_credential_manager(resolwe.GENESIS_PLATFORM)
        del cm.username
        del cm.password
        # Use public credentials when user signs out
        self.res = connect(**self.genapi_pub_auth, server_type=resolwe.GENESIS_PLATFORM)

    def update_user_status(self):
        cm = get_credential_manager(resolwe.GENESIS_PLATFORM)

        if cm.username:
            user_info = f"User: {cm.username}"
            self.sign_in_btn.setEnabled(False)
            self.sign_out_btn.setEnabled(True)
        else:
            user_info = 'User: Anonymous'
            self.sign_in_btn.setEnabled(True)
            self.sign_out_btn.setEnabled(False)

        self.user_info.setText(user_info)
        self.server_info.setText(f'Server: {self.res._gen.url[8:]}')

    def clear_cache(self):
        resolwe.GenAPI.clear_cache()
        self.reset()
        self.load_experiments()

    def reset(self):
        self.experimentsWidget.clear()  # clear QTreeWidget
        self.items = []
        # self.lastSelected = None
        self.searchString = ""

    def search_update(self):
        parts = self.searchString.split()
        for item in self.items:
            item.setHidden(not all(s in item for s in parts))

    def on_exception(self, ex):
        if isinstance(ex, ConnectionError) or isinstance(ex, ValueError):
            self.Error.unreachable_host()

        print(ex)

    def on_done(self, results):
        if isinstance(results, list):
            self.load_tree_items(results)
        elif isinstance(results, tuple):
            self.send_to_output(results)

    def load_experiments(self):
        if self.res:
            self.start(self.res.fetch_etc_objects)

    def load_tree_items(self, list_of_exp):
        self.items = [
            CustomTreeItem(self.experimentsWidget, item) for item in list_of_exp
        ]

        for i in range(len(self.headerLabels)):
            self.experimentsWidget.resizeColumnToContents(i)

        self.set_cached_indicator()
        self.set_selected()

    def set_selected(self):
        for item in self.items:
            if self.selected_item and item.gen_data_id == self.selected_item:
                self.experimentsWidget.setCurrentItem(item)

    def on_selection_changed(self):
        self.invalidate()

    def invalidate(self):
        self.commit()

    def send_to_output(self, result):
        etc_json, table_name = result

        # convert to table
        data = etc_to_table(etc_json, bool(self.gene_as_attr_name))
        # set table name
        data.name = table_name

        # match genes
        gene_matcher = GeneMatcher(str(self.organism))

        if not bool(self.gene_as_attr_name):
            if 'Gene' in data.domain:
                data = gene_matcher.match_table_column(
                    data, 'Gene', StringVariable(ENTREZ_ID)
                )
            data.attributes[GENE_ID_COLUMN] = ENTREZ_ID
        else:
            data = gene_matcher.match_table_attributes(data)
            data.attributes[GENE_ID_ATTRIBUTE] = ENTREZ_ID

        # add table attributes
        data.attributes[TAX_ID] = str(self.organism)
        data.attributes[GENE_AS_ATTRIBUTE_NAME] = bool(self.gene_as_attr_name)

        # reset cache indicators
        self.set_cached_indicator()
        # send data to the output signal
        self.Outputs.etc_data.send(data)

    def commit(self):
        self.Error.clear()
        selected_items = self.experimentsWidget.selectedItems()  # get selected TreeItem

        if len(selected_items) < 1:
            self.Outputs.etc_data.send(None)
            return

        selected_item = selected_items[0]
        self.selected_item = selected_item.gen_data_id
        self.start(
            self.res.download_etc_data,
            selected_item.gen_data_id,
            table_name=selected_item.data_name,
        )

    def set_cached_indicator(self):
        cached = self.res.get_cached_ids()

        for item in self.items:
            if item.gen_data_id in cached:
                item.setData(0, Qt.DisplayRole, " ")
            else:
                item.setData(0, Qt.DisplayRole, "")


class CustomTreeItem(QTreeWidgetItem):
    def __init__(self, parent, gen_data):
        super(CustomTreeItem, self).__init__(
            parent
        )  # Init super class (QtGui.QTreeWidgetItem )

        self._gen_data = gen_data  # GenData object
        self.set_rows(self._gen_data.annotation)  # set rows in QTreeWidget

    def __contains__(self, text):
        return any(
            text.upper() in str(self.text(i)).upper() for i in range(self.columnCount())
        )

    @property
    def gen_data_id(self):
        return self._gen_data.id

    @property
    def data_name(self):
        try:
            project = self._gen_data.var['project']
            experiment = self._gen_data.static['name']
        except (AttributeError, KeyError):
            project = ''
            experiment = ''

        return '{} ({})'.format(project, experiment)

    def set_rows(self, row):
        for index, label in enumerate(Labels):
            if index > 0:
                try:
                    if isinstance(row[label[0]]["value"], list):
                        self.setText(index, row[label[0]]["value"][0]["name"])
                    else:
                        self.setText(index, row[label[0]]["value"])
                except (IndexError, KeyError):
                    self.setText(index, 'No data')


if __name__ == "__main__":
    from orangewidget.utils.widgetpreview import WidgetPreview

    WidgetPreview(OWdictyExpress).run()
