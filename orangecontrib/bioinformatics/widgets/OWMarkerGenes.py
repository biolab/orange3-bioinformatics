""" OWMarkerGenes """
import sys
import numpy as np


from typing import List
from collections import defaultdict
from functools import partial


from AnyQt.QtCore import (
    Qt,
    QSize,
    QTimer,
    QModelIndex,
    QItemSelection,
    QItemSelectionModel,
    QItemSelectionRange,
)
from AnyQt.QtGui import QFont, QColor
from AnyQt.QtWidgets import QTreeView, QLineEdit

from Orange.data import MISSING_VALUES, Table, Domain
from Orange.widgets import widget, gui, settings
from Orange.widgets.utils.itemmodels import TableModel

from orangecontrib.bioinformatics.utils import serverfiles
from orangecontrib.bioinformatics.widgets.utils.data import (
    GENE_AS_ATTRIBUTE_NAME,
    TAX_ID,
    GENE_ID_COLUMN,
)

serverfiles_domain = 'marker_genes'


def get_available_db_sources():
    found_sources = dict()

    try:
        found_sources.update(
            serverfiles.ServerFiles().allinfo(serverfiles_domain)
        )
    except ConnectionError as e:
        raise ConnectionError(
            'Can not connect to {}. Using only local files.'.format(
                serverfiles.server_url
            )
        )
    finally:
        found_sources.update(serverfiles.allinfo(serverfiles_domain))
        return {item.get('title').split(': ')[-1]: item for item in found_sources.values()}


# todo: this is ugly, refactor this in the future.
class HeaderIndex:
    NAME = 0
    GENE = 1
    CELL_TYPE = 2
    FUNCTION = 3
    REFERENCE = 4
    URL = 5


HeaderLabels = {
    HeaderIndex.GENE: 'Entrez ID',
    HeaderIndex.REFERENCE: 'Reference',
    HeaderIndex.URL: 'URL',
}


class SearchableTableModel(TableModel):

    ClassVar, Meta, Attribute = range(3)

    ColorForRole = {ClassVar: None, Meta: None, Attribute: None}

    def __init__(self, data, parent):
        TableModel.__init__(self, data, parent)
        self._data = data
        self._roleData = {Qt.DisplayRole: self.source}
        self._roleData = partial(
            defaultdict,
            partial(defaultdict, partial(defaultdict, lambda: None)),
        )(self._roleData)
        self.set_column_links()

    def set_column_links(self):
        domain = self._data.domain
        ref_col = domain.metas.index(
            domain[HeaderLabels[HeaderIndex.REFERENCE]]
        )
        font = QFont()
        font.setUnderline(True)
        color = QColor(Qt.blue)
        for i, row in enumerate(self.source):
            link = row[HeaderLabels[HeaderIndex.URL]].value
            if len(link) and 'http' in link:
                self._roleData[gui.LinkRole][i][ref_col] = link
                self._roleData[Qt.FontRole][i][ref_col] = font
                self._roleData[Qt.ForegroundRole][i][ref_col] = color

    def data(self, index, role=Qt.DisplayRole):
        if role == Qt.DisplayRole:
            cell_data = super().data(index, role)
            return "" if cell_data in MISSING_VALUES else cell_data
        elif role in (gui.LinkRole, Qt.FontRole, Qt.ForegroundRole):
            row, col = index.row(), index.column()
            return self._roleData[role][row][col]
        return super().data(index, role)

    def rowCount(self, parent=QModelIndex()):
        return 0 if parent.isValid() else len(self.source)

    def wrap(self, table):
        self.beginResetModel()
        self.source = table
        self._roleData.clear()
        self.resetSorting()
        self.endResetModel()

    def filter_table(self, filter_pattern):
        selection = np.full(self._data.metas.shape, True)
        for search_word in filter_pattern.split():
            match_result = (
                np.core.defchararray.find(
                    np.char.lower(self._data.metas.astype(str)),
                    search_word.lower(),
                )
                >= 0
            )
            selection = selection & match_result
        return selection

    def update_model(self, filter_pattern=None):
        # clear cache if model changes
        self._row_instance.cache_clear()
        self.wrap(self._data[self.filter_table(filter_pattern).any(axis=1), :])
        self.set_column_links()


class MarkerGroupContextHandler(settings.ContextHandler):
    def __init__(self):
        super().__init__()

    def match(self, context, group, *args):
        if not context.group == group:
            return self.NO_MATCH

        return self.PERFECT_MATCH

    def new_context(self, group):
        context = super().new_context()
        context.group = group
        return context


class OWMarkerGenes(widget.OWWidget):
    name = "Marker Genes"
    icon = 'icons/OWMarkerGenes.svg'
    priority = 170

    replaces = ['orangecontrib.single_cell.widgets.owmarkergenes.OWMarkerGenes']

    class Outputs:
        genes = widget.Output("Genes", Table)

    want_main_area = True

    selected_group: str = settings.Setting('')
    selected_db_source: str = settings.Setting('')
    filter_text: str = settings.Setting('')
    header_state: bytes = settings.Setting(b'')

    settingsHandler = MarkerGroupContextHandler()
    selected_genes: List[tuple] = settings.ContextSetting([])

    def __init__(self):
        super().__init__()
        self._data = None
        self._available_db_sources = None

        self._timer = QTimer()
        self._timer.timeout.connect(self._filter_table)
        self._timer.setSingleShot(True)

        box = gui.widgetBox(self.controlArea, 'Database', margin=0)
        self.db_source_index = -1
        self.db_source_cb = gui.comboBox(box, self, 'db_source_index')
        self.db_source_cb.activated[int].connect(self.handle_source_changed)

        box = gui.widgetBox(self.controlArea, 'Organism', margin=0)
        self.group_index = -1
        self.group_cb = gui.comboBox(box, self, 'group_index')
        self.group_cb.activated[int].connect(self.set_group_index)
        gui.rubber(self.controlArea)

        # TODO: to avoid this, marker genes table should have 'tax_id' column
        self.map_group_to_taxid = {'Human': '9606', 'Mouse': '10090'}

        filter_line_edit = gui.lineEdit(
            self.mainArea, self, "filter_text"
        )  # type: QLineEdit
        filter_line_edit.setPlaceholderText("Filter...")
        filter_line_edit.textEdited.connect(self.call_filter_timer)

        self.view = view = QTreeView(
            rootIsDecorated=False,
            uniformRowHeights=True,
            selectionMode=QTreeView.ExtendedSelection,
            sortingEnabled=True,
        )

        view.viewport().setMouseTracking(True)
        self.mainArea.layout().addWidget(view)

        self._load_data()
        if self.header_state:
            view.header().restoreState(self.header_state)

    @property
    def available_db_sources(self) -> dict:
        return self._available_db_sources

    @available_db_sources.setter
    def available_db_sources(self, value: dict):
        self._available_db_sources = value

        items = list(value.keys())
        try:
            idx = items.index(self.selected_db_source)
        except ValueError:
            idx = -1

        self.db_source_cb.clear()
        self.db_source_cb.addItems(items)

        if idx != -1:
            self.db_source_index = idx
            self.selected_db_source = items[idx]
        elif items:
            self.db_source_index = min(
                max(self.db_source_index, 0), len(items) - 1
            )

        self.set_db_source_index(self.db_source_index)

    @property
    def data(self) -> Table:
        return self._data

    @data.setter
    def data(self, value: Table):
        """ Set the source data.

        The data is then filtered on the first meta column (group)
        """
        self._data = value
        domain = value.domain

        if domain.metas:
            group = domain.metas[0]
            groupcol, _ = value.get_column_view(group)

            if group.is_string:
                group_values = list(set(groupcol))
            elif group.is_discrete:
                group_values = group.values
            else:
                raise TypeError("Invalid column type")
            try:
                idx = group_values.index(self.selected_group)
            except ValueError:
                idx = -1

            self.group_cb.clear()
            self.group_cb.addItems(group_values)

            if idx != -1:
                self.group_index = idx
                self.selected_group = group_values[idx]
            elif group_values:
                self.group_index = min(
                    max(self.group_index, 0), len(group_values) - 1
                )

            self.set_group_index(self.group_index)

    def _load_data(self):
        self.available_db_sources = get_available_db_sources()
        file_name = self.available_db_sources[self.selected_db_source][
            'filename'
        ]

        try:
            serverfiles.update(serverfiles_domain, file_name)
        except ConnectionError as e:
            raise ConnectionError(
                'Can not connect to {}. '
                'Using only local files.'.format(serverfiles.server_url)
            )
        finally:
            file_path = serverfiles.localpath_download(
                serverfiles_domain, file_name
            )
            data = Table(file_path)

            # enforce order
            old_domain = data.domain
            new_domain = Domain(
                [],
                metas=[
                    old_domain['Organism'],
                    old_domain['Name'],
                    old_domain['Entrez ID'],
                    old_domain['Cell Type'],
                    old_domain['Function'],
                    old_domain['Reference'],
                    old_domain['URL'],
                ],
            )
            data = data.transform(new_domain)
            self.data = data

    def set_selection(self):
        selected = self.selected_rows()

        if len(selected):
            header_count = self.view.header().count() - 1

            if self.view.model().rowCount() <= selected[-1]:
                return

            selection = QItemSelection()

            for row_index in selected:
                selection.append(
                    QItemSelectionRange(
                        self.view.model().index(row_index, 0),
                        self.view.model().index(row_index, header_count),
                    )
                )

            self.view.selectionModel().select(
                selection, QItemSelectionModel.ClearAndSelect
            )

    def handle_source_changed(self, source_index):
        self.set_db_source_index(source_index)
        self._load_data()

    def set_db_source_index(self, source_index):
        self.closeContext()
        self.db_source_index = source_index
        self.selected_db_source = self.db_source_cb.itemText(source_index)

    def set_group_index(self, group_index):
        self.closeContext()
        self.group_index = group_index
        self.selected_group = self.group_cb.itemText(group_index)
        self._setup()

    def call_filter_timer(self, search_string):
        self._timer.stop()
        if search_string != self.filter_text:
            self.filter_text = search_string
        self._timer.start(500)

    def _filter_table(self):
        model = self.view.model()
        assert isinstance(model, SearchableTableModel)
        model.update_model(str(self.filter_text))

    def _setup(self):
        self.closeContext()
        data = self.data
        group = data.domain.metas[0]
        gvec = data.get_column_view(group)[0]

        if group.is_string:
            mask = gvec == self.group_cb.itemData(
                self.group_index, Qt.DisplayRole
            )
        else:
            mask = gvec == self.group_index

        data = data[mask]
        rest = data[:, data.domain.metas[1:]]
        model = SearchableTableModel(rest, parent=self)
        ref_col = rest.domain.metas.index(
            rest.domain[HeaderLabels[HeaderIndex.REFERENCE]]
        )
        self.view.setItemDelegateForColumn(
            ref_col, gui.LinkStyledItemDelegate(self.view)
        )

        self.view.setModel(model)
        self.view.selectionModel().selectionChanged.connect(
            self._on_selection_changed
        )

        self.openContext(self.selected_group)
        self.call_filter_timer(self.filter_text)
        self.view.hideColumn(HeaderIndex.URL)
        self.set_selection()

        self.commit()

    def _on_selection_changed(self):
        self.commit()

    def selected_rows(self):
        """ Return row index for selected genes
        """
        if not self.selected_genes:
            return []

        model = self.view.model()
        return [
            row_index
            for row_index in range(model.rowCount())
            if (
                model.index(row_index, HeaderIndex.GENE).data(),
                model.index(row_index, HeaderIndex.CELL_TYPE).data(),
                model.index(row_index, HeaderIndex.REFERENCE).data(),
            )
            in self.selected_genes
        ]

    def commit(self):
        model = self.view.model()
        assert isinstance(model, SearchableTableModel)
        rows = [mi.row() for mi in self.view.selectionModel().selectedRows(0)]

        if rows:
            rows = model.mapToSourceRows(rows)
            output = model.source[rows]
        else:
            output = model.source

        gene_id = self.view.selectionModel().selectedRows(HeaderIndex.GENE)
        cell_type = self.view.selectionModel().selectedRows(
            HeaderIndex.CELL_TYPE
        )
        ref = self.view.selectionModel().selectedRows(HeaderIndex.REFERENCE)

        self.selected_genes = [
            (entrez.data(), cell.data(), ref.data())
            for entrez, cell, ref in zip(gene_id, cell_type, ref)
        ]

        # always false for marker genes data tables in single cell
        output.attributes[GENE_AS_ATTRIBUTE_NAME] = False
        # set taxonomy id in data.attributes
        output.attributes[TAX_ID] = self.map_group_to_taxid.get(
            self.selected_group, ''
        )
        # set column id flag
        output.attributes[GENE_ID_COLUMN] = HeaderLabels[HeaderIndex.GENE]
        output.name = 'Marker Genes'

        self.Outputs.genes.send(output)

    def closeEvent(self, event):
        self.header_state = bytes(self.view.header().saveState())
        super().closeEvent(event)

    def sizeHint(self):
        return super().sizeHint().expandedTo(QSize(600, 500))


def main(argv=None):
    from AnyQt.QtWidgets import QApplication

    app = QApplication(argv or sys.argv)
    w = OWMarkerGenes()
    w.show()
    w.activateWindow()
    rv = app.exec_()
    w.saveSettings()
    w.onDeleteWidget()
    return rv


if __name__ == "__main__":
    sys.exit(main())
