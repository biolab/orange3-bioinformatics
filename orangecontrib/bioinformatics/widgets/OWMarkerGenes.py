""" OWMarkerGenes """
from typing import List, Tuple, Iterable, Optional
from functools import partial

import numpy as np
import requests
from orangewidget.settings import Setting, ContextHandler

from AnyQt import QtGui, QtCore
from AnyQt.QtCore import Qt, QSize, QMimeData, QModelIndex, QAbstractItemModel, QSortFilterProxyModel, pyqtSignal
from AnyQt.QtWidgets import QLayout, QWidget, QLineEdit, QTreeView, QGridLayout, QTextBrowser, QAbstractItemView

from Orange.data import Table, RowInstance
from Orange.widgets import gui, widget, settings

from orangecontrib.bioinformatics.utils import serverfiles
from orangecontrib.bioinformatics.widgets.utils.data import TAX_ID, GENE_ID_COLUMN, GENE_AS_ATTRIBUTE_NAME
from orangecontrib.bioinformatics.widgets.utils.settings import MarkerGroupContextHandler

SERVER_FILES_DOMAIN = 'marker_genes'
GROUP_BY_ITEMS = ["Cell Type", "Function"]
FILTER_COLUMNS_DEFAULT = ["Name", "Entrez ID"]
NUM_LINES_TEXT = 5
MAP_GROUP_TO_TAX_ID = {'Human': '9606', 'Mouse': '10090'}


class TreeItem(object):
    def __init__(
        self, name: str, is_type_item: bool, data_row: Optional[RowInstance], parent: Optional["TreeItem"] = None
    ):
        self.name: str = name
        self.parentItem: Optional["TreeItem"] = None
        self.childItems: List = []
        self.change_parent(parent)
        self.isTypeItem: bool = is_type_item

        self.data_row: Optional[RowInstance] = data_row
        # holds the information whether element is expanded in the tree
        self.expanded = False

    def change_parent(self, parent: Optional["TreeItem"] = None) -> None:
        """
        Function replaces the item parent.

        Parameters
        ----------
        parent
            Parent which is set as an items parent
        """
        self.parentItem = parent
        if parent is not None:
            parent.append_child(self)

    @staticmethod
    def change_parents(items: Iterable["TreeItem"], parent: "TreeItem") -> None:
        """
        This function assigns all items a parent

        Parameters
        ----------
        items
            Items that will get parent assigned
        parent
            Parent that will be assigned
        """
        for it in items:
            it.change_parent(parent)

    def append_child(self, item: "TreeItem") -> None:
        """
        This function appends child to self. This function just append a child
        and do not set parent of child item (it is set by function calling this function.

        Parameters
        ----------
        item
            Child item.
        """
        self.childItems.append(item)

    def child(self, row: int) -> "TreeItem":
        """
        Returns row-th child the item.

        Parameters
        ----------
        row
            Index of the child

        Returns
        -------
        Child item
        """
        return self.childItems[row]

    def remove_children(self, position: int, count: int) -> bool:
        """
        This function removes children form position to position + count

        Parameters
        ----------
        position
            Position of the first removed child
        count
            Number of removed children
        Returns
        -------
        Success of removing
        """
        # we could check for every index separately but in case when all items cannot be removed there is something
        # wrong with the cal - just refuse to do the removal
        if position < 0 or position + count > len(self.childItems):
            return False
        # remove from the back otherwise you get index out of range when removing more than one
        for i in range(position + count - 1, position - 1, -1):
            self.childItems.pop(i)
        return True

    def children_count(self) -> int:
        """
        Function return number fo children of an item.

        Returns
        -------
        Number of children.
        """
        return len(self.childItems)

    def row(self) -> int:
        """
        Function returns the index of the item in a parents list. If it does not have parent returns 0.

        Returns
        -------
        The index of an item.
        """
        if self.parentItem:
            return self.parentItem.childItems.index(self)
        return 0

    def child_from_name_index(self, name: str) -> Tuple[Optional[int], Optional["TreeItem"]]:
        """
        This function returns the index (in the array of children) and item with the name specified.

        Parameters
        ----------
        name
            The name of requested item.

        Returns
        -------
        The index and the searched item. If it does not exist return None.
        """
        for i, c in enumerate(self.childItems):
            if c.name == name:
                return i, c
        return None, None

    def contains_text(self, text: str, filter_columns: List[str]) -> bool:
        """
        Function indicates whether the text is present in any of columns in FILTER_COLUMNS.

        Parameters
        ----------
        text
            Tested text.
        filter_columns
            Columns used in filtering

        Returns
        -------
        Boolean that indicates that text is present in the item.
        """
        if self.data_row is not None:
            return any(text.lower() in str(self.data_row[col]).lower() for col in filter_columns)
        else:
            return False

    def get_data_rows(self) -> List[RowInstance]:
        """
        This function returns a list of rows present in self and any children node.

        Returns
        -------
        Rows with data items.
        """
        rows = [self.data_row] if self.data_row is not None else []
        for c in self.childItems:
            rows += c.get_data_rows()
        return rows

    def count_marker_genes(self) -> int:
        """
        This function counts number of marker genes in its tree

        Returns
        -------
        int
            Number of marker genes in the tree including itself.
        """
        return sum(c.count_marker_genes() for c in self.childItems) + (0 if self.data_row is None else 1)


class TreeModel(QAbstractItemModel):

    MIME_TYPE = "application/x-Orange-TreeModelData"
    # it was slow to relay on the qt model signals since they are emitted every time endInsertRows is called
    # this one will be emitted just a the end
    data_added = pyqtSignal()
    data_removed = pyqtSignal()
    expand_item = pyqtSignal(QModelIndex, bool)

    def __init__(self, data: Table, parent_column: str, parent: QModelIndex = None) -> None:
        super().__init__(parent)

        self.rootItem = TreeItem("Genes", False, None)
        self.setup_model_data(data, self.rootItem, parent_column)
        self.filter_columns = FILTER_COLUMNS_DEFAULT + [parent_column]

    def flags(self, index: QModelIndex) -> Qt.ItemFlag:
        if not index.isValid():
            return Qt.ItemIsEnabled | Qt.ItemIsSelectable | Qt.ItemIsDropEnabled
        return Qt.ItemIsEnabled | Qt.ItemIsSelectable | Qt.ItemIsDragEnabled | Qt.ItemIsDragEnabled

    def supportedDropActions(self) -> Qt.DropAction:  # noqa: N802
        return Qt.MoveAction

    def supportedDragActions(self) -> Qt.DropAction:  # noqa: N802
        return Qt.MoveAction

    def columnCount(self, parent: QModelIndex = None) -> int:  # noqa: N802
        """
        Our TreeView has only one column
        """
        return 1

    def rowCount(self, parent: QModelIndex = None) -> int:  # noqa: N802
        """
        Function returns the number children of an item.
        """
        if not parent or not parent.isValid():
            parent_item = self.rootItem
        else:
            parent_item = self.node_from_index(parent)

        count = parent_item.children_count()
        return count

    def data(self, index: QModelIndex, role=None) -> Optional[str]:
        """
        Function returns the index's data that are shown in the view.
        """
        if not index.isValid() or role != Qt.DisplayRole:
            return None
        item = self.node_from_index(index)
        rd = item.name
        return rd

    def headerData(self, section, orientation, role=None):  # noqa: N802
        """
        Return header data (they are no used so function is just here to make Qt happy).
        """
        if orientation == Qt.Horizontal and role == Qt.DisplayRole:
            return self.rootItem.name
        return None

    def index(self, row: int, column: int, parent: QModelIndex = None) -> QModelIndex:
        """
        Function create QModelIndex from the parent item and row and column index.

        Parameters
        ----------
        row
            Row of the child
        column
            Column of th child - not used since we have only one column.
        parent
            Parent of the child we are making index of.

        Returns
        -------
        Index of the child at row-th row of the parent.
        """
        if not self.hasIndex(row, column, parent):
            return QModelIndex()

        if not parent.isValid():
            parent_item = self.rootItem
        else:
            parent_item = self.node_from_index(parent)

        child_item = parent_item.child(row)
        if child_item:
            return self.createIndex(row, column, child_item)
        else:
            return QModelIndex()

    def parent(self, index: QModelIndex) -> QModelIndex:
        """
        Function returns parent of the item in index.
        """
        if not index.isValid():
            return QModelIndex()

        child_item = self.node_from_index(index)
        parent_item = child_item.parentItem if child_item else None

        if not parent_item or parent_item == self.rootItem:
            return QModelIndex()

        return self.createIndex(parent_item.row(), 0, parent_item)

    def removeRows(  # noqa: N802
        self, position: int, rows: int, parent: QModelIndex = QModelIndex(), emit_data_removed: bool = True
    ) -> bool:
        """
        This function implements removing rows form position to position + rows form the model. It also removes
        a parent when it becomes empty. If remove is called over the parent it remove children together with the
        parent.

        Parameters
        ----------
        emit_data_removed
            It indicates whether we emit dataRemoved signal. It makes possible that it is emitted when function called
            by Qt and not when called by removeNodeList. When qt call this function indices are optimized to be in
            range and this function is called just once for consecutive data. It is not true for data from
            removeNodeList. In this case we call removeRows for each removed item and it can become slow when many
            items are removed.
        """
        node = self.node_from_index(parent)
        if not node:
            return False
        if not (position == 0 and node.children_count() == rows) or not node.isTypeItem:
            # when not all children are removed or when type item is removed - when removed node is a type item its
            # parent is not
            self.beginRemoveRows(QModelIndex() if node is self.rootItem else parent, position, position + rows - 1)
            success = node.remove_children(position, rows)
        else:
            # when node is a children of a type item and all children of node needs to be removed, remove node itself
            self.beginRemoveRows(QModelIndex(), node.row(), node.row())
            success = self.rootItem.remove_children(node.row(), 1)
        self.endRemoveRows()
        if emit_data_removed:
            self.data_removed.emit()
        return success

    def mimeTypes(self) -> List[str]:  # noqa: N802
        return [self.MIME_TYPE]

    def mimeData(self, index_list: List[QModelIndex]) -> QMimeData:  # noqa: N802
        """
        This function packs data in QMimeData in order to move them to other view with drag and drop functionality.
        """
        items = [self.node_from_index(index) for index in index_list]
        mime = QMimeData()
        # the encoded 'data' is empty, variables are passed by properties
        mime.setData(self.MIME_TYPE, b'')
        mime.setProperty("_items", items)
        mime.setProperty("_source_id", id(self))
        return mime

    def dropMimeData(  # noqa: N802
        self, mime: QMimeData, action: Qt.DropAction, row: int, column: int, parent: QModelIndex = QModelIndex
    ) -> bool:
        """
        This function unpacks QMimeData and insert them to the view where they are dropped.
        Inserting is a bit complicated since sometimes we need to create a parent when it is not in the view yet:
        - when child (marker gene) is dropped we insert it under the current parent (group) if it exist othervise
          we crate a new parent.
        - when group item is inserted it merge it with current same parent item if exist or insert it if it does
          not exists
        """
        if action == Qt.IgnoreAction:
            return True
        if not mime.hasFormat(self.MIME_TYPE):
            return False
        items = mime.property("_items")
        if items is None:
            return False
        # when parent item in items remove its child since all group will be moved
        items = [it for it in items if it.parentItem not in items]
        self.insert_items(items)
        self.data_added.emit()
        return True

    def insert_items(self, items: List[TreeItem]) -> None:
        """
        This function goes through all items and insert them in the view.
        Parameters
        ----------
        items
            List of items that will be inserted.
        """
        for item in items:
            if item.isTypeItem:
                # inserting group item
                self.insert_group(item)
            else:
                # inserting children (genes)
                self.insert_child_item(item)

    def insert_group(self, item: TreeItem) -> None:
        """
        This function makes insert of the item in the model in case that it is a group item.

        Parameters
        ----------
        item
            Item that is inserted - group item (not gene).
        """
        expand = None

        item_curr_idx, curr_item = self.rootItem.child_from_name_index(item.name)
        if item_curr_idx is None:
            # it is not in the view yet - just insert it
            self.beginInsertRows(QModelIndex(), self.rootItem.children_count(), self.rootItem.children_count())
            item.change_parent(self.rootItem)
            expand = item.expanded
        else:
            # it is in the view already - insert only children
            tr_idx = self.createIndex(item_curr_idx, 0, curr_item)
            self.beginInsertRows(tr_idx, self.rowCount(tr_idx), self.rowCount(tr_idx) + len(item.childItems) - 1)
            TreeItem.change_parents(item.childItems, curr_item)
        self.endInsertRows()

        # emit the signal if item must be expanded
        if expand is not None:
            self.expand_item.emit(self.index_from_node(item), expand)

    def insert_child_item(self, item: TreeItem) -> None:
        """
        This function makes insert of the item in the model in case that it is a child (gene) item.

        Parameters
        ----------
        item
            Item that is inserted - child item (gene).
        """
        expand_parent = None
        parent_curr_idx, curr_parent = self.rootItem.child_from_name_index(item.parentItem.name)
        if parent_curr_idx is None:
            # parent (group) do not exists yet - crate it
            self.beginInsertRows(QModelIndex(), self.rootItem.children_count(), self.rootItem.children_count())
            curr_parent = TreeItem(item.parentItem.name, True, None, self.rootItem)
            # when creating new parent because child is inserted always expend since it is expanded in the other view
            expand_parent = True
        else:
            # parent (group) exists - insert into it
            tr_idx = self.createIndex(parent_curr_idx, 0, curr_parent)
            self.beginInsertRows(tr_idx, self.rowCount(tr_idx), self.rowCount(tr_idx))
        item.change_parent(curr_parent)
        self.endInsertRows()

        # emit the signal if parent must be expanded
        if expand_parent is not None:
            self.expand_item.emit(self.index_from_node(curr_parent), expand_parent)

    def canDropMimeData(  # noqa: N802
        self, data: QMimeData, action: Qt.DropAction, row: int, column: int, parent: QModelIndex
    ) -> bool:
        """
        This function enable or disable drop action to the view. With current implementation we disable drop action to
        the same view.
        """
        if data.property("_source_id") == id(self):
            return False
        else:
            return True

    def remove_node_list(self, indices: List[QModelIndex]) -> None:
        """
        Sometimes we need to remove items provided as a list of QModelIndex. Since it is not possible directly with
        removeRows it is implemented with this function.

        Parameters
        ----------
        indices
            List of nodes as QModelIndex
        """
        for idx in indices:
            node = self.node_from_index(idx)
            if node in node.parentItem.childItems:  # it is possible that item already removed since no children left
                parent_idx = self.index_from_node(node.parentItem)
                self.removeRows(node.row(), 1, parent_idx, emit_data_removed=False)
        self.data_removed.emit()

    def node_from_index(self, index: QModelIndex) -> TreeItem:  # noqa: N802
        """
        Function returns a TreeItem saved in the index.

        Parameters
        ----------
        index
            Index of the TreeItem

        Returns
        -------
            Tree item from the index.
        """
        if index.isValid():
            return index.internalPointer()
        else:
            return self.rootItem

    def index_from_node(self, node: TreeItem) -> QModelIndex:  # noqa: N802
        """
        This function is similar to the index function. We use it in case when we create index but not base on the
        parent index. We could call index but it is faster to create it directly.

        Parameters
        ----------
        node
            TreeItem for which we crate index.

        Returns
        -------
        QModelIndex for TreeItem
        """
        return self.createIndex(node.row(), 0, node)

    def index_from_name(self, name: str, parent: TreeItem) -> Optional[QModelIndex]:  # noqa: N802
        """
        Function create index from tree item name.

        Parameters
        ----------
        name
            Name of the item
        parent
            Parent TreeItem

        Returns
        -------
        QModelIndex for an item
        """
        i, node = parent.child_from_name_index(name)
        if i is None:
            return None
        return self.createIndex(i, 0, node)

    def setup_model_data(self, data_table: Table, parent: TreeItem, parent_column: str) -> None:  # noqa: N802
        """
        Function populates the view with the items from the data table. Items are inserted as a tree with the depth 2
        each row (gene from the table) get inserted under the parent item (group) which is defined with the
        paren_column parameter (e.g. marker gene).

        Parameters
        ----------
        data_table
            Data table with marker genes.
        parent
            Parent under which items are inserted (usually root item of the model).
        parent_column
            Column which is the group for items in the tree.
        """
        parents_dict = {}
        names = data_table.get_column_view("Name")[0]
        types = data_table.get_column_view(parent_column)[0]
        for n, pt, row in zip(names, types, data_table):
            if pt not in parents_dict:
                parents_dict[pt] = TreeItem(pt, True, None, parent)

            TreeItem(n, False, row, parents_dict[pt])

    def set_expanded(self, index: QModelIndex, expanded: bool) -> None:
        """
        Sets the expanded flag of the item. This flag tell whether the item in the view is expanded.

        Parameters
        ----------
        index
            Index that is expanded/collapsed
        expanded
            If true expanded otherwise collapsed
        """
        node = self.node_from_index(index)
        node.expanded = expanded

    def __len__(self):
        """
        Len function in this case returns number of marker genes in the view.

        Returns
        -------
        int
            Number of marker genes in the view.
        """
        return self.rootItem.count_marker_genes()


class FilterProxyModel(QSortFilterProxyModel):
    """
    FilterProxyModel is used to sort items in alphabetical order, and to make filter on the view of available genes.

    Parameters
    ----------
    filter_line_edit
        Line edit used to filter the view.
    """

    def __init__(self, filter_line_edit: QLineEdit) -> None:
        super().__init__()
        self.filter_line_edit = filter_line_edit
        if filter_line_edit is not None:
            filter_line_edit.textChanged.connect(self.setFilterFixedString)

        # have marker genes in the alphabetical order
        self.sort(0)

    def filterAcceptsRow(self, p_int: int, index: QModelIndex) -> bool:  # noqa: N802
        """
        Indicates whether the item should be shown based on the the filter string.
        """
        if self.filter_line_edit is None:
            return True
        model = self.sourceModel()
        idx = model.index(p_int, 0, index)
        res = model.node_from_index(idx).contains_text(self.filter_line_edit.text(), self.sourceModel().filter_columns)

        if model.hasChildren(idx):
            num_items = model.rowCount(idx)
            for i in range(num_items):
                res = res or self.filterAcceptsRow(i, idx)
        return res


class TreeView(QTreeView):
    """
    QTreeView with some additional settings.
    """

    def __init__(self, *args, **kwargs) -> None:
        super().__init__(*args, **kwargs)
        self.setHeaderHidden(True)

        # drag and drop enabled
        self.setDragEnabled(True)
        self.setAcceptDrops(True)
        self.setDropIndicatorShown(True)
        self.setSelectionMode(QAbstractItemView.ExtendedSelection)
        self.setUniformRowHeights(True)

        # the view needs to be aware of the other view for the expanding purposes
        self.otherView = None

    def mousePressEvent(self, e: QtGui.QMouseEvent) -> None:  # noqa: N802
        super().mousePressEvent(e)
        self.handle_expand(e)

    def mouseDoubleClickEvent(self, e: QtGui.QMouseEvent) -> None:  # noqa: N802
        super().mouseDoubleClickEvent(e)
        self.handle_expand(e)

    def handle_expand(self, e: QtGui.QMouseEvent) -> None:
        """
        Check if click happened on the part with cause expanding the item in the view. If so also expand/collapse
        the same item in the other view.
        """
        clicked_index = self.indexAt(e.pos())
        if clicked_index.isValid():
            # if click is left from text then it is expand event
            vrect = self.visualRect(clicked_index)
            item_identation = vrect.x() - self.visualRect(self.rootIndex()).x()
            if e.pos().x() < item_identation:
                self.expand_in_other_view(clicked_index, not self.isExpanded(clicked_index))

    def expand_in_other_view(self, clicked_index: QModelIndex, is_expanded: bool) -> None:
        """
        Expand/collapse item which is same to the clicked_index in the other view. The other view is one of the left
        if click happened on the right and vice versa.

        Parameters
        ----------
        clicked_index
            Index of the item that was expanded/collapsed
        is_expanded
            This flag is true if the item was expanded otherwise if false (collapsed).
        """
        source_model_this = self.model().sourceModel()
        source_model_other = self.otherView.model().sourceModel()
        clicked_model_index = self.model().mapToSource(clicked_index)

        # register expanding change
        source_model_this.set_expanded(clicked_model_index, is_expanded)

        clicked_name = source_model_this.node_from_index(clicked_model_index).name

        index_other = source_model_other.index_from_name(clicked_name, source_model_other.rootItem)
        if index_other is not None:
            # if same item (item with the same name) exists in the other view
            self.otherView.setExpanded(index_other, is_expanded)

    def setModel(self, model: QtCore.QAbstractItemModel) -> None:  # noqa: N802
        """
        Override setModel that we can connect expandItem signal. This signal is emitted when model
        expand the item.
        """
        super().setModel(model)
        self.model().sourceModel().expand_item.connect(self.setExpanded)

    def setExpanded(self, index: QModelIndex, expand: bool) -> None:  # noqa: N802
        """
        Set item expanded/collapsed.

        Parameters
        ----------
        index
            Item's index
        expand
            This flag is true if the item was expanded otherwise if false (collapsed).
        """
        self.model().sourceModel().set_expanded(index, expand)
        index = self.model().mapFromSource(index)
        if expand:
            self.expand(index)
        else:
            self.collapse(index)


class OWMarkerGenes(widget.OWWidget):
    name = "Marker Genes"
    icon = 'icons/OWMarkerGenes.svg'
    priority = 130

    replaces = ['orangecontrib.single_cell.widgets.owmarkergenes.OWMarkerGenes']

    class Warning(widget.OWWidget.Warning):
        using_local_files = widget.Msg("Can't connect to serverfiles. Using cached files.")

    class Outputs:
        genes = widget.Output("Genes", Table)

    want_main_area = True
    want_control_area = True

    auto_commit = Setting(True)
    selected_source = Setting("")
    selected_organism = Setting("")
    selected_root_attribute = Setting(0)

    settingsHandler = MarkerGroupContextHandler()  # noqa: N815
    selected_genes = settings.ContextSetting([])

    settings_version = 2

    _data = None
    _available_sources = None

    def __init__(self) -> None:
        super().__init__()

        # define the layout
        main_area = QWidget(self.mainArea)
        self.mainArea.layout().addWidget(main_area)
        layout = QGridLayout()
        main_area.setLayout(layout)
        layout.setContentsMargins(4, 4, 4, 4)

        # filter line edit
        self.filter_line_edit = QLineEdit()
        self.filter_line_edit.setPlaceholderText("Filter marker genes")
        layout.addWidget(self.filter_line_edit, 0, 0, 1, 3)

        # define available markers view
        self.available_markers_view = TreeView()
        box = gui.vBox(self.mainArea, "Available markers", addToLayout=False)
        box.layout().addWidget(self.available_markers_view)
        layout.addWidget(box, 1, 0, 2, 1)

        # create selected markers view
        self.selected_markers_view = TreeView()
        box = gui.vBox(self.mainArea, "Selected markers", addToLayout=False)
        box.layout().addWidget(self.selected_markers_view)
        layout.addWidget(box, 1, 2, 2, 1)

        self.available_markers_view.otherView = self.selected_markers_view
        self.selected_markers_view.otherView = self.available_markers_view

        # buttons
        box = gui.vBox(self.mainArea, addToLayout=False, margin=0)
        layout.addWidget(box, 1, 1, 1, 1)
        self.move_button = gui.button(box, self, ">", callback=self._move_selected)

        self._init_description_area(layout)
        self._init_control_area()
        self._load_data()

    def _init_description_area(self, layout: QLayout) -> None:
        """
        Function define an info area with description of the genes and add it to the layout.
        """
        box = gui.widgetBox(self.mainArea, "Description", addToLayout=False)
        self.descriptionlabel = QTextBrowser(
            openExternalLinks=True, textInteractionFlags=(Qt.TextSelectableByMouse | Qt.LinksAccessibleByMouse)
        )
        box.setMaximumHeight(self.descriptionlabel.fontMetrics().height() * (NUM_LINES_TEXT + 3))

        # description filed
        self.descriptionlabel.setText("Select a gene to see information.")
        self.descriptionlabel.setFrameStyle(QTextBrowser.NoFrame)
        # no (white) text background
        self.descriptionlabel.viewport().setAutoFillBackground(False)
        box.layout().addWidget(self.descriptionlabel)
        layout.addWidget(box, 3, 0, 1, 3)

    def _init_control_area(self) -> None:
        """
        Function defines dropdowns and the button in the control area.
        """
        box = gui.widgetBox(self.controlArea, 'Database', margin=0)
        self.source_index = -1
        self.db_source_cb = gui.comboBox(box, self, 'source_index')
        self.db_source_cb.activated[int].connect(self._set_db_source_index)

        box = gui.widgetBox(self.controlArea, 'Organism', margin=0)
        self.organism_index = -1
        self.group_cb = gui.comboBox(box, self, 'organism_index')
        self.group_cb.activated[int].connect(self._set_group_index)

        box = gui.widgetBox(self.controlArea, 'Group by', margin=0)
        self.group_by_cb = gui.comboBox(
            box, self, 'selected_root_attribute', items=GROUP_BY_ITEMS, callback=self._setup
        )

        gui.rubber(self.controlArea)
        gui.auto_commit(self.controlArea, self, "auto_commit", "Commit")

    def sizeHint(self):
        return super().sizeHint().expandedTo(QSize(900, 500))

    @property
    def available_sources(self) -> dict:
        return self._available_sources

    @available_sources.setter
    def available_sources(self, value: dict) -> None:
        """
        Set _available_sources variable, add them to dropdown, and select the source that was previously.
        """
        self._available_sources = value

        items = sorted(list(value.keys()), reverse=True)  # panglao first
        try:
            idx = items.index(self.selected_source)
        except ValueError:
            idx = -1

        self.db_source_cb.clear()
        self.db_source_cb.addItems(items)

        if idx != -1:
            self.source_index = idx
            self.selected_source = items[idx]
        elif items:
            self.source_index = min(max(self.source_index, 0), len(items) - 1)

        self._set_db_source_index(self.source_index)

    @property
    def data(self) -> Table:
        return self._data

    @data.setter
    def data(self, value: Table):
        """
        Set the source data. The data is then filtered on the first meta column (group).
        Select set dropdown with the groups and select the one that was selected previously.
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
            group_values = sorted(group_values)  # human first

            try:
                idx = group_values.index(self.selected_organism)
            except ValueError:
                idx = -1

            self.group_cb.clear()
            self.group_cb.addItems(group_values)

            if idx != -1:
                self.organism_index = idx
                self.selected_organism = group_values[idx]
            elif group_values:
                self.organism_index = min(max(self.organism_index, 0), len(group_values) - 1)

            self._set_group_index(self.organism_index)

    def _load_data(self) -> None:
        """
        Collect available data sources (marker genes data sets).
        """
        self.Warning.using_local_files.clear()

        found_sources = {}
        try:
            found_sources.update(serverfiles.ServerFiles().allinfo(SERVER_FILES_DOMAIN))
        except requests.exceptions.ConnectionError:
            found_sources.update(serverfiles.allinfo(SERVER_FILES_DOMAIN))
            self.Warning.using_local_files()

        self.available_sources = {item.get('title').split(': ')[-1]: item for item in found_sources.values()}

    def _source_changed(self) -> None:
        """
        Respond on change of the source and download the data.
        """
        if self.available_sources:
            file_name = self.available_sources[self.selected_source]['filename']

            try:
                serverfiles.update(SERVER_FILES_DOMAIN, file_name)
            except requests.exceptions.ConnectionError:
                # try to update file. Ignore network errors.
                pass

            try:
                file_path = serverfiles.localpath_download(SERVER_FILES_DOMAIN, file_name)
            except requests.exceptions.ConnectionError as err:
                # Unexpected error.
                raise err
            self.data = Table.from_file(file_path)

    def _setup(self) -> None:
        """
        Setup the views with data.
        """
        self.closeContext()
        self.selected_genes = []
        self.openContext((self.selected_organism, self.selected_source))
        data_not_selected, data_selected = self._filter_data_group(self.data)

        # add model to available markers view
        group_by = GROUP_BY_ITEMS[self.selected_root_attribute]
        tree_model = TreeModel(data_not_selected, group_by)
        proxy_model = FilterProxyModel(self.filter_line_edit)
        proxy_model.setSourceModel(tree_model)

        self.available_markers_view.setModel(proxy_model)
        self.available_markers_view.selectionModel().selectionChanged.connect(
            partial(self._on_selection_changed, self.available_markers_view)
        )

        tree_model = TreeModel(data_selected, group_by)
        proxy_model = FilterProxyModel(self.filter_line_edit)
        proxy_model.setSourceModel(tree_model)
        self.selected_markers_view.setModel(proxy_model)

        self.selected_markers_view.selectionModel().selectionChanged.connect(
            partial(self._on_selection_changed, self.selected_markers_view)
        )
        self.selected_markers_view.model().sourceModel().data_added.connect(self._selected_markers_changed)
        self.selected_markers_view.model().sourceModel().data_removed.connect(self._selected_markers_changed)

        # update output and messages
        self._selected_markers_changed()

    def _filter_data_group(self, data: Table) -> Tuple[Table, Tuple]:
        """
        Function filter the table based on the selected group (Mouse, Human) and divide them in two groups based on
        selected_data variable.

        Parameters
        ----------
        data
            Table to be filtered

        Returns
        -------
        data_not_selected
            Data that will initially be in available markers view.
        data_selected
            Data that will initially be in selected markers view.
        """
        group = data.domain.metas[0]
        gvec = data.get_column_view(group)[0]

        if group.is_string:
            mask = gvec == self.selected_organism
        else:
            mask = gvec == self.organism_index
        data = data[mask]

        # divide data based on selected_genes variable (context)
        unique_gene_names = np.core.defchararray.add(
            data.get_column_view("Entrez ID")[0].astype(str), data.get_column_view("Cell Type")[0].astype(str)
        )
        mask = np.isin(unique_gene_names, self.selected_genes)
        data_not_selected = data[~mask]
        data_selected = data[mask]
        return data_not_selected, data_selected

    def commit(self) -> None:
        rows = self.selected_markers_view.model().sourceModel().rootItem.get_data_rows()
        if len(rows) > 0:
            metas = [r.metas for r in rows]
            data = Table.from_numpy(self.data.domain, np.empty((len(metas), 0)), metas=np.array(metas))
            # always false for marker genes data tables in single cell
            data.attributes[GENE_AS_ATTRIBUTE_NAME] = False
            # set taxonomy id in data.attributes
            data.attributes[TAX_ID] = MAP_GROUP_TO_TAX_ID.get(self.selected_organism, '')
            # set column id flag
            data.attributes[GENE_ID_COLUMN] = "Entrez ID"
            data.name = 'Marker Genes'
        else:
            data = None
        self.Outputs.genes.send(data)

    def _update_description(self, view: TreeView) -> None:
        """
        Upate the description about the gene. Only in case when one gene is selected.
        """
        selection = self._selected_rows(view)
        qmodel = view.model().sourceModel()

        if len(selection) > 1 or len(selection) == 0 or qmodel.node_from_index(selection[0]).data_row is None:
            self.descriptionlabel.setText("Select a gene to see information.")
        else:
            data_row = qmodel.node_from_index(selection[0]).data_row
            self.descriptionlabel.setHtml(
                f"<b>Gene name:</b> {data_row['Name']}<br/>"
                f"<b>Entrez ID:</b> {data_row['Entrez ID']}<br/>"
                f"<b>Cell Type:</b> {data_row['Cell Type']}<br/>"
                f"<b>Function:</b> {data_row['Function']}<br/>"
                f"<b>Reference:</b> <a href='{data_row['URL']}'>{data_row['Reference']}</a>"
            )

    def _update_data_info(self) -> None:
        """
        Updates output info in the control area.
        """
        sel_model = self.selected_markers_view.model().sourceModel()
        self.info.set_output_summary(f"Selected: {str(len(sel_model))}")

    # callback functions

    def _selected_markers_changed(self) -> None:
        """
        This function is called when markers in the selected view are added or removed.
        """
        rows = self.selected_markers_view.model().sourceModel().rootItem.get_data_rows()
        self.selected_genes = [row["Entrez ID"].value + row["Cell Type"].value for row in rows]
        self._update_data_info()
        self.commit()

    def _on_selection_changed(self, view: TreeView) -> None:
        """
        When selection in one of the view changes in a view button should change a sign in the correct direction and
        other view should reset the selection. Also gene description is updated.
        """
        self.move_button.setText(">" if view is self.available_markers_view else "<")
        if view is self.available_markers_view:
            self.selected_markers_view.clearSelection()
        else:
            self.available_markers_view.clearSelection()
        self._update_description(view)

    def _set_db_source_index(self, source_index: int) -> None:
        """
        Set the index of selected database source - index in a dropdown.
        """
        self.source_index = source_index
        self.selected_source = self.db_source_cb.itemText(source_index)
        self._source_changed()

    def _set_group_index(self, group_index: int) -> None:
        """
        Set the index of organism - index in a dropdown.
        """
        self.organism_index = group_index
        self.selected_organism = self.group_cb.itemText(group_index)
        self._setup()

    def _move_selected(self) -> None:
        """
        Move selected genes when button clicked.
        """
        if self._selected_rows(self.selected_markers_view):
            self._move_selected_from_to(self.selected_markers_view, self.available_markers_view)
        elif self._selected_rows(self.available_markers_view):
            self._move_selected_from_to(self.available_markers_view, self.selected_markers_view)

    # support functions for callbacks

    def _move_selected_from_to(self, src: TreeView, dst: TreeView) -> None:
        """
        Function moves items from src model to dst model.
        """
        selected_items = self._selected_rows(src)

        src_model = src.model().sourceModel()
        dst_model = dst.model().sourceModel()

        # move data as mimeData from source to destination tree view
        mime_data = src_model.mimeData(selected_items)
        # remove nodes from the source view
        src_model.remove_node_list(selected_items)
        dst_model.dropMimeData(mime_data, Qt.MoveAction, -1, -1)

    @staticmethod
    def _selected_rows(view: TreeView) -> List[QModelIndex]:
        """
        Return the selected rows in the view.
        """
        rows = view.selectionModel().selectedRows()
        return list(map(view.model().mapToSource, rows))

    @classmethod
    def migrate_settings(cls, settings, version=0):
        def migrate_to_version_2():
            settings["selected_source"] = settings.pop("selected_db_source", "")
            settings["selected_organism"] = settings.pop("selected_group", "")
            if "context_settings" in settings:
                for co in settings["context_settings"]:
                    co.values["selected_genes"] = [g[0] + g[1] for g in co.values["selected_genes"]]

        if version < 2:
            migrate_to_version_2()


if __name__ == "__main__":
    from orangewidget.utils.widgetpreview import WidgetPreview

    WidgetPreview(OWMarkerGenes).run()
