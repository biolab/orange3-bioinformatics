import warnings
from functools import singledispatch
from itertools import chain, starmap
from collections import namedtuple, defaultdict
from xml.sax.saxutils import escape

import numpy as np

from AnyQt.QtGui import QStandardItem, QStandardItemModel
from AnyQt.QtCore import Qt, QSize, QItemSelection, QItemSelectionModel, QPersistentModelIndex
from AnyQt.QtCore import pyqtSignal as Signal
from AnyQt.QtWidgets import QWidget, QComboBox, QGroupBox, QListView, QSizePolicy, QVBoxLayout

from Orange.data import Variable, StringVariable, DiscreteVariable, ContinuousVariable
from Orange.widgets import gui
from Orange.widgets.utils import itemmodels


@singledispatch
def standarditem_from(obj):
    item = QStandardItem(str(obj))
    item.setData(obj, Qt.UserRole)
    item.setFlags(item.flags() & ~Qt.ItemIsEditable)
    return item


@standarditem_from.register(Variable)
def standarditem_from_var(var):
    item = QStandardItem(var.name)
    _, icon = gui.attributeItem(var)
    item.setIcon(icon)
    item.setToolTip(itemmodels.Variable)


def variable_tooltip(var):
    if isinstance(var, DiscreteVariable):
        return discrete_variable_tooltip(var)
    elif isinstance(var, ContinuousVariable):
        return continuous_variable_toltip(var)
    elif isinstance(var, StringVariable):
        return string_variable_tooltip(var)


def variable_labels_tooltip(var):
    text = ""
    if var.attributes:
        items = [(escape(key), escape(value)) for key, value in var.attributes.items()]
        labels = list(starmap("{!s} = {!s}".format, items))
        text += "<br/>Variable Labels:<br/>"
        text += "<br/>".join(labels)
    return text


def discrete_variable_tooltip(var):
    text = "<b>%s</b><br/>Discrete with %i values: " % (escape(var.name), len(var.values))
    text += ", ".join("%r" % escape(v) for v in var.values)
    text += variable_labels_tooltip(var)
    return text


def continuous_variable_toltip(var):
    text = "<b>%s</b><br/>Continuous" % escape(var.name)
    text += variable_labels_tooltip(var)
    return text


def string_variable_tooltip(var):
    text = "<b>%s</b><br/>String" % escape(var.name)
    text += variable_labels_tooltip(var)
    return text


def itemselection(modelindexlist):
    """
    Return an QtCore.QItemSelection from QModelIndex list

    Parameters
    ----------
    modelindexlist : list of QtCore.QModelIndex
        Selected model indices.

    Returns
    -------
    selection : QtCore.QItemSelection
    """
    selection = QItemSelection()
    for index in modelindexlist:
        selection.select(index, index)
    return selection


ColumnGroup = namedtuple("ColumnGroup", ["name", "key", "values"])
RowGroup = namedtuple("RowGroup", ["name", "var", "values"])


def group_candidates(data):
    items = [attr.attributes.items() for attr in data.domain.attributes]
    items = list(chain(*items))

    targets = defaultdict(set)
    for label, value in items:
        targets[label].add(value)

    # Need at least 2 distinct values or key
    targets = [(key, sorted(vals)) for key, vals in targets.items() if len(vals) >= 2]

    column_groups = [ColumnGroup(key, key, values) for key, values in sorted(targets)]

    disc_vars = [
        var
        for var in data.domain.class_vars + data.domain.metas
        if isinstance(var, DiscreteVariable) and len(var.values) >= 2
    ]

    row_groups = [RowGroup(var.name, var, var.values) for var in disc_vars]
    return column_groups, row_groups


@standarditem_from.register(ColumnGroup)
def standarditem_from_columngroup(colgroup):
    item = QStandardItem(colgroup.name)
    item.setToolTip("Split by column label: '{!s}'".format(escape(colgroup.name)))
    item.setFlags(item.flags() & ~Qt.ItemIsEditable)
    item.setData(colgroup, Qt.UserRole)
    children = [standarditem_from(val) for val in colgroup.values]
    item.appendRows(children)
    return item


@standarditem_from.register(RowGroup)
def standarditem_from_rowgroup(rowgroup):
    item = QStandardItem(rowgroup.name)
    icon, _ = gui.attributeItem(rowgroup.var)
    item.setIcon(icon)
    item.setToolTip(variable_tooltip(rowgroup.var))
    item.setData(rowgroup, Qt.UserRole)
    item.setFlags(item.flags() & ~Qt.ItemIsEditable)
    children = [standarditem_from(val) for val in rowgroup.values]
    item.appendRows(children)
    return item


def group_selection_mask(data, group, indices):
    """ Return the selection masks for the group.
    """
    if isinstance(group, ColumnGroup):
        selected = [group.values[i] for i in indices]
        target = {(group.key, value) for value in selected}
        _i = [bool(set(var.attributes.items()).intersection(target)) for var in data.domain.attributes]
        return np.array(_i, dtype=bool)
    elif isinstance(group, RowGroup):
        target = set(indices)
        x, _ = data.get_column_view(group.var)
        _i = np.zeros_like(x, dtype=bool)
        for i in target:
            _i |= x == i
        return _i
    else:
        raise TypeError("ColumnGroup or RowGroup expected, got {}".format(type(group).__name__))


class LabelSelectionWidget(QWidget):
    """
    A widget for selection of label values.

    The widget displays the contents of the model with two widgets:

    * The top level model items are displayed in a combo box.
    * The children of the current selected top level item are
      displayed in a subordinate list view.

    .. note:: This is not a QAbstractItemView subclass.
    """

    # Current group/root index has changed.
    groupChanged = Signal(int)
    # Selection for the current group/root has changed.
    groupSelectionChanged = Signal(int)

    def __init__(self):
        super().__init__()
        self.__model = None
        self.__selectionMode = QListView.ExtendedSelection

        self.__currentIndex = -1
        self.__selections = {}

        self.__parent = None

        self.targets = []

        layout = QVBoxLayout()
        layout.setContentsMargins(0, 0, 0, 0)

        def group_box(title):
            box = QGroupBox(title)
            box.setFlat(True)
            lay = QVBoxLayout()
            lay.setContentsMargins(0, 0, 0, 0)
            box.setLayout(lay)
            return box

        self.labels_combo = QComboBox()
        self.values_view = QListView(selectionMode=self.__selectionMode)

        self.labels_combo.currentIndexChanged.connect(self.__onCurrentIndexChanged)

        l_box = group_box(self.tr("Label"))
        v_box = group_box(self.tr("Values"))

        l_box.layout().addWidget(self.labels_combo)
        v_box.layout().addWidget(self.values_view)

        layout.addWidget(l_box)
        layout.addWidget(v_box)

        self.setLayout(layout)

        self.setSizePolicy(QSizePolicy.Expanding, QSizePolicy.Expanding)

    def set_selection(self):
        model = self.model()

        # FIX: This assumes fixed target key/values order.
        indices = [
            model.index(i, 0, model.index(keyind, 0))
            for keyind, selected in enumerate(self.__parent.stored_selections)
            for i in selected
        ]

        selection = itemselection(indices)
        self.setCurrentGroupIndex(self.__parent.current_group_index)
        self.setSelection(selection)

    def selected_split(self):
        group_index = self.currentGroupIndex()
        if not (0 <= group_index < len(self.targets)):
            return None, []

        group = self.targets[group_index]
        selection = [ind.row() for ind in self.currentGroupSelection().indexes()]
        return group, selection

    def set_data(self, parent, data):
        """ Initialize widget state after receiving new data.
        """
        self.__parent = parent
        # self.__currentIndex = parent.current_group_index
        if data is not None:
            column_groups, row_groups = group_candidates(data)

            self.targets = column_groups + row_groups

            model = QStandardItemModel()
            for item in [standarditem_from(desc) for desc in self.targets]:
                model.appendRow(item)

            self.setModel(model)
        else:
            self.targets = []
            self.setModel(None)

    def clear(self):
        """ Clear the widget/model (same as ``setModel(None)``).
        """
        if self.__model is not None:
            self.values_view.selectionModel().clearSelection()
            self.values_view.selectionModel().selectionChanged.disconnect(self.__onSelectionChanged)

            self.values_view.setModel(None)
            self.labels_combo.setModel(QStandardItemModel(self.labels_combo))
            self.__currentIndex = -1
            self.__selections = {}
            self.__model = None

    def setModel(self, model):
        """
        Set the source model for display.

        The model should be a tree model with depth 2.

        Parameters
        ----------
        model : QtCore.QAbstractItemModel
        """
        if model is self.__model:
            return

        self.clear()

        if model is None:
            return

        self.__model = model
        self.values_view.setModel(model)
        self.values_view.setRootIndex(model.index(0, 0))

        self.values_view.selectionModel().selectionChanged.connect(self.__onSelectionChanged)
        # will emit the currentIndexChanged (if the model is not empty)
        self.labels_combo.setModel(model)

    def model(self):
        """
        Return the current model.

        Returns
        -------
        model : QtCore.QAbstractItemModel
        """
        return self.__model

    def setCurrentGroupIndex(self, index):
        """
        Set the current selected group/root index row.

        Parameters
        ----------
        index : int
            Group index.
        """
        self.labels_combo.setCurrentIndex(index)

    def currentGroupIndex(self):
        """
        Return the current selected group/root index row.

        Returns
        -------
        row : index
            Current group row index (-1 if there is no current index)
        """
        return self.labels_combo.currentIndex()

    def setSelection(self, selection):
        """
        Set the model item selection.

        Parameters
        ----------
        selection : QtCore.QItemSelection
            Item selection.
        """
        if self.values_view.selectionModel() is not None:
            indices = selection.indexes()
            pind = defaultdict(list)

            for index in indices:
                parent = index.parent()
                if parent.isValid():
                    if parent == self.__model.index(parent.row(), parent.column()):
                        pind[parent.row()].append(QPersistentModelIndex(index))
                    else:
                        warnings.warn("Die Die Die")
                else:
                    # top level index
                    pass

            self.__selections = pind
            self.__restoreSelection()

    def selection(self):
        """
        Return the item selection.

        Returns
        -------
        selection : QtCore.QItemSelection
        """
        selection = QItemSelection()
        if self.__model is None:
            return selection

        for pind in chain(*self.__selections.values()):
            ind = self.__model.index(pind.row(), pind.column(), pind.parent())
            if ind.isValid():
                selection.select(ind, ind)
        return selection

    def currentGroupSelection(self):
        """
        Return the item selection for the current group only.
        """
        if self.values_view.selectionModel() is not None:
            return self.values_view.selectionModel().selection()
        else:
            return QItemSelection()

    def __onCurrentIndexChanged(self, index):
        self.__storeSelection(self.__currentIndex, self.values_view.selectedIndexes())

        self.__currentIndex = index
        if self.__model is not None:
            root = self.__model.index(index, 0)
            self.values_view.setRootIndex(root)

            self.__restoreSelection()
        self.groupChanged.emit(index)

    def __onSelectionChanged(self, old, new):
        self.__storeSelection(self.__currentIndex, self.values_view.selectedIndexes())

        self.groupSelectionChanged.emit(self.__currentIndex)

    def __storeSelection(self, groupind, indices):
        # Store current values selection for the current group
        groupind = self.__currentIndex
        indices = [QPersistentModelIndex(ind) for ind in self.values_view.selectedIndexes()]
        self.__selections[groupind] = indices

    def __restoreSelection(self):
        # Restore previous selection for root (if available)
        assert self.__model is not None
        groupind = self.__currentIndex
        root = self.__model.index(groupind, 0)
        sel = self.__selections.get(groupind, [])
        indices = [
            self.__model.index(pind.row(), pind.column(), root)
            for pind in sel
            if pind.isValid() and pind.parent() == root
        ]

        selection = QItemSelection()
        for ind in indices:
            selection.select(ind, ind)
        self.values_view.selectionModel().select(selection, QItemSelectionModel.ClearAndSelect)

    def sizeHint(self):
        """Reimplemented from QWidget.sizeHint"""
        return QSize(100, 200)
