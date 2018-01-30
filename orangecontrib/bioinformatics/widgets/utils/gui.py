import warnings
import re

from itertools import chain
from collections import defaultdict
from functools import singledispatch

from AnyQt import QtGui, QtCore, QtWidgets
from AnyQt.QtWidgets import QCompleter
from AnyQt.QtCore import Qt, QObject, QStringListModel, pyqtSignal as Signal


@singledispatch
def standarditem_from(obj):
    item = QtGui.QStandardItem(str(obj))
    item.setData(obj, Qt.UserRole)
    item.setFlags(item.flags() & ~Qt.ItemIsEditable)
    return item


class LabelSelectionWidget(QtWidgets.QWidget):
    """
    A widget for selection of label values.

    The widget displays the contents of the model with two widgets:

    * The top level model items are displayed in a combo box.
    * The children of the current selected top level item are
      displayed in a subordinate list view.

    .. note:: This is not a QAbstractItemView subclass.
    """
    #: Current group/root index has changed.
    groupChanged = Signal(int)
    #: Selection for the current group/root has changed.
    groupSelectionChanged = Signal()

    def __init__(self, parent=None, **kwargs):
        super().__init__(parent, **kwargs)
        self.__model = None
        self.__selectionMode = QtWidgets.QListView.ExtendedSelection

        self.__currentIndex = -1
        self.__selections = {}

        layout = QtWidgets.QVBoxLayout()
        layout.setContentsMargins(0, 0, 0, 0)

        def group_box(title):
            box = QtWidgets.QGroupBox(title)
            box.setFlat(True)
            lay = QtWidgets.QVBoxLayout()
            lay.setContentsMargins(0, 0, 0, 0)
            box.setLayout(lay)
            return box

        self.labels_combo = QtWidgets.QComboBox()
        self.values_view = QtWidgets.QListView(
            selectionMode=self.__selectionMode
        )

        self.labels_combo.currentIndexChanged.connect(
            self.__onCurrentIndexChanged)

        l_box = group_box(self.tr("Label"))
        v_box = group_box(self.tr("Values"))

        l_box.layout().addWidget(self.labels_combo)
        v_box.layout().addWidget(self.values_view)

        layout.addWidget(l_box)
        layout.addWidget(v_box)

        self.setLayout(layout)

        self.setSizePolicy(QtWidgets.QSizePolicy.Expanding,
                           QtWidgets.QSizePolicy.Expanding)

    def clear(self):
        """ Clear the widget/model (same as ``setModel(None)``).
        """
        if self.__model is not None:
            self.values_view.selectionModel().clearSelection()
            self.values_view.selectionModel().selectionChanged.disconnect(
                self.__onSelectionChanged)

            self.values_view.setModel(None)
            self.labels_combo.setModel(
                QtGui.QStandardItemModel(self.labels_combo))
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
        self.values_view.selectionModel().selectionChanged.connect(
            self.__onSelectionChanged)

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
                    if parent == self.__model.index(parent.row(),
                                                    parent.column()):
                        pind[parent.row()].append(
                            QtCore.QPersistentModelIndex(index))
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
        selection = QtCore.QItemSelection()
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
            return QtCore.QItemSelection()

    def __onCurrentIndexChanged(self, index):
        self.__storeSelection(self.__currentIndex,
                              self.values_view.selectedIndexes())

        self.__currentIndex = index
        if self.__model is not None:
            root = self.__model.index(index, 0)
            self.values_view.setRootIndex(root)

            self.__restoreSelection()
        self.groupChanged.emit(index)

    def __onSelectionChanged(self, old, new):
        self.__storeSelection(self.__currentIndex,
                              self.values_view.selectedIndexes())

        self.groupSelectionChanged.emit()

    def __storeSelection(self, groupind, indices):
        # Store current values selection for the current group
        groupind = self.__currentIndex
        indices = [QtCore.QPersistentModelIndex(ind)
                   for ind in self.values_view.selectedIndexes()]
        self.__selections[groupind] = indices

    def __restoreSelection(self):
        # Restore previous selection for root (if available)
        assert self.__model is not None
        groupind = self.__currentIndex
        root = self.__model.index(groupind, 0)
        sel = self.__selections.get(groupind, [])
        indices = [self.__model.index(pind.row(), pind.column(), root)
                   for pind in sel if pind.isValid() and pind.parent() == root]

        selection = QtCore.QItemSelection()
        for ind in indices:
            selection.select(ind, ind)
        self.values_view.selectionModel().select(
            selection, QtCore.QItemSelectionModel.ClearAndSelect)

    def sizeHint(self):
        """Reimplemented from QWidget.sizeHint"""
        return QtCore.QSize(100, 200)


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
    selection = QtCore.QItemSelection()
    for index in modelindexlist:
        selection.select(index, index)
    return selection


class TokenListCompleter(QCompleter):
    """
    A completer allowing completion of multiple words (tokens)

    Example
    -------
    completer = TokenListCompleter()
    completer.setTokenList(["foo", "bar", "baz"])
    completer.setCompletionPrefix("foo b")
    for i in range(completer.completionCount()):
         completer.setCurrentRow(i)
         print(completer.currentCompletion())

    """
    def __init__(self, *args, **kwargs):
        separator = kwargs.pop("separator", " ")
        super().__init__(*args, **kwargs)
        self.__tokenList = []
        self.__completerModel = None
        self.__separator = separator
        # The current 'known' completion prefix (tracked in splitPath)
        self.__currentKnownPrefix = ""
        self.setModelSorting(QCompleter.CaseSensitivelySortedModel)
        super().setModel(QStringListModel(self))

    def setTokenList(self, tokenList):
        """
        Set the token word completion list

        Equivalent to
        `completer.setModel(QStringListModel(tokenList, completer))`

        Note that this is the preferred method for setting the completion
        model.

        Parameters
        ----------
        tokenList : List[str]
        """
        tokenList = list(sorted(set(tokenList)))
        self.setModel(QStringListModel(tokenList, self))

    def tokenList(self):
        """
        Return the current token word completion list.

        Return
        ------
        tokens : List[str]
        """
        return list(self.__tokenList)

    def setModel(self, model):
        """
        Reimplemented.

        Parameters
        ----------
        model : QAbstractItemModel
        """
        if model is self.__completerModel:
            return

        if self.__completerModel is not None:
            self.__completerModel.dataChanged.disconnect(
                self.__initDynamicModel)
            self.__completerModel.rowsInserted.disconnect(
                self.__initDynamicModel)
            self.__completerModel.rowsRemoved.disconnect(
                self.__initDynamicModel)

            if QObject.parent(self.__completerModel) is self:
                self.__completerModel.deleteLater()
            self.__completerModel = None

        self.__completerModel = model

        if self.__completerModel is not None:
            self.__completerModel.dataChanged.connect(
                self.__initDynamicModel)
            self.__completerModel.rowsInserted.connect(
                self.__initDynamicModel)
            self.__completerModel.rowsRemoved.connect(
                self.__initDynamicModel)

        self.__initDynamicModel()

    def model(self):
        return self.__completerModel

    def setSeparator(self, separator):
        if self.__separator != separator:
            self.__separator = separator

    def separator(self):
        return self.__separator

    # QCompleter::setCompleterPrefix is not a virtual method, and is called
    # directly by the QLineEdit to set the prefix. There does not appear to
    # be any other direct way that a QCompleter subclass could be notified
    # of the completion prefix change.
    # However in QCompleter::setCompleterPrefix
    # the prefix is immediately split by splitPath which is virtual.
    # So we use this to track changes in the current prefix.
    def splitPath(self, prefix):
        """reimplemented."""
        items = super().splitPath(prefix)
        if self.__currentKnownPrefix != self.completionPrefix():
            self.__currentKnownPrefix = self.completionPrefix()
            self.__updateEffectiveCompleterModel()
        return items

    def __updateEffectiveCompleterModel(self):
        """
        Prefix all items in the effective completer model with the current
        prefix to enable the completion of multiple keywords.
        """
        prefix = self.completionPrefix()
        model = super().model()
        assert isinstance(model, QStringListModel)
        separator = self.__separator
        if not prefix.rstrip().endswith(separator) and separator in prefix:
            # We are in the middle of a second (or later) word completion
            prefix, rest = prefix.rsplit(separator, 1)
            # preserve any whitespace after the seperator
            match = re.match(r"^(?P<whitespace>\s*)(?P<rest>.*)", rest)
            if match:
                whitespace, rest = match.groups()
            else:
                whitespace = ""
            items = [prefix + separator + whitespace + item
                     for item in self.__tokenList]
        else:
            # Either completing the the first word or immediately after
            # the separator. In both cases the original tokens list is
            # the active completion list
            items = self.__tokenList
        itemset = set(items)
        current = model.stringList()

        if itemset != set(current):
            model.setStringList(list(sorted(itemset)))

    def __initDynamicModel(self):
        """
        [Re]Initialize the effective completion model from the client
        supplied model
        """
        model = self.__completerModel
        if model is not None:
            if isinstance(model, QStringListModel):
                tokens = model.stringList()
            else:
                tokens = [model.data(model.index(row, 0), Qt.DisplayRole)
                          for row in range(model.rowCount())]
                tokens = [str(token) for token in filter(None, tokens)]
        else:
            tokens = []

        tokenList = list(sorted(set(tokens)))
        if tokenList != self.__tokenList:
            self.__tokenList = tokenList
            self.__updateEffectiveCompleterModel()
