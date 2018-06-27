""" GUI utils for widgets """
from collections import namedtuple
from numbers import Integral, Real
from typing import Sequence

from AnyQt.QtWidgets import (
    QFrame, QStyledItemDelegate
)

from AnyQt.QtCore import (
    Qt, QSortFilterProxyModel,
)

from .list_completer import TokenListCompleter
from .label_selection import (
    standarditem_from, group_candidates, itemselection, RowGroup, ColumnGroup, group_selection_mask,
    LabelSelectionWidget
)


# Creates line separator
def horizontal_line():
    line = QFrame()
    line.setFrameShape(QFrame.HLine)
    line.setFrameShadow(QFrame.Sunken)
    return line


class FilterProxyModel(QSortFilterProxyModel):
    """
    A simple filter proxy model with settable filter predicates
    Example
    -------
    >>> proxy = FilterProxyModel()
    >>> proxy.set_filters([
    ...     FilterProxyModel.Filter(0, Qt.DisplayRole, lambda value: value < 1)
    ... ])
    """
    Filter = namedtuple("Filter", ["column", "role", "predicate"])

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.__filters = []

    def reset_filters(self):
        self.__filters = []
        self.invalidateFilter()

    def set_filters(self, filters):
        # type: (Sequence[FilterProxyModel.Filter]) -> None

        filters = [FilterProxyModel.Filter(f.column, f.role, f.predicate) for f in filters]
        self.__filters = filters
        self.invalidateFilter()

    def filterAcceptsRow(self, row, parent):
        source = self.sourceModel()

        def apply(f):
            index = source.index(row, f.column, parent)
            data = source.data(index, f.role)
            try:
                return f.predicate(data)
            except (TypeError, ValueError):
                return False

        return all(apply(f) for f in self.__filters)


class NumericalColumnDelegate(QStyledItemDelegate):
    """
    An Item delegate for displaying numerical columns
    """
    def __init__(self, parent=None, precision=4, notation='f'):
        super().__init__(parent)
        self.precision = precision
        self.notation = notation

    def displayText(self, value, locale):
        if isinstance(value, Integral):
            return locale.toString(int(value))
        elif isinstance(value, Real):
            return locale.toString(float(value), self.notation, self.precision)
        else:
            return super().displayText(value, locale)

    def initStyleOption(self, option, index):
        super().initStyleOption(option, index)
        align = index.data(Qt.TextAlignmentRole)
        data = index.data(Qt.DisplayRole)
        if align is None and isinstance(data, Real):
            # Right align if the model does not specify otherwise
            option.displayAlignment = Qt.AlignRight | Qt.AlignVCenter
