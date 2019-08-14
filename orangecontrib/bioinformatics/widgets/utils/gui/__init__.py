""" GUI utils for widgets """
from typing import Sequence
from numbers import Real, Integral
from collections import namedtuple

from AnyQt.QtGui import QTextDocument, QAbstractTextDocumentLayout
from AnyQt.QtCore import Qt, QSize, QSortFilterProxyModel
from AnyQt.QtWidgets import QFrame, QStyle, QApplication, QStyledItemDelegate, QStyleOptionViewItem

from .gene_sets import GeneSetsSelection
from .gene_scoring import GeneScoringWidget, gene_scoring_method
from .list_completer import TokenListCompleter
from .label_selection import (
    RowGroup,
    ColumnGroup,
    LabelSelectionWidget,
    itemselection,
    group_candidates,
    standarditem_from,
    group_selection_mask,
)

__all__ = (
    'GeneSetsSelection',
    'GeneScoringWidget',
    'gene_scoring_method',
    'TokenListCompleter',
    'RowGroup',
    'ColumnGroup',
    'LabelSelectionWidget',
    'itemselection',
    'group_candidates',
    'standarditem_from',
    'group_selection_mask',
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


class HTMLDelegate(QStyledItemDelegate):
    """
    https://stackoverflow.com/questions/1956542/how-to-make-item-view-render-rich-html-text-in-qt
    https://stackoverflow.com/questions/2375763/how-to-open-an-url-in-a-qtableview
    """

    def sizeHint(self, option, index):
        options = QStyleOptionViewItem(option)
        gene_obj = index.data(Qt.DisplayRole)
        self.initStyleOption(options, index)

        doc = QTextDocument()
        doc.setHtml(gene_obj.to_html())
        doc.setTextWidth(options.rect.width() - 10)

        return QSize(doc.idealWidth(), doc.size().height())

    def paint(self, painter, option, index):
        options = QStyleOptionViewItem(option)
        row_obj = index.data(Qt.DisplayRole)
        self.initStyleOption(options, index)
        # print(option.rect.width())
        style = QApplication.style() if options.widget is None else options.widget.style()

        doc = QTextDocument()
        doc.setHtml(row_obj.to_html())
        doc.setTextWidth(option.rect.width() - 10)

        # doc.setPageSize(300)
        # print(doc.loadResource(3))

        options.text = ""
        style.drawControl(QStyle.CE_ItemViewItem, options, painter)

        ctx = QAbstractTextDocumentLayout.PaintContext()

        text_rect = style.subElementRect(QStyle.SE_ItemViewItemText, options)
        painter.save()
        painter.translate(text_rect.topLeft())
        painter.setClipRect(text_rect.translated(-text_rect.topLeft()))
        doc.documentLayout().draw(painter, ctx)

        painter.restore()
