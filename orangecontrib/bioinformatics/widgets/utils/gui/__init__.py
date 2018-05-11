""" GUI utils for widgets """
from AnyQt.QtWidgets import (
    QFrame
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
