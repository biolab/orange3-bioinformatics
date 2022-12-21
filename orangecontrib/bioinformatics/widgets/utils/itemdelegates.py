from AnyQt.QtCore import QModelIndex
from AnyQt.QtWidgets import QStyleOptionViewItem

from Orange.widgets.gui import LinkRole, LinkStyledItemDelegate


class LinkStyledItemDelegate(LinkStyledItemDelegate):
    def initStyleOption(self, option: QStyleOptionViewItem, index: QModelIndex) -> None:
        super().initStyleOption(option, index)
        if index.data(LinkRole) is not None:
            option.font.setUnderline(True)

    def paint(self, painter, option, index):
        self.initStyleOption(option, index)
        super().paint(painter, option, index)
