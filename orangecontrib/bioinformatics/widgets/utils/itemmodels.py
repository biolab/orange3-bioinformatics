from typing import Any, Union, Mapping, Callable, Sequence, NamedTuple

from AnyQt.QtCore import Qt, QModelIndex

from orangewidget.utils.itemmodels import AbstractSortTableModel

ItemDataRole = Union[int, Qt.ItemDataRole]


class TableModel(AbstractSortTableModel):
    class Column(NamedTuple):
        #: The column's header data
        header: Mapping[ItemDataRole, Any]
        data: Mapping[ItemDataRole, Callable[[int], Any]]
        setData: Mapping[ItemDataRole, Callable[[int, Any], bool]] = {}

    __slots__ = ("__columns", "__rowCount", "__columnCount")

    def __init__(self, length: int, columns: Sequence[Column], *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.__columns = tuple(columns)
        self.__rowCount = length
        self.__columnCount = len(self.__columns)

    def setColumns(self, lenght: int, columns: Sequence[Column]):
        self.beginResetModel()
        self.__columns = tuple(columns)
        self.__columnCount = len(columns)
        self.__rowCount = lenght
        self.endResetModel()

    def columns(self) -> Sequence[Column]:
        return list(self.__columns)

    def rowCount(self, parent: QModelIndex = QModelIndex()) -> int:
        return 0 if parent.isValid() else self.__rowCount

    def columnCount(self, parent: QModelIndex = QModelIndex()) -> int:
        return 0 if parent.isValid() else self.__columnCount

    def data(self, index: QModelIndex, role: int = Qt.DisplayRole) -> Any:
        if not index.isValid():
            return None
        row, col = index.row(), index.column()
        N, M = self.__rowCount, self.__columnCount
        if not 0 <= row < N and 0 <= col < M:
            return None
        row = self.mapToSourceRows(row)
        column = self.__columns[col]
        dispatch = column.data.get(role, None)
        if dispatch is not None:
            return dispatch(row)
        else:
            return None

    def setData(self, index: QModelIndex, value: Any, role: int = Qt.EditRole) -> bool:
        if not index.isValid():
            return False
        row, col = index.row(), index.column()
        N, M = self.__rowCount, self.__columnCount
        if not 0 <= row < N and 0 <= col < M:
            return False
        row = self.mapToSourceRows(row)
        column = self.__columns[col]
        dispatch = column.setData.get(role, None)
        if dispatch is not None:
            return dispatch(row, value)
        else:
            return False

    def headerData(self, section: int, orientation: Qt.Orientation, role: int = Qt.DisplayRole) -> Any:
        if orientation == Qt.Horizontal:
            if not 0 <= section < self.__columnCount:
                return None
            column = self.__columns[section]
            return column.header.get(role, None)
        return super().headerData(section, orientation, role)
