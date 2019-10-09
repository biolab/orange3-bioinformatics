from typing import Callable
from collections import namedtuple, defaultdict

import numpy as np

from AnyQt.QtCore import QObject, QItemSelection, QItemSelectionModel
from AnyQt.QtCore import pyqtSignal as Signal
from AnyQt.QtWidgets import QListView, QListWidget

from Orange.data import Table, DiscreteVariable
from Orange.widgets.gui import comboBox, widgetBox, doubleSpin
from Orange.widgets.widget import OWComponent
from Orange.widgets.settings import ContextSetting

from orangecontrib.bioinformatics.utils.statistics import score_t_test, score_mann_whitney, score_hypergeometric_test
from orangecontrib.bioinformatics.widgets.utils.data import TableAnnotation

column_group = namedtuple('ColumnGroup', ['name', 'key', 'values'])
row_group = namedtuple('RowGroup', ['name', 'var', 'values'])


class GeneScoringComponent(OWComponent, QObject):

    # Current group/root index has changed.
    group_changed = Signal(int)
    # Selection for the current group/root has changed.
    selection_changed = Signal(int)
    # Scoring method has changed.
    score_method_changed = Signal(int)
    # Expression threshold changed
    expression_threshold_changed = Signal(float)
    # have one general signal if other are not needed.
    controls_changed = Signal()

    # component settings
    current_method_index: int
    current_method_index = ContextSetting(0)
    current_group_index: int
    current_group_index = ContextSetting(0)
    stored_selections: dict
    stored_selections = ContextSetting({})

    # default threshold defining expressed genes for Hypergeometric Test
    expression_threshold_value: float
    expression_threshold_value = ContextSetting(1.0)

    score_methods = [
        ('T-test', score_t_test),
        ('Mann-Whitney', score_mann_whitney),
        ('Hypergeometric Test', score_hypergeometric_test),
    ]

    feature_name = 'Feature name'

    def __init__(self, parent_widget, parent_component):
        QObject.__init__(self)
        OWComponent.__init__(self, widget=parent_widget)

        self.score_method_combo = comboBox(
            parent_component, self, 'current_method_index', label='Method', callback=self.__on_score_method_changed
        )
        self.score_method_combo.addItems([name for name, _ in self.score_methods])

        self.expression_threshold_box = widgetBox(parent_component, 'Expression threshold', margin=0)
        self.expression_threshold_box.setFlat(True)
        self.expression_threshold_spinbox = doubleSpin(
            self.expression_threshold_box,
            self,
            'expression_threshold_value',
            minv=0,
            maxv=1e2,
            step=1e-2,
            callback=self.__on_expression_threshold_changed,
            callbackOnReturn=True,
        )
        self.__show_expression_threshold_spinbox()

        self.group_combo = comboBox(
            parent_component, self, 'current_group_index', label='Group', callback=self.__on_group_index_changed
        )

        self.list_widget = QListWidget()
        self.list_widget.setSelectionMode(QListView.ExtendedSelection)
        self.list_widget.selectionModel().selectionChanged.connect(self.__on_selection_changed)

        box = widgetBox(parent_component, 'Values', margin=0)
        box.setFlat(True)
        box.layout().addWidget(self.list_widget)

        self.groups = {}
        self.data = None

    def initialize(self, data):
        """ Initialize widget state after receiving new data.
        """

        if data is not None:
            self.data = data

            column_groups, row_groups = self.group_candidates(data)
            self.groups = {index: value for index, value in enumerate(column_groups + row_groups)}

            self.group_combo.clear()
            self.group_combo.addItems([str(x.name) for x in self.groups.values()])
            self.group_combo.setCurrentIndex(self.current_group_index)

            self.__populate_list_widget()

    def get_selection_mask(self):
        """ Return the selection masks for the group.
        """
        if self.data is not None:
            group, indices = self.selected_split()

            if isinstance(group, column_group):
                selected = [group.values[i] for i in indices]
                target = {(group.key, value) for value in selected}
                _i = [
                    bool({*col.attributes.items(), (self.feature_name, col.name)}.intersection(target))
                    for col in self.data.domain.attributes
                ]
                return np.array(_i, dtype=bool)

            elif isinstance(group, row_group):
                target = set(indices)
                x, _ = self.data.get_column_view(group.var)
                _i = np.zeros_like(x, dtype=bool)
                for i in target:
                    _i |= x == i
                return _i
            else:
                raise TypeError("column_group or row_group expected, got {}".format(type(group).__name__))

    def selected_split(self):
        group_index = self.group_combo.currentIndex()
        if not (0 <= group_index < len(self.groups)):
            return None, []

        group = self.groups[group_index]
        selection = [model_idx.row() for model_idx in self.list_widget.selectedIndexes()]

        return group, selection

    def get_expression_threshold(self) -> float:
        return self.expression_threshold_value

    def get_score_method(self) -> Callable:
        _, method = self.score_methods[self.current_method_index]
        return method

    def __populate_list_widget(self) -> None:
        if self.list_widget is not None:
            self.list_widget.selectionModel().selectionChanged.disconnect(self.__on_selection_changed)
            self.list_widget.clear()

            target = self.groups.get(self.current_group_index, None)
            if target is not None:
                self.list_widget.addItems(target.values)

            self.list_widget.setSizeAdjustPolicy(QListWidget.AdjustToContents)

            self.__set_selection()
            self.list_widget.selectionModel().selectionChanged.connect(self.__on_selection_changed)

    def __show_expression_threshold_spinbox(self):
        self.expression_threshold_box.setHidden(self.get_score_method().__name__ != score_hypergeometric_test.__name__)

    def __on_expression_threshold_changed(self):
        self.controls_changed.emit()
        self.expression_threshold_changed.emit(self.expression_threshold_value)

    def __on_score_method_changed(self) -> None:
        self.__show_expression_threshold_spinbox()
        self.controls_changed.emit()
        self.score_method_changed.emit(self.current_method_index)

    def __on_group_index_changed(self) -> None:
        self.__populate_list_widget()
        self.controls_changed.emit()
        self.group_changed.emit(self.current_group_index)

    def __on_selection_changed(self) -> None:
        self.__store_selection()
        self.controls_changed.emit()
        self.selection_changed.emit(self.current_group_index)

    def __store_selection(self) -> None:
        self.stored_selections[self.current_group_index] = tuple(
            model_idx.row() for model_idx in self.list_widget.selectedIndexes()
        )

    def __set_selection(self) -> None:
        # Restore previous selection for root (if available)
        group_index = self.current_group_index
        indices = [
            self.list_widget.indexFromItem(self.list_widget.item(index))
            for index in self.stored_selections.get(group_index, [])
        ]
        selection = QItemSelection()
        for ind in indices:
            selection.select(ind, ind)
        self.list_widget.selectionModel().select(selection, QItemSelectionModel.ClearAndSelect)

    @staticmethod
    def group_candidates(data: Table) -> (column_group, row_group):
        genes_in_columns = data.attributes.get(TableAnnotation.gene_as_attr_name, None)

        # this should never happen. Widget should handle  table annotation warnings.
        if genes_in_columns is None:
            raise ValueError('Table annotation incorrect or missing.')

        domain = data.domain
        targets = defaultdict(set)

        for col, (label, value) in ((col, attrs) for col in domain.attributes for attrs in col.attributes.items()):
            targets[label].add(value)

            if not genes_in_columns:
                targets[GeneScoringComponent.feature_name].add(col.name)

        # Need at least 2 distinct values or key
        targets = [(key, sorted(vals)) for key, vals in targets.items() if len(vals) >= 2]
        column_groups = [column_group(key, key, values) for key, values in sorted(targets)]

        disc_vars = [
            var
            for var in data.domain.class_vars + data.domain.metas
            if isinstance(var, DiscreteVariable) and len(var.values) >= 2
        ]

        row_groups = [row_group(var.name, var, var.values) for var in disc_vars]
        return column_groups, row_groups


if __name__ == "__main__":
    from Orange.widgets.utils.widgetpreview import WidgetPreview
    from Orange.widgets.widget import OWWidget
    from Orange.data import Table
    from Orange.widgets.settings import SettingProvider

    class MockWidget(OWWidget):
        name = "Mock"
        scoring_component = SettingProvider(GeneScoringComponent)

        def __init__(self):
            self.scoring_component = GeneScoringComponent(self, self.mainArea)
            iris = Table('iris')
            iris.attributes[TableAnnotation.gene_as_attr_name] = False
            self.scoring_component.initialize(iris)

    WidgetPreview(MockWidget).run()
