import unittest

import numpy as np

from AnyQt.QtCore import QItemSelection, QItemSelectionModel
from AnyQt.QtTest import QSignalSpy

from Orange.data import Table
from Orange.preprocess import Remove
from Orange.widgets.widget import OWWidget
from Orange.widgets.settings import SettingProvider
from Orange.widgets.tests.base import WidgetTest
from Orange.widgets.tests.utils import simulate

from orangecontrib.bioinformatics.utils.statistics import score_hypergeometric_test
from orangecontrib.bioinformatics.widgets.components import GeneScoringComponent
from orangecontrib.bioinformatics.widgets.utils.data import TableAnnotation


class MockWidget(OWWidget):
    name = "Mock"
    scoring_component = SettingProvider(GeneScoringComponent)

    def __init__(self):
        self.scoring_component = GeneScoringComponent(self, self.mainArea)


def iris_test_case(data: Table):
    class TestGeneScoringComponent(WidgetTest):
        def setUp(self):
            self.widget = MockWidget()
            self.component = self.widget.scoring_component

        def test_expression_threshold_spinbox(self):
            # find index of item in combobox for hypergeometric test
            method_index, *_ = [
                index
                for index, (name, method) in enumerate(self.component.score_methods)
                if method == score_hypergeometric_test
            ]

            # check if spinbox appears after hypergeometric test is selected
            self.assertTrue(self.component.expression_threshold_box.isHidden())
            simulate.combobox_activate_index(self.component.score_method_combo, method_index)
            self.assertFalse(self.component.expression_threshold_box.isHidden())

        def test_scoring_methods_combobox(self):
            combo_box_values = [
                self.component.score_method_combo.itemText(i) for i in range(self.component.score_method_combo.count())
            ]
            self.assertTrue(len(combo_box_values) > 0)
            self.assertEqual([name for name, _ in self.component.score_methods], combo_box_values)

            signals_cb_emits = QSignalSpy(self.component.score_method_changed)
            simulate.combobox_run_through_all(self.component.score_method_combo)

            self.assertEqual(self.component.score_method_combo.currentIndex(), self.component.current_method_index)
            self.assertEqual(self.component.current_method_index, len(combo_box_values) - 1)

            # number of signals combobox emits should be equal to the length of available scoring methods
            self.assertEqual(len(combo_box_values), len(signals_cb_emits))

        def test_selected_group_values(self):
            self.assertIsNone(self.component.data)
            self.component.initialize(data)
            self.assertIsNotNone(self.component.data)

            # we expect only one value 'iris'
            combo_box_value, *_ = [
                self.component.group_combo.itemText(i) for i in range(self.component.group_combo.count())
            ]
            self.assertEqual(combo_box_value, 'iris')

            group_values = [
                self.component.list_widget.item(i).text() for i in range(self.component.list_widget.count())
            ]
            self.assertEqual(group_values, ['Iris-setosa', 'Iris-versicolor', 'Iris-virginica'])

        def test_selection(self):
            self.component.initialize(data)
            list_widget = self.component.list_widget
            signals_cb_emits = QSignalSpy(self.component.selection_changed)

            # get modelIndex from list widget
            iris_setosa_index = list_widget.indexFromItem(list_widget.item(0))
            iris_versicolor_index = list_widget.indexFromItem(list_widget.item(1))

            # set selection
            selection = QItemSelection()
            selection.select(iris_setosa_index, iris_setosa_index)
            selection.select(iris_versicolor_index, iris_versicolor_index)
            list_widget.selectionModel().select(selection, QItemSelectionModel.ClearAndSelect)

            # test if correct number of signals is emited
            self.assertEqual(1, len(signals_cb_emits))
            # test if selection is OK
            self.assertEqual(2, len(list_widget.selectedItems()))

            selection_mask = self.component.get_selection_mask()
            _selection_mask = ~selection_mask

            self.assertIsInstance(selection_mask, np.ndarray)
            self.assertTrue(selection_mask.dtype == np.bool)

            if 'iris' in data.domain:
                # test selection mask
                self.assertEqual(data.X[selection_mask, :].shape, (100, 4))

                remover = Remove(class_flags=Remove.RemoveUnusedValues)
                x1, x2 = remover(data[selection_mask, :]), remover(data[_selection_mask, :])

                selected_row_values = x1.domain['iris'].values
                unselected_row_values = x2.domain['iris'].values

                self.assertTrue(len(selected_row_values) == 2)
                self.assertIn('Iris-setosa', selected_row_values)
                self.assertIn('Iris-versicolor', selected_row_values)
                self.assertNotIn('Iris-virginica', selected_row_values)

                self.assertTrue(len(unselected_row_values) == 1)
                self.assertIn('Iris-virginica', unselected_row_values)
            else:
                # test selection mask
                x = data.X
                self.assertEqual(x[:, selection_mask].shape, (4, 100))
                self.assertEqual(x[:, _selection_mask].shape, (4, 50))

                selected_col_values = {
                    col.attributes.get('iris')
                    for col, selected in zip(data.domain.variables, selection_mask)
                    if selected
                }
                unselected_col_values = {
                    col.attributes.get('iris')
                    for col, selected in zip(data.domain.variables, _selection_mask)
                    if selected
                }

                self.assertTrue(len(selected_col_values) == 2)
                self.assertIn('Iris-setosa', selected_col_values)
                self.assertIn('Iris-versicolor', selected_col_values)
                self.assertNotIn('Iris-virginica', selected_col_values)

                self.assertTrue(len(unselected_col_values) == 1)
                self.assertIn('Iris-virginica', unselected_col_values)

    return TestGeneScoringComponent


iris = Table('iris')
iris.attributes[TableAnnotation.gene_as_attr_name] = False

iris_transposed = Table.transpose(Table('iris'))
iris_transposed.attributes[TableAnnotation.gene_as_attr_name] = True


class TestRowGroup(iris_test_case(iris)):
    pass


class TestColumnGroup(iris_test_case(iris_transposed)):
    pass


if __name__ == "__main__":
    unittest.main()
