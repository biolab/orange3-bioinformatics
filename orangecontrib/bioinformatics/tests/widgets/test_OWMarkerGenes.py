import unittest

from Orange.widgets.tests.base import WidgetTest
from orangecontrib.bioinformatics.widgets.OWMarkerGenes import OWMarkerGenes


class TestOWMarkerGenes(WidgetTest):
    def setUp(self):
        self.widget = self.create_widget(OWMarkerGenes)

        self.cell_marker_db = 'CellMarker'
        self.panglao_db = 'PanglaoDB'

    def test_data_not_empty(self):
        w = self.widget
        self.assertIsNotNone(w.data)

    def test_change_db_source(self):
        w = self.widget

        self.assertTrue(w.db_source_index == 0)
        self.assertTrue(w.selected_db_source == self.cell_marker_db)

        w.handle_source_changed(1)


        print(w.selected_db_source)
        print(self.panglao_db)


        self.assertTrue(w.db_source_index == 1)
        self.assertTrue(w.selected_db_source == self.panglao_db)

        self.assertIsNotNone(w.data)


if __name__ == '__main__':
    unittest.main()
