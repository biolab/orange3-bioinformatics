import unittest

from orangewidget.tests.base import WidgetTest

from Orange.data import Table

from orangecontrib.bioinformatics.widgets.utils.data import TableAnnotation
from orangecontrib.bioinformatics.widgets.OWGEODatasets import OWGEODatasets


class TestOWGEODatasets(WidgetTest):
    def setUp(self):
        self.test_sample = 'GDS1001'
        self.test_organism = '10090'
        self.widget = self.create_widget(OWGEODatasets)

    def test_minimum_size(self):
        pass

    def test_output_data(self):
        self.widget.auto_commit = True

        gds_info = self.widget.gds_info

        # check if GDSInfo loaded successfully
        self.assertIsNotNone(gds_info.get)

        # set selected gds
        self.widget.selected_gds = gds_info.get(self.test_sample)
        self.widget._set_selection()
        self.widget.on_gds_selection_changed()

        output_data = self.get_output(self.widget.Outputs.gds_data, wait=20000)

        self.assertIsInstance(output_data, Table)
        # test if table on the output is properly annotated
        table_annotations = output_data.attributes

        self.assertTrue(len(table_annotations) > 0)
        self.assertTrue(TableAnnotation.tax_id in table_annotations)
        self.assertTrue(TableAnnotation.gene_as_attr_name in table_annotations)

        # by default genes are in columns not in rows
        self.assertTrue(TableAnnotation.gene_id_attribute in table_annotations)
        # check if taxonomy is correct
        self.assertTrue(table_annotations[TableAnnotation.tax_id] == self.test_organism)

    def test_filter_no_match(self):
        filter_input = 'this will not match'
        self.widget.filter.setText(filter_input)
        self.assertEqual(self.widget.table_model._table.size, 0)
        self.assertEqual(self.widget.search_pattern, filter_input)

    def test_filter_match(self):
        filter_input = 'Embryonic stem cell-derived cardiomyocytes'
        gds_id = 'GDS3513'

        # set filter
        self.widget.filter.setText(filter_input)
        self.assertEqual(len(self.widget.table_model._table), 1)
        self.assertEqual(self.widget.search_pattern, filter_input)
        self.assertEqual(self.widget.table_model.get_row_index(gds_id), 0)
        self.assertTrue(gds_id in self.widget.table_model._table)

        # remove filter
        self.widget.filter.clear()
        self.assertEqual(self.widget.table_model._table.shape, self.widget.table_model.table.shape)


if __name__ == '__main__':
    unittest.main()
