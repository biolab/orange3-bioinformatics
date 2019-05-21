import unittest

from orangecontrib.bioinformatics.ncbi.gene import OrangeTableAnnotations
from orangecontrib.bioinformatics.widgets.OWGEODatasets import OWGEODatasets
from Orange.widgets.tests.base import WidgetTest
from Orange.data import Table


class TestOWGEODatasets(WidgetTest):

    def setUp(self):
        self.test_sample = 'GDS1001'
        self.test_organism = '10090'
        self.widget = self.create_widget(OWGEODatasets)

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
        self.assertTrue(OrangeTableAnnotations.tax_id in table_annotations)
        self.assertTrue(OrangeTableAnnotations.gene_as_attribute_name in table_annotations)

        # by default genes are in columns not in rows
        self.assertTrue(OrangeTableAnnotations.gene_id_attribute in table_annotations)
        # check if taxonomy is correct
        self.assertTrue(table_annotations[OrangeTableAnnotations.tax_id] == self.test_organism)


if __name__ == '__main__':
    unittest.main()
