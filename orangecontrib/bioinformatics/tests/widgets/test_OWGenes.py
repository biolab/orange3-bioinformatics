import os
import unittest

from orangewidget.tests.base import WidgetTest

from Orange.data import Table

from orangecontrib.bioinformatics.widgets.OWGenes import OWGenes
from orangecontrib.bioinformatics.widgets.utils.data import TAX_ID, GENE_ID_COLUMN, GENE_AS_ATTRIBUTE_NAME


class TestOWGenes(WidgetTest):
    def setUp(self):
        self._path = os.path.join(os.path.dirname(__file__), 'data')
        self.widget = self.create_widget(OWGenes)
        self.markers_file = 'genes.tab'

        file_name = os.path.join(self._path, self.markers_file)

        # gene names are stored in "Name" column.
        # we expect widget to detect this on its own.
        self.markers_column = 'Name'
        self.target_db = 'Entrez ID'

        # get input data
        self.send_signal(self.widget.Inputs.data_table, Table(file_name))
        self.assertIsNotNone(self.widget.input_data)

        # check if settings were set correctly
        self.assertEqual(self.markers_column, self.widget.selected_gene_col.name)
        self.assertFalse(self.widget.use_attr_names)

        # check if genes were extracted from the table
        self.assertIsNotNone(self.widget.input_genes)
        genes, _ = self.widget.input_data.get_column_view(self.widget.selected_gene_col)
        self.assertEqual(len(self.widget.input_genes), len(genes))

        self.wait_until_finished()
        # wait for gene matcher to finish
        # self.widget.threadpool.waitForDone()
        # self.process_events()

        # gene matcher should map all of the marker names into corresponding Entrez ID
        self.gene_match_result = self.widget.gene_matcher.get_known_genes()
        self.assertEqual(len(genes), len(self.gene_match_result))

    def test_output_data(self):
        # get output data
        out_data = self.get_output(self.widget.Outputs.data_table)

        # check if 'target_db' column exists. Target database is Entrez
        self.assertTrue(self.target_db in out_data.domain)

        # test if genes in the output data is the same as known genes from gene matcher.
        ids, _ = out_data.get_column_view(self.target_db)
        self.assertTrue((ids == [str(g.gene_id) for g in self.gene_match_result]).all())

        # test if table on the output is properly annotated
        table_annotations = out_data.attributes
        self.assertTrue(len(table_annotations) > 0)
        self.assertTrue(GENE_ID_COLUMN in table_annotations)
        self.assertTrue(TAX_ID in table_annotations)
        self.assertTrue(GENE_AS_ATTRIBUTE_NAME in table_annotations)
        self.assertFalse(table_annotations[GENE_AS_ATTRIBUTE_NAME])
        self.assertTrue(table_annotations[GENE_ID_COLUMN] == self.target_db)
        self.assertTrue(table_annotations[TAX_ID] == self.widget.tax_id)

    def test_filtered_output_data(self):
        self.widget.search_pattern = 'tubulin gamma'
        self.widget._apply_filter()

        # get output data
        out_data = self.get_output(self.widget.Outputs.data_table)

        # check if 'target_db' column exists. Target database is Entrez
        self.assertTrue(self.target_db in out_data.domain)

        # we expect filtered data table on the output
        ids, _ = out_data.get_column_view(self.target_db)
        self.assertTrue(len(ids) == 2)

        # test if table on the output is properly annotated
        table_annotations = out_data.attributes
        self.assertTrue(len(table_annotations) > 0)
        self.assertTrue(GENE_ID_COLUMN in table_annotations)
        self.assertTrue(TAX_ID in table_annotations)
        self.assertTrue(GENE_AS_ATTRIBUTE_NAME in table_annotations)
        self.assertFalse(table_annotations[GENE_AS_ATTRIBUTE_NAME])
        self.assertTrue(table_annotations[GENE_ID_COLUMN] == self.target_db)
        self.assertTrue(table_annotations[TAX_ID] == self.widget.tax_id)


if __name__ == '__main__':
    unittest.main()
