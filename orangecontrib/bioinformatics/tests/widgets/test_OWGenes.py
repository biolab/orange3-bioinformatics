import unittest
import os

from orangecontrib.bioinformatics.ncbi.gene import GENE_AS_ATTRIBUTE_NAME, GENE_ID_COLUMN, TAX_ID
from orangecontrib.bioinformatics.widgets.OWGenes import OWGenes
from Orange.widgets.tests.base import WidgetTest
from Orange.data import Table


class TestOWGenes(WidgetTest):

    def setUp(self):
        self._path = os.path.join(os.path.dirname(__file__), 'data')
        self.widget = self.create_widget(OWGenes)
        self.markers_file = 'markers.tab.gz'

    def test_input_markers(self):
        file_name = os.path.join(self._path, self.markers_file)

        # gene names are stored in "Name" attribute.
        # we expect widget to detect this on its own.
        markers_column = 'Name'
        target_db = 'Entrez ID'

        # get input data
        self.send_signal(self.widget.Inputs.data_table, Table(file_name))
        self.assertIsNotNone(self.widget.input_data)

        # check if settings were set correctly
        self.assertEqual(markers_column, self.widget.selected_gene_col.name)
        self.assertFalse(self.widget.use_attr_names)

        # check if genes were extracted from the table
        self.assertIsNotNone(self.widget.input_genes)
        genes, _ = self.widget.input_data.get_column_view(self.widget.selected_gene_col)
        self.assertEqual(len(self.widget.input_genes), len(genes))

        # wait for gene matcher to finish
        self.widget.threadpool.waitForDone()

        # gene matcher should map all of the marker names into corresponding Entrez ID
        gene_match_result = self.widget.gene_matcher.get_known_genes()
        self.assertEqual(len(genes), len(gene_match_result))

        # get output data
        self.widget.commit()
        out_data = self.get_output(self.widget.Outputs.data_table)

        # check if 'target_db' column exists. Target database is Entrez
        self.assertTrue(target_db in out_data.domain)

        # test if genes in the output data is the same as known genes from gene matcher.
        ids, _ = out_data.get_column_view(target_db)
        self.assertTrue((ids == [str(g.ncbi_id) for g in gene_match_result]).all())

        # test if table on the output is properly annotated
        table_annotations = out_data.attributes
        self.assertTrue(len(table_annotations) > 0)
        self.assertTrue(GENE_ID_COLUMN in table_annotations)
        self.assertTrue(TAX_ID in table_annotations)
        self.assertTrue(GENE_AS_ATTRIBUTE_NAME in table_annotations)
        self.assertFalse(table_annotations[GENE_AS_ATTRIBUTE_NAME])
        self.assertTrue(table_annotations[GENE_ID_COLUMN] == target_db)
        self.assertTrue(table_annotations[TAX_ID] == self.widget.tax_id)


if __name__ == '__main__':
    unittest.main()
