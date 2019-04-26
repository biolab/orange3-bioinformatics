""" GEO module tests """
import unittest

from Orange.data import Table

from orangecontrib.bioinformatics.geo.dataset import GDSInfo, GDS
from orangecontrib.bioinformatics.ncbi.gene import OrangeTableAnnotations


class TestGEO(unittest.TestCase):
    test_sample = 'GDS1001'
    test_organism = '10090'

    def test_gds_info(self):
        gds_info = GDSInfo()
        self.assertIsNotNone(gds_info)
        self.assertGreater(len(gds_info.keys()), 0)
        self.assertGreater(len(gds_info.items()), 0)
        self.assertGreater(len(gds_info.values()), 0)

        self.assertIsNotNone(gds_info[self.test_sample])
        self.assertEqual(gds_info[self.test_sample]['genes'], 9561)
        self.assertEqual(int(gds_info[self.test_sample]['sample_count']), 4)
        self.assertEqual(len(gds_info[self.test_sample]['subsets']), 2)

    def test_gds_data(self):
        gds_info = GDSInfo()
        gds_table = GDS(self.test_sample)

        # test if data is downloaded
        self.assertIsNotNone(gds_table)
        self.assertIsInstance(gds_table, Table)

        # test data table values
        rows, columns = gds_table.X.shape
        self.assertEqual(int(gds_info[self.test_sample]['sample_count']), rows)
        self.assertEqual(int(gds_info[self.test_sample]['genes']), columns)

        # test data table annotations
        self.assertTrue(OrangeTableAnnotations.gene_as_attribute_name in gds_table.attributes)
        self.assertTrue(OrangeTableAnnotations.gene_id_attribute in gds_table.attributes)
        self.assertTrue(OrangeTableAnnotations.tax_id in gds_table.attributes)

        self.assertTrue(gds_table.attributes[OrangeTableAnnotations.gene_as_attribute_name])

    def test_gds_data_transposed(self):
        gds_info = GDSInfo()
        gds_table = GDS(self.test_sample, transpose=True)

        # test if data is downloaded
        self.assertIsNotNone(gds_table)
        self.assertIsInstance(gds_table, Table)

        # test data table values
        rows, columns = gds_table.X.shape
        self.assertEqual(int(gds_info[self.test_sample]['sample_count']), columns)
        self.assertEqual(int(gds_info[self.test_sample]['genes']), rows)

        # test data table annotations
        self.assertTrue(OrangeTableAnnotations.gene_as_attribute_name in gds_table.attributes)
        self.assertTrue(OrangeTableAnnotations.gene_id_column in gds_table.attributes)
        self.assertTrue(OrangeTableAnnotations.tax_id in gds_table.attributes)

        self.assertFalse(gds_table.attributes[OrangeTableAnnotations.gene_as_attribute_name])


if __name__ == '__main__':
    unittest.main()
