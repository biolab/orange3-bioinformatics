import unittest
from os import makedirs, path

from Orange.data import Table

from orangecontrib.bioinformatics.geo.utils import (
    gds_is_cached, gds_download_url, gds_ensure_downloaded, gds_download,
    spots_mean, spots_median, spots_min, spots_max, unknown, parse_attribute_line, DOMAIN
)
from orangecontrib.bioinformatics.geo.dataset import GDSInfo, GDS
from orangecontrib.bioinformatics.utils import serverfiles


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
        self.assertEqual(gds_info[self.test_sample]['gene_count'], 9561)
        self.assertEqual(len(gds_info[self.test_sample]['samples']), 4)
        self.assertEqual(len(gds_info[self.test_sample]['subsets']), 2)

    def test_gds_utils(self):
        self.assertEqual(spots_mean([2, 2, 8, 8]), 5)
        self.assertEqual(spots_mean([8, unknown, 4, unknown]), 6)
        self.assertIsInstance(spots_mean([unknown, unknown]), float)

        self.assertEqual(spots_median([2, 2, 4, 8, 8]), 4)
        self.assertEqual(spots_median([8, unknown, 4, unknown]), 6)
        self.assertIsInstance(spots_median([unknown, unknown]), float)

        self.assertEqual(spots_min([2, 2, 8, 8]), 2)
        self.assertEqual(spots_min([8, unknown, 4, unknown]), 4)
        self.assertIsInstance(spots_min([unknown, unknown]), float)

        self.assertEqual(spots_max([2, 2, 8, 8]), 8)
        self.assertEqual(spots_max([8, unknown, 4, unknown]), 8)
        self.assertIsInstance(spots_max([unknown, unknown]), float)

        self.assertEqual(parse_attribute_line('!dataset_title = test_title'), ('title', 'test_title'))
        self.assertEqual(parse_attribute_line('!dataset_description = test_desc'),  ('description', 'test_desc'))
        self.assertEqual(parse_attribute_line('!dataset_gds_type = test_type'), ('gds_type', 'test_type'))

        # header like line
        self.assertRaises(AttributeError, parse_attribute_line, 'ID  SAMPLE1 SAMPLE2 SAMPLE3')

    def test_gds_data(self):
        # test url
        self.assertIsNotNone(gds_download_url(self.test_sample))

        # file not in cache
        self.assertFalse(gds_is_cached(self.test_sample))

        # download gds from serverfiles
        try:
            makedirs(serverfiles.localpath(DOMAIN))
        except OSError:
            if path.exists(serverfiles.localpath(DOMAIN)):
                pass
            else:
                # There was an error on creation, so make sure we know about it
                raise
        gds_download(self.test_sample)

        # file in cache
        self.assertIsNone(gds_ensure_downloaded(self.test_sample))
        self.assertTrue(gds_is_cached(self.test_sample))

        gds = GDS(self.test_sample)
        self.assertIsNotNone(gds.info)
        self.assertEqual(gds.info['gene_count'], 9561)
        self.assertEqual(len(gds.info['samples']), 4)
        self.assertEqual(len(gds.info['subsets']), 2)

        self.assertEqual(gds.info['taxid'], self.test_organism)

        self.assertIsInstance(gds.get_data(), Table)
        self.assertIsInstance(gds.get_data(transpose=True), Table)


if __name__ == '__main__':
    unittest.main()
