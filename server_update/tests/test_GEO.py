from unittest import TestCase
from orangecontrib.bioinformatics.geo.dataset import GDSInfo


class GEOTest(TestCase):

    def setUp(self):
        self.info = GDSInfo()

    def test_GDS16(self):
        """ Test content of gds_info file. Check if GDS16 is correctly stored.
        """
        gds16 = self.info['GDS16']
        self.assertEqual(gds16['taxid'], '4932')
        self.assertEqual(gds16['feature_count'], 9216)
        self.assertEqual(gds16['sample_count'], 8)

    def test_size(self):
        self.assertGreater(len(self.info), 4000)
