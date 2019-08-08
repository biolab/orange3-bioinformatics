import unittest

from os.path import normpath, basename

from orangecontrib.bioinformatics.ncbi.homologene import HomoloGene


class TestHomoloGene(unittest.TestCase):

    def setUp(self) -> None:
        self.homology = HomoloGene()
        self.assertEqual(basename(normpath(self.homology.file_path)), 'homologene.tab')

    def test_homology(self):
        self.assertEqual(self.homology.find_homolog('920', '9913').gene_id, '407098')
        self.assertEqual(self.homology.find_homolog('920', '10090').gene_id, '12504')
        self.assertEqual(self.homology.find_homolog('920', '10116').gene_id, '24932')
