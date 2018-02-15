from unittest import TestCase
from orangecontrib.bioinformatics.dicty import DictyBase
from orangecontrib.bioinformatics.dicty import phenotypes as DictyMutants


class DictyBaseTest(TestCase):

    def setUp(self):
        self.dicty_base = DictyBase()

    def test_size(self):
        self.assertGreater(len(self.dicty_base.info), 14000)
        self.assertGreater(len(self.dicty_base.mappings), 27000)


class DictyMutantsTest(TestCase):

    def test_size(self):
        self.assertGreater(len(DictyMutants.gene_mutants()), 1000)
