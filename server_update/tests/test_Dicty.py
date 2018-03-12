from unittest import TestCase
from orangecontrib.bioinformatics.dicty import phenotypes as DictyMutants


class DictyMutantsTest(TestCase):

    def test_size(self):
        self.assertGreater(len(DictyMutants.gene_mutants()), 1000)
