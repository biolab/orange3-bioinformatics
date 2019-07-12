import unittest

from orangecontrib.bioinformatics.ncbi.gene import GeneMatcher, GeneInfo, Gene


class TestGeneMatcher(unittest.TestCase):

    def test_synonym_multiple_matches(self):
        gm = GeneMatcher('9606')
        gm.genes = ['HB1']
        gene = gm.genes[0]
        self.assertEqual(gene.input_identifier, 'HB1')
        # Gene matcher should not find any unique match
        self.assertEqual(gene.gene_id, None)

    def test_symbol_match_scenario(self):
        gm = GeneMatcher('9606')
        gm.genes = ['SCN5A']
        gene = gm.genes[0]

        self.assertEqual(gene.input_identifier, 'SCN5A')
        self.assertEqual(gene.symbol, 'SCN5A')
        self.assertEqual(gene.gene_id, '6331')

    def test_different_input_identifier_types(self):
        gm = GeneMatcher('9606')
        gm.genes = ['CD4', '614535', 'HB-1Y', 'ENSG00000205426']

        for gene in gm.genes:
            self.assertIsNotNone(gene.description)
            self.assertIsNotNone(gene.tax_id)
            self.assertIsNotNone(gene.species)
            self.assertIsNotNone(gene.gene_id)


class TestGeneInfo(unittest.TestCase):

    def test_gene_info(self):
        gi = GeneInfo('9606')
        gene = gi['6331']

        self.assertTrue('6331' in gi)
        self.assertIsInstance(gene, Gene)
        self.assertIsNotNone(gene.description)
        self.assertIsNotNone(gene.tax_id)
        self.assertIsNotNone(gene.species)
        self.assertIsNotNone(gene.gene_id)
        # must be None
        self.assertIsNone(gene.input_identifier)


if __name__ == '__main__':
    unittest.main()
