import unittest

from os.path import normpath, basename


from Orange.data import Table
from orangecontrib.bioinformatics.ncbi.gene import GeneMatcher, GeneInfo, Gene, ENTREZ_ID


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

    def test_taxonomy_change(self):
        gm = GeneMatcher('4932')
        self.assertEqual(gm.tax_id, '4932')
        self.assertEqual(basename(normpath(gm.gene_db_path)), '4932.sqlite')

        gm.tax_id = '9606'
        self.assertEqual(gm.tax_id, '9606')
        self.assertEqual(basename(normpath(gm.gene_db_path)), '9606.sqlite')

    def test_match_table_column(self):
        gm = GeneMatcher('4932')

        data = gm.match_table_column(Table('brown-selected.tab'), 'gene')
        self.assertTrue(ENTREZ_ID in data.domain)

    def test_match_table_attributes(self):
        gm = GeneMatcher('4932')

        data = Table('brown-selected.tab')
        data = Table.transpose(data, feature_names_column='gene')
        gm.match_table_attributes(data)

        for column in data.domain.attributes:
            self.assertTrue(ENTREZ_ID in column.attributes)


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
