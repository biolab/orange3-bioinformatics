import unittest
from os.path import basename, normpath

from Orange.data import Table

from orangecontrib.bioinformatics.ncbi.gene import ENTREZ_ID, Gene, GeneInfo, GeneMatcher


class TestGene(unittest.TestCase):
    def test_load_attributes(self):
        g = Gene()
        g.load_attributes(('Human', '9606', '920', 'CD4'), attributes=('species', 'tax_id', 'gene_id', 'symbol'))
        self.assertEqual(g.species, 'Human')
        self.assertEqual(g.tax_id, '9606')
        self.assertEqual(g.gene_id, '920')
        self.assertEqual(g.symbol, 'CD4')
        self.assertEqual(str(g), '<Gene symbol=CD4, tax_id=9606, gene_id=920>')

        self.assertIsNone(g.input_identifier)
        self.assertIsNone(g.synonyms)

    def test_homologs(self):
        gm = GeneMatcher('9606')
        gm.genes = ['920']
        g = gm.genes[0]

        self.assertIsNotNone(g.homologs)
        self.assertTrue(len(g.homologs))
        self.assertIn('10090', g.homologs)
        self.assertEqual(g.homology_group_id, '513')

        self.assertEqual(g.homolog_gene('10090'), '12504')
        self.assertIsNone(g.homolog_gene('Unknown_taxonomy'))


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
        data = gm.match_table_attributes(data, rename=True, source_name='FooBar')

        for column in data.domain.attributes:
            self.assertTrue(ENTREZ_ID in column.attributes)
            self.assertTrue('FooBar' in column.attributes)


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
