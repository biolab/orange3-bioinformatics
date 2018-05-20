import unittest

from orangecontrib.bioinformatics.ncbi import gene


class GeneMatcher(unittest.TestCase):

    def test_types(self):
        with self.assertRaises(TypeError):
            gene.GeneMatcher(9606)

    def test_multiple_hits_scenario(self):
        input_gene_name = 'HB1'
        organism = '9606'

        gene_matcher = gene.GeneMatcher(organism)
        gene_matcher.genes = [input_gene_name]
        gene_matcher.run_matcher()
        result = gene_matcher.genes[0]

        self.assertEqual(result.input_name, input_gene_name)
        self.assertEqual(result.type_of_match, None)
        self.assertEqual(result.ncbi_id, None)
        self.assertGreater(len(result.possible_hits), 0)

    def test_symbol_match_scenario(self):
        input_gene_name = 'SCN5A'
        ncbi_id = 6331
        organism = '9606'

        gene_matcher = gene.GeneMatcher(organism)
        gene_matcher.genes = [input_gene_name]
        gene_matcher.run_matcher()
        result = gene_matcher.genes[0]

        self.assertEqual(result.input_name, input_gene_name)
        self.assertEqual(result.type_of_match, gene._symbol)
        self.assertEqual(result.ncbi_id, ncbi_id)

        result.load_ncbi_info()
        for tag in gene.GENE_INFO_TAGS:
            self.assertIsNotNone(getattr(result, tag))

    def test_case_insensitive(self):
        organism = '10090'
        original_name = 'Pou4f1'
        ncbi_id = 18996

        gene_matcher = gene.GeneMatcher(organism)
        gene_matcher.genes = ['Pou4F1', 'pou4F1', 'POU4F1', 'pou4f1', original_name]
        gene_matcher.run_matcher()

        self.assertEqual(len([g.ncbi_id for g in gene_matcher.genes if g.ncbi_id]), 1)

        gene_matcher = gene.GeneMatcher(organism, case_insensitive=True)
        gene_matcher.genes = ['Pou4F1', 'pou4F1', 'POU4F1', 'pou4f1', original_name]
        gene_matcher.run_matcher()

        self.assertEqual(len([g.ncbi_id for g in gene_matcher.genes if g.ncbi_id]), 5)
        self.assertEqual(set([g.ncbi_id for g in gene_matcher.genes if g.ncbi_id]).pop(), ncbi_id)


class GeneInfo(unittest.TestCase):

    def test_gene_info(self):
        ncbi_id = 6331
        organism = '9606'

        gene_info_obj = gene.GeneInfo(organism)
        gene_info = gene_info_obj.get_gene_by_id(ncbi_id)

        self.assertIsNone(gene_info.input_name)

        for tag in gene.GENE_INFO_TAGS:
            self.assertIsNotNone(getattr(gene_info, tag))


if __name__ == '__main__':
    unittest.main()
