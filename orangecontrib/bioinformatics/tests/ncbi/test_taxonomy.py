import unittest

from orangecontrib.bioinformatics.ncbi import taxonomy


class TestTaxonomy(unittest.TestCase):

    human = '9606'
    dicty = '352472'

    def test_common_taxonomy(self):
        self.assertGreater(len(taxonomy.common_taxids()), 0)

        self.assertEqual(taxonomy.name(self.human), 'Homo sapiens')
        self.assertEqual(taxonomy.name(self.dicty), 'Dictyostelium discoideum AX4')

        self.assertEqual(taxonomy.taxname_to_taxid('Homo sapiens'), self.human)
        self.assertEqual(taxonomy.taxname_to_taxid('Dictyostelium discoideum AX4'), self.dicty)

        self.assertGreater(len(taxonomy.shortname(self.human)), 0)
        self.assertGreater(len(taxonomy.shortname(self.dicty)), 0)

    def test_human(self):
        tax_obj = taxonomy.Taxonomy()
        self.assertGreater(len(tax_obj.taxids()), 1500000)

        self.assertTrue(isinstance(tax_obj.get_entry(self.human), taxonomy.utils.Taxon))
        self.assertRaises(taxonomy.utils.UnknownSpeciesIdentifier, tax_obj.get_entry, 'unknown_tax')

        self.assertTrue(('man', 'common name') in tax_obj.other_names(self.human))
        self.assertEqual(tax_obj.rank(self.human), 'species')
        self.assertEqual(tax_obj.parent(self.human), '9605')
        self.assertGreater(len(tax_obj.search('Homo sapiens', exact=True)), 0)
        self.assertGreater(len((tax_obj.lineage(self.human))), 0)

        self.assertGreater(len(tax_obj.get_all_strains(self.human)), 0)


if __name__ == '__main__':
    unittest.main()
