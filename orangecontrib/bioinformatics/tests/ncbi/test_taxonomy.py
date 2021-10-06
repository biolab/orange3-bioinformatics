import unittest

from orangecontrib.bioinformatics.ncbi import taxonomy


class TestTaxonomy(unittest.TestCase):

    human = '9606'
    dicty = '44689'
    dog = '9615'

    def setUp(self) -> None:
        self.tax_obj = taxonomy.Taxonomy()
        self.assertGreater(len(self.tax_obj.taxids()), 1500000)

    def test_common_taxonomy(self):
        self.assertGreater(len(taxonomy.common_taxids()), 0)

        self.assertEqual(taxonomy.name(self.human), 'Homo sapiens')
        self.assertEqual(taxonomy.name(self.dicty), 'Dictyostelium discoideum')

        self.assertEqual(taxonomy.species_name_to_taxid('Homo sapiens'), self.human)
        self.assertEqual(taxonomy.species_name_to_taxid('Dictyostelium discoideum'), self.dicty)

        self.assertGreater(len(taxonomy.shortname(self.human)), 0)
        self.assertGreater(len(taxonomy.shortname(self.dicty)), 0)

    def test_uncommon_taxonomy(self):
        self.assertTrue(self.dog not in taxonomy.common_taxids())
        self.assertEqual(taxonomy.name(self.dog), 'Canis lupus familiaris')

        # not supported yet.
        self.assertIsNone(taxonomy.species_name_to_taxid('Canis lupus familiaris'))
        self.assertFalse(len(taxonomy.shortname(self.dog)))

    def test_human(self):
        self.assertTrue(isinstance(self.tax_obj.get_entry(self.human), taxonomy.utils.Taxon))

        self.assertRaises(taxonomy.utils.UnknownSpeciesIdentifier, self.tax_obj.get_entry, 'unknown_tax')

        self.assertTrue(('man', 'common name') in self.tax_obj.other_names(self.human))
        self.assertEqual(self.tax_obj.rank(self.human), 'species')
        self.assertEqual(self.tax_obj.parent(self.human), '9605')

        self.assertGreater(len(self.tax_obj.search('Homo sapiens', exact=True)), 0)
        self.assertGreater(len(self.tax_obj.lineage(self.human)), 0)
        self.assertGreater(len(self.tax_obj.get_all_strains(self.human)), 0)

        subnodes = self.tax_obj.subnodes(self.human)
        self.assertTrue(len(subnodes) >= 2)
        self.assertTrue('63221' in subnodes)
        self.assertTrue('741158' in subnodes)

        neanderthal = self.tax_obj.get_entry('63221')
        self.assertTrue(neanderthal.parent_tax_id == self.tax_obj.get_species('63221'))
        self.assertEqual(neanderthal.name, 'Homo sapiens neanderthalensis')

        denisovan = self.tax_obj.get_entry('741158')
        self.assertTrue(denisovan.parent_tax_id == self.tax_obj.get_species('741158'))
        self.assertEqual(denisovan.name, "Homo sapiens subsp. 'Denisova'")

        self.assertTrue(len(taxonomy.search('Homo sapiens', exact=True)) == 1)
        self.assertIn(self.human, taxonomy.search('Homo sapiens', exact=True))

        linage = taxonomy.lineage(self.human)
        self.assertEqual(linage[0], '1')
        self.assertEqual(linage[-1], self.tax_obj.parent(self.human))

        search_result = taxonomy.search('Homo sapiens')
        self.assertTrue(len(search_result))
        self.assertIn(self.human, search_result)

        # unclassified Mammalia: Homo sapiens x Mus musculus hybrid cell line
        self.assertIn('1131344', search_result)

        # subnodes are included
        self.assertIn(neanderthal.tax_id, taxonomy.search('Homo sapiens', only_species=False))
        self.assertIn(denisovan.tax_id, taxonomy.search('Homo sapiens', only_species=False))


if __name__ == '__main__':
    unittest.main()
