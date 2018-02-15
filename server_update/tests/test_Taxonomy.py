from unittest import TestCase
from orangecontrib.bioinformatics.ncbi import taxonomy


class TaxonomyTest(TestCase):

    def setUp(self):
        self.common_ids = taxonomy.common_taxids()
        self.organisms = [(taxonomy.name(tax_id), tax_id) for tax_id in self.common_ids]
        self.taxon = taxonomy.Taxonomy()

    def test_ids_count(self):
        self.assertGreater(len(self.taxon.taxids()), 1413000)

    def test_common_organisms(self):
        for id in self.common_ids:
            # create taxon from common organisms
            self.taxon.get_entry(id)
