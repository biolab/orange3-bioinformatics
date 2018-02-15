from unittest import TestCase
from server_update import sf_local
from orangecontrib.bioinformatics.geneset.utils import filename_parse


class GeneSetsTest(TestCase):

    def setUp(self):
        self.files = [fn for domain, fn in sf_local.listfiles('gene_sets') if fn != 'index.pck']

    def test_gene_sets(self):
        for file in self.files:
            hierarchy, organism = filename_parse(file)
            # TODO: fix tests
            # load_serverfiles(hierarchy, organism)
