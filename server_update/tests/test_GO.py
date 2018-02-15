from unittest import TestCase

from server_update import sf_local
from orangecontrib.bioinformatics.go import Ontology, Annotations


class GOTest(TestCase):

    def setUp(self):
        print(sf_local.listfiles('GO'))
        self.organisms = [fname.split('.')[1] for domain,temp, fname in sf_local.listfiles('GO')
                          if fname.startswith('gene_association')]
        self.onto = Ontology()

    def test_ontology(self):
        self.assertGreater(len(self.onto), 46000)

    def test_associations(self):
        for organism in self.organisms:
            Annotations(organism, ontology=self.onto)
