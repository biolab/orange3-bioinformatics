import unittest

import numpy as np
from Orange.data import Table, Domain, StringVariable, ContinuousVariable
from scipy.stats import hypergeom

from orangecontrib.bioinformatics.annotation.annotate_samples import \
    AnnotateSamples


class TestAnnotateSamples(unittest.TestCase):

    def setUp(self):
        m_domain = Domain(
            [], None, [
                StringVariable("Cell Type"), StringVariable("Entrez ID")])
        m_data = [["Type 1", 111], ["Type 1", 112], ["Type 1", 113],
                  ["Type 2", 211], ["Type 2", 212], ["Type 2", 213]]
        self.markers = Table(m_domain, np.empty((len(m_data), 0)), None, m_data)

        genes = [111, 112, 113, 211, 212, 213]
        domain = Domain([ContinuousVariable(
            str(g)) for g in genes])
        for v, g in zip(domain.attributes, genes):
            v.attributes = {"Entrez ID": g}
        self.data = Table(domain, np.array(
            [[1, 1, 1, 0, 0, 0],
             [1, .8, .9, 0, 0, 0],
             [.7, 1.1, 1, 0, 0, 0],
             [.8, .7, 1.1, 0, .1, 0],
             [0, 0, 0, 1.05, 1.05, 1.1],
             [0, 0, 0, 1.1, 1.0, 1.05],
             [0, 0, 0, 1.05, .9, 1.1],
             [0, 0, 0, .9, .9, 1.2]
             ]))

    def test_artificial_data(self):
        annotator = AnnotateSamples(p_value_th=0.05)
        annotations = annotator.annotate_samples(self.data, self.markers)

        self.assertEqual(type(annotations), Table)
        self.assertEqual(len(annotations), len(self.data))
        self.assertEqual(len(annotations[0]), 2)  # two types in the data
        self.assertGreater(annotations.X.sum(), 0)
        self.assertLessEqual(annotations.X.max(), 1)
        self.assertGreaterEqual(annotations.X.min(), 0)

    def test_sf(self):
        """
        Test annotations with hypergeom.sf
        """
        annotator = AnnotateSamples(p_value_th=0.05, p_value_fun=hypergeom.sf)
        annotations = annotator.annotate_samples(self.data, self.markers)

        self.assertEqual(type(annotations), Table)
        self.assertEqual(len(annotations), len(self.data))
        self.assertEqual(len(annotations[0]), 2)  # two types in the data
        self.assertGreater(annotations.X.sum(), 0)
        self.assertLessEqual(annotations.X.max(), 1)
        self.assertGreaterEqual(annotations.X.min(), 0)


    def test_two_example(self):
        self.data.X = self.data.X[:2]

        annotator = AnnotateSamples(p_value_th=0.05)
        annotations = annotator.annotate_samples(self.data, self.markers)

        self.assertEqual(type(annotations), Table)
        self.assertEqual(len(annotations), len(self.data))
