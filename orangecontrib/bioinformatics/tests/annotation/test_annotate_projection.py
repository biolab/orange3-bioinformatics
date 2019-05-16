import unittest

import numpy as np
from Orange.clustering import DBSCAN
from Orange.data import Domain, ContinuousVariable, Table

from orangecontrib.bioinformatics.annotation.annotate_projection import annotate_projection


class TestAnnotateSamples(unittest.TestCase):

    def setUp(self):
        domain = Domain([ContinuousVariable("x"), ContinuousVariable("y")])
        data = np.array([[1., 1.], [1., 1.5], [1.5, 1.], [1.2, 1.2], [1.7, 1.8],
                         [10., 10.], [10., 10.5], [10.5, 11.], [10.5, 10.7],
                         [10.5, 10.6]])
        self.data = Table(domain, data)

        domain_ann = Domain([ContinuousVariable("a"), ContinuousVariable("b"),
                            ContinuousVariable("c")])
        data_ann = np.array([[1, 0.5, 0.4], [0.9, 0.4, 0.2], [0.9, 0.4, 0.2],
                             [0.9, 0.4, 0.2], [0.9, 0.4, 0.2],
                             [0.2, 0.4, 0.7], [0.1, 0.4, 0.6], [0.2, 0.9, 0.2],
                             [0.1, 0.4, 0.6], [0.1, 0.8, 0.6]])
        self.annotations = Table(domain_ann, data_ann)

    def test_annotate_projection(self):
        clusters, annotations_cl = annotate_projection(
            self.annotations, self.data, clustering_algorithm=DBSCAN(eps=2))
        np.testing.assert_array_equal(
            clusters, np.array([0, 0, 0, 0, 0, 1, 1, 1, 1, 1]))

        self.assertEqual(len(annotations_cl), 2)
        self.assertEqual(len(annotations_cl[0]), 1)
        self.assertEqual(len(annotations_cl[1]), 2)
        self.assertEqual(annotations_cl[1][0][0], 'c')
        self.assertAlmostEqual(annotations_cl[1][0][1], 0.6, 5)
        self.assertEqual(annotations_cl[1][1][0], 'b')
        self.assertAlmostEqual(annotations_cl[1][1][1], 0.4, 5)

    def test_example_not_clustered(self):
        self.data[-1] = [23, 23]
        clusters, annotations_cl = annotate_projection(
            self.annotations, self.data, clustering_algorithm=DBSCAN(
                eps=2, min_samples=3))
        np.testing.assert_array_equal(
            clusters, np.array([0, 0, 0, 0, 0, 1, 1, 1, 1, -1]))

        self.assertEqual(len(annotations_cl), 2)
        self.assertEqual(len(annotations_cl[0]), 1)
        self.assertEqual(len(annotations_cl[1]), 2)
        self.assertEqual(annotations_cl[1][0][0], 'c')
        self.assertAlmostEqual(annotations_cl[1][0][1], 0.75, 5)
        self.assertEqual(annotations_cl[1][1][0], 'b')
        self.assertAlmostEqual(annotations_cl[1][1][1], 0.25, 5)
