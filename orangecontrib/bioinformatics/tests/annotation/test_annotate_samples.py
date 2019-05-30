import unittest

import numpy as np
from Orange.data import Table, Domain, StringVariable, ContinuousVariable

from orangecontrib.bioinformatics.annotation.annotate_samples import \
    AnnotateSamples
from orangecontrib.bioinformatics.widgets.utils.data import TAX_ID


class TestAnnotateSamples(unittest.TestCase):

    def setUp(self):
        m_domain = Domain(
            [], None, [
                StringVariable("Cell Type"), StringVariable("Entrez ID")])
        m_data = [["Type 1", "111"], ["Type 1", "112"], ["Type 1", "113"],
                  ["Type 2", "211"], ["Type 2", "212"], ["Type 2", "213"]]
        self.markers = Table(m_domain, np.empty((len(m_data), 0)), None, m_data)

        genes = ["111", "112", "113", "211", "212", "213"]
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
        self.data.attributes[TAX_ID] = "9606"  # id for a human

        self.annotator = AnnotateSamples()

    def test_artificial_data(self):
        annotations = self.annotator.annotate_samples(self.data, self.markers)

        self.assertEqual(type(annotations), Table)
        self.assertEqual(len(annotations), len(self.data))
        self.assertEqual(len(annotations[0]), 2)  # two types in the data
        self.assertGreater(annotations.X.sum(), 0)
        self.assertLessEqual(annotations.X.max(), 1)
        self.assertGreaterEqual(annotations.X.min(), 0)

    def test_remove_empty_column(self):
        m_domain = Domain(
            [], None, [
                StringVariable("Cell Type"), StringVariable("Entrez ID")])
        m_data = [["Type 1", "111"], ["Type 1", "112"], ["Type 1", "113"],
                  ["Type 2", "211"], ["Type 2", "212"], ["Type 2", "213"],
                  ["Type 3", "311"], ["Type 3", "312"], ["Type 3", "313"]]
        markers = Table(m_domain, np.empty((len(m_data), 0)), None, m_data)

        annotations = self.annotator.annotate_samples(self.data, markers)

        self.assertEqual(type(annotations), Table)
        self.assertEqual(len(annotations), len(self.data))
        self.assertEqual(len(annotations[0]), 2)  # two types in the data
        self.assertGreater(annotations.X.sum(), 0)
        self.assertLessEqual(annotations.X.max(), 1)
        self.assertGreaterEqual(annotations.X.min(), 0)

        annotations = self.annotator.annotate_samples(
            self.data, markers, return_nonzero_annotations=False)
        self.assertEqual(len(annotations), len(self.data))
        self.assertEqual(len(annotations[0]), 3)  # two types in the data
        self.assertGreater(annotations.X.sum(), 0)
        self.assertLessEqual(annotations.X.max(), 1)
        self.assertGreaterEqual(annotations.X.min(), 0)

    def test_sf(self):
        """
        Test annotations with hypergeom.sf
        """
        annotator = AnnotateSamples(p_value_fun="TEST_HYPERGEOMETRIC")
        annotations = annotator.annotate_samples(self.data, self.markers)

        self.assertEqual(type(annotations), Table)
        self.assertEqual(len(annotations), len(self.data))
        self.assertEqual(len(annotations[0]), 2)  # two types in the data
        self.assertGreater(annotations.X.sum(), 0)
        self.assertLessEqual(annotations.X.max(), 1)
        self.assertGreaterEqual(annotations.X.min(), 0)

    def test_two_example(self):
        self.data.X = self.data.X[:2]

        annotations = self.annotator.annotate_samples(self.data, self.markers)

        self.assertEqual(type(annotations), Table)
        self.assertEqual(len(annotations), len(self.data))

    def test_markers_without_entrez_id(self):
        self.markers[1, "Entrez ID"] = "?"
        annotations = self.annotator.annotate_samples(self.data, self.markers)

        self.assertEqual(type(annotations), Table)
        self.assertEqual(len(annotations), len(self.data))
        self.assertEqual(len(annotations[0]), 2)  # two types in the data
        self.assertGreater(annotations.X.sum(), 0)
        self.assertLessEqual(annotations.X.max(), 1)
        self.assertGreaterEqual(annotations.X.min(), 0)

    def test_select_attributes(self):
        selected_attrs, z = self.annotator.select_attributes(self.data)

        self.assertEqual(len(self.data), len(selected_attrs))
        self.assertSetEqual({"112", "111"}, selected_attrs[0])

    def test_assign_annotations(self):
        attrs = [{"111", "112", "113"}, {"111", "112", "211"},
                 {"211", "212", "111"}, {"211", "212", "213"}]
        exp_ann = np.array([
            [1, 0],
            [2/3, 1/3],
            [1/3, 2/3],
            [0, 1]])
        annotations, fdrs = self.annotator.assign_annotations(
            attrs, self.markers, "9606")

        self.assertEqual(len(attrs), len(annotations))
        self.assertEqual(len(attrs), len(fdrs))
        self.assertEqual(2, annotations.X.shape[1])  # only two types in markers
        self.assertEqual(2, fdrs.X.shape[1])
        np.testing.assert_array_almost_equal(exp_ann, annotations)

        exp_fdrs_smaller = np.array([
            [0.05, 2],
            [0.05, 2],
            [2, 0.05],
            [2, 0.05]])

        np.testing.assert_array_less(fdrs, exp_fdrs_smaller)
