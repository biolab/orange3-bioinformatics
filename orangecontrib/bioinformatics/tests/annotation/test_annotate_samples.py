import unittest

import numpy as np
from Orange.data import Table, Domain, StringVariable, ContinuousVariable

from orangecontrib.bioinformatics.annotation.annotate_samples import \
    AnnotateSamplesMeta, SCORING_EXP_RATIO, SCORING_MARKERS_SUM, \
    SCORING_LOG_FDR, \
    SCORING_LOG_PVALUE, PFUN_BINOMIAL, PFUN_HYPERGEOMETRIC
from orangecontrib.bioinformatics.widgets.utils.data import TAX_ID


class TestAnnotateSamples(unittest.TestCase):

    def setUp(self):
        m_domain = Domain(
            [], None, [
                StringVariable("Cell Type"), StringVariable("Entrez ID")])
        m_data = [["Type 1", "111"], ["Type 1", "112"], ["Type 1", "113"],
                  ["Type 1", "114"],
                  ["Type 2", "211"], ["Type 2", "212"], ["Type 2", "213"],
                  ["Type 2", "214"],
                  ]
        self.markers = Table(m_domain, np.empty((len(m_data), 0)), None, m_data)

        genes = ["111", "112", "113", "114", "211", "212", "213", "214"]
        self.domain = Domain([ContinuousVariable(
            str(g)) for g in genes])
        for v, g in zip(self.domain.attributes, genes):
            v.attributes = {"Entrez ID": g}
        self.data = Table(self.domain, np.array(
            [[1, 1, 1, 1.1, 0, 0, 0, 0],
             [1, .8, .9, 1, 0, 0, 0, 0],
             [.7, 1.1, 1, 1.2, 0, 0, 0, 0],
             [.8, .7, 1.1, 1, 0, .1, 0, 0],
             [0, 0, 0, 0, 1.05, 1.05, 1.1, 1],
             [0, 0, 0, 0, 1.1, 1.0, 1.05, 1.1],
             [0, 0, 0, 0, 1.05, .9, 1.1, 1.1],
             [0, 0, 0, 0, .9, .9, 1.2, 1]
             ]))
        self.data.attributes[TAX_ID] = "9606"  # id for a human

        self.annotator = AnnotateSamplesMeta()

    def test_mann_whitney_test(self):
        d = self.annotator.mann_whitney_test(self.data)
        self.assertEqual(type(d), Table)
        self.assertTupleEqual(self.data.X.shape, d.X.shape)

    def annotate_samples(
            self, data, markers, return_nonzero_annotations=True,
            scoring=SCORING_EXP_RATIO, p_value_fun=PFUN_BINOMIAL):
        """
        Helper method that performs all operations and returns final results.
        """
        z_values = self.annotator.mann_whitney_test(data)
        scores, fdrs = self.annotator.assign_annotations(
            z_values, markers, data, scoring=scoring, p_value_fun=p_value_fun)
        return self.annotator.filter_annotations(
            scores, fdrs, return_nonzero_annotations)

    def test_artificial_data(self):
        annotations = self.annotate_samples(self.data, self.markers)

        self.assertEqual(type(annotations), Table)
        self.assertEqual(len(annotations), len(self.data))
        self.assertEqual(len(annotations[0]), 2)  # two types in the data
        self.assertGreater(np.nansum(annotations.X), 0)
        self.assertLessEqual(np.nanmax(annotations.X), 1)
        self.assertGreaterEqual(np.nanmin(annotations.X), 0)

    def test_remove_empty_column(self):
        """
        Type 3 column must be removed here
        """
        m_domain = Domain(
            [], None, [
                StringVariable("Cell Type"), StringVariable("Entrez ID")])
        m_data = [["Type 1", "111"], ["Type 1", "112"], ["Type 1", "113"],
                  ["Type 1", "114"],
                  ["Type 2", "211"], ["Type 2", "212"], ["Type 2", "213"],
                  ["Type 2", "214"],
                  ["Type 3", "311"], ["Type 3", "312"], ["Type 3", "313"]]
        markers = Table(m_domain, np.empty((len(m_data), 0)), None, m_data)

        annotations = self.annotate_samples(self.data, markers)

        self.assertEqual(type(annotations), Table)
        self.assertEqual(len(annotations), len(self.data))
        self.assertEqual(len(annotations[0]), 2)  # two types in the data
        self.assertGreater(np.nansum(annotations.X), 0)
        self.assertLessEqual(np.nanmax(annotations.X), 1)
        self.assertGreaterEqual(np.nanmin(annotations.X), 0)

        annotations = self.annotate_samples(
            self.data, markers, return_nonzero_annotations=False)
        self.assertEqual(len(annotations), len(self.data))
        self.assertEqual(len(annotations[0]), 3)  # two types in the data
        self.assertGreater(np.nansum(annotations.X), 0)
        self.assertLessEqual(np.nanmax(annotations.X), 1)
        self.assertGreaterEqual(np.nanmin(annotations.X), 0)

    def test_sf(self):
        """
        Test annotations with hypergeom.sf
        """
        annotations = self.annotate_samples(
            self.data, self.markers, p_value_fun=PFUN_HYPERGEOMETRIC)

        self.assertEqual(type(annotations), Table)
        self.assertEqual(len(annotations), len(self.data))
        self.assertEqual(len(annotations[0]), 2)  # two types in the data
        self.assertGreater(np.nansum(annotations.X), 0)
        self.assertLessEqual(np.nanmax(annotations.X), 1)
        self.assertGreaterEqual(np.nanmin(annotations.X), 0)

    def test_two_example(self):
        self.data.X = self.data.X[:2]

        annotations = self.annotate_samples(self.data, self.markers)

        self.assertEqual(type(annotations), Table)
        self.assertEqual(len(annotations), len(self.data))

    def test_markers_without_entrez_id(self):
        self.markers[1, "Entrez ID"] = "?"
        annotations = self.annotate_samples(
            self.data, self.markers, return_nonzero_annotations=False)

        self.assertEqual(type(annotations), Table)
        self.assertEqual(len(annotations), len(self.data))
        self.assertEqual(len(annotations[0]), 2)  # two types in the data
        self.assertGreater(np.nansum(annotations.X), 0)
        self.assertLessEqual(np.nanmax(annotations.X), 1)
        self.assertGreaterEqual(np.nanmin(annotations.X), 0)

    def test_select_attributes(self):
        z = self.annotator.mann_whitney_test(self.data)

        self.assertEqual(z.X.shape, self.data.X.shape)
        self.assertGreaterEqual(z.X[0, 0], 1)
        self.assertGreaterEqual(z.X[0, 1], 1)
        self.assertGreaterEqual(z.X[0, 3], 1)

    def test_assign_annotations(self):
        z = np.array([
            [1.1, 1.1, 1.1, 1.1, 0, 0, 0, 0],
            [1.1, 1.1, 0, 0, 1.1, 0, 0, 0],
            [1.1, 0, 0, 0, 1.1, 1.1, 0, 0],
            [0, 0, 0, 0, 1.1, 1.1, 1.1, 1.1]
        ])
        z_table = Table(self.domain, z)
        attrs = [{"111", "112", "113", "114"}, {"111", "112", "211"},
                 {"211", "212", "111"}, {"211", "212", "213", "214"}]
        exp_ann = np.array([
            [1, 0],
            [1/2, 1/4],
            [1/4, 1/2],
            [0, 1]])
        annotations, fdrs = self.annotator.assign_annotations(
            z_table, self.markers, self.data)

        self.assertEqual(len(attrs), len(annotations))
        self.assertEqual(len(attrs), len(fdrs))
        self.assertEqual(2, annotations.X.shape[1])  # only two types in markers
        self.assertEqual(2, fdrs.X.shape[1])
        np.testing.assert_array_almost_equal(exp_ann, annotations)

        exp_fdrs_smaller = np.array([
            [0.05, 2],
            [2, 2],
            [2, 2],
            [2, 0.05]])

        np.testing.assert_array_less(fdrs, exp_fdrs_smaller)

    def test_scoring(self):
        # scoring SCORING_EXP_RATIO
        annotations = self.annotate_samples(
            self.data, self.markers, scoring=SCORING_EXP_RATIO)

        self.assertEqual(type(annotations), Table)
        self.assertEqual(len(annotations), len(self.data))
        self.assertEqual(len(annotations[0]), 2)  # two types in the data
        self.assertGreater(np.nansum(annotations.X), 0)
        self.assertLessEqual(np.nanmax(annotations.X), 1)
        self.assertGreaterEqual(np.nanmin(annotations.X), 0)

        # scoring SCORING_MARKERS_SUM
        annotations = self.annotate_samples(
            self.data, self.markers, scoring=SCORING_MARKERS_SUM)

        self.assertEqual(type(annotations), Table)
        self.assertEqual(len(annotations), len(self.data))
        self.assertEqual(len(annotations[0]), 2)  # two types in the data

        # based on provided data it should match
        # the third row is skipped, since it is special
        self.assertEqual(annotations[0, 0].value, self.data.X[0].sum())
        self.assertEqual(annotations[5, 1].value, self.data.X[5].sum())

        # scoring SCORING_LOG_FDR
        annotations = self.annotate_samples(
            self.data, self.markers, scoring=SCORING_LOG_FDR)

        self.assertEqual(type(annotations), Table)
        self.assertEqual(len(annotations), len(self.data))
        self.assertEqual(len(annotations[0]), 2)  # two types in the data

        # scoring SCORING_LOG_PVALUE
        annotations = self.annotate_samples(
            self.data, self.markers, scoring=SCORING_LOG_PVALUE)

        self.assertEqual(type(annotations), Table)
        self.assertEqual(len(annotations), len(self.data))
        self.assertEqual(len(annotations[0]), 2)  # two types in the data

    def test_entrez_id_not_string(self):
        """
        It seems that some datasets (e.g. AML dataset) have Entrez ID as int
        although they should be strings. Here we add the test for those cases.
        """
        # change Entrez IDs to int
        for i, att in enumerate(self.data.domain.attributes):
            att.attributes["Entrez ID"] = int(att.attributes["Entrez ID"])

        annotations = self.annotate_samples(self.data,
                                                      self.markers)

        self.assertEqual(type(annotations), Table)
        self.assertEqual(len(annotations), len(self.data))
        self.assertEqual(len(annotations[0]), 2)  # two types in the data
        self.assertGreater(np.nansum(annotations.X), 0)
        self.assertLessEqual(np.nanmax(annotations.X), 1)
        self.assertGreaterEqual(np.nanmin(annotations.X), 0)

