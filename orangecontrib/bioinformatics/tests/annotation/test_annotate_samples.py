import unittest
import functools

import numpy as np
from pointannotator.annotate_samples import (
    SCORING_LOG_FDR,
    SCORING_LOG_PVALUE,
    PFUN_HYPERGEOMETRIC,
    SCORING_MARKERS_SUM,
)

from Orange.data import Table, Domain, StringVariable, DiscreteVariable, ContinuousVariable
from Orange.tests.test_statistics import dense_sparse

from orangecontrib.bioinformatics.widgets.utils.data import TAX_ID
from orangecontrib.bioinformatics.annotation.annotate_samples import (
    PFUN_BINOMIAL,
    SCORING_EXP_RATIO,
    AnnotateSamplesMeta,
)


class TestAnnotateSamples(unittest.TestCase):
    def setUp(self):
        m_domain = Domain([], None, [StringVariable("Cell Type"), StringVariable("Entrez ID")])
        m_data = [
            ["Type 1", "111"],
            ["Type 1", "112"],
            ["Type 1", "113"],
            ["Type 1", "114"],
            ["Type 2", "211"],
            ["Type 2", "212"],
            ["Type 2", "213"],
            ["Type 2", "214"],
        ]
        self.markers = Table(m_domain, np.empty((len(m_data), 0)), None, m_data)

        genes = ["111", "112", "113", "114", "211", "212", "213", "214"]
        self.domain = Domain([ContinuousVariable(str(g)) for g in genes])
        for v, g in zip(self.domain.attributes, genes):
            v.attributes = {"Entrez ID": g}
        self.data = Table(
            self.domain,
            np.array(
                [
                    [1, 1, 1, 1.1, 0, 0, 0, 0],
                    [1, 0.8, 0.9, 1, 0, 0, 0, 0],
                    [0.7, 1.1, 1, 1.2, 0, 0, 0, 0],
                    [0.8, 0.7, 1.1, 1, 0, 0.1, 0, 0],
                    [0, 0, 0, 0, 1.05, 1.05, 1.1, 1],
                    [0, 0, 0, 0, 1.1, 1.0, 1.05, 1.1],
                    [0, 0, 0, 0, 1.05, 0.9, 1.1, 1.1],
                    [0, 0, 0, 0, 0.9, 0.9, 1.2, 1],
                ]
            ),
        )
        self.data.attributes[TAX_ID] = "9606"  # id for a human

        self.annotator = AnnotateSamplesMeta()

    @dense_sparse
    def test_mann_whitney_test(self, array):
        self.data.X = array(self.data.X)
        d = self.annotator.mann_whitney_test(self.data)
        self.assertEqual(type(d), Table)
        self.assertTupleEqual(self.data.X.shape, d.X.shape)

    def annotate_samples(
        self, data, markers, return_nonzero_annotations=True, scoring=SCORING_EXP_RATIO, p_value_fun=PFUN_BINOMIAL
    ):
        """
        Helper method that performs all operations and returns final results.
        """
        z_values = self.annotator.mann_whitney_test(data)
        scores, fdrs = self.annotator.assign_annotations(
            z_values, markers, data, scoring=scoring, p_value_fun=p_value_fun
        )
        return self.annotator.filter_annotations(scores, fdrs, return_nonzero_annotations)

    @dense_sparse
    def test_artificial_data(self, array):
        self.data.X = array(self.data.X)
        annotations = self.annotate_samples(self.data, self.markers)

        self.assertEqual(type(annotations), Table)
        self.assertEqual(len(annotations), len(self.data))
        self.assertEqual(len(annotations[0]), 2)  # two types in the data
        self.assertGreater(np.nansum(annotations.X), 0)
        self.assertLessEqual(np.nanmax(annotations.X), 1)
        self.assertGreaterEqual(np.nanmin(annotations.X), 0)

    @dense_sparse
    def test_remove_empty_column(self, array):
        """
        Type 3 column must be removed here
        """
        m_domain = Domain([], None, [StringVariable("Cell Type"), StringVariable("Entrez ID")])
        m_data = [
            ["Type 1", "111"],
            ["Type 1", "112"],
            ["Type 1", "113"],
            ["Type 1", "114"],
            ["Type 2", "211"],
            ["Type 2", "212"],
            ["Type 2", "213"],
            ["Type 2", "214"],
            ["Type 3", "311"],
            ["Type 3", "312"],
            ["Type 3", "313"],
        ]
        markers = Table(m_domain, np.empty((len(m_data), 0)), None, m_data)

        self.data.X = array(self.data.X)
        annotations = self.annotate_samples(self.data, markers)

        self.assertEqual(type(annotations), Table)
        self.assertEqual(len(annotations), len(self.data))
        self.assertEqual(len(annotations[0]), 2)  # two types in the data
        self.assertGreater(np.nansum(annotations.X), 0)
        self.assertLessEqual(np.nanmax(annotations.X), 1)
        self.assertGreaterEqual(np.nanmin(annotations.X), 0)

        annotations = self.annotate_samples(self.data, markers, return_nonzero_annotations=False)
        self.assertEqual(len(annotations), len(self.data))
        self.assertEqual(len(annotations[0]), 3)  # two types in the data
        self.assertGreater(np.nansum(annotations.X), 0)
        self.assertLessEqual(np.nanmax(annotations.X), 1)
        self.assertGreaterEqual(np.nanmin(annotations.X), 0)

    @dense_sparse
    def test_sf(self, array):
        """
        Test annotations with hypergeom.sf
        """
        self.data.X = array(self.data.X)
        annotations = self.annotate_samples(self.data, self.markers, p_value_fun=PFUN_HYPERGEOMETRIC)

        self.assertEqual(type(annotations), Table)
        self.assertEqual(len(annotations), len(self.data))
        self.assertEqual(len(annotations[0]), 2)  # two types in the data
        self.assertGreater(np.nansum(annotations.X), 0)
        self.assertLessEqual(np.nanmax(annotations.X), 1)
        self.assertGreaterEqual(np.nanmin(annotations.X), 0)

    @dense_sparse
    def test_two_example(self, array):
        self.data = self.data[:2]

        self.data.X = array(self.data.X)
        annotations = self.annotate_samples(self.data, self.markers)

        self.assertEqual(type(annotations), Table)
        self.assertEqual(len(annotations), len(self.data))

    @dense_sparse
    def test_markers_without_entrez_id(self, array):
        self.markers[1, "Entrez ID"] = "?"
        self.data.X = array(self.data.X)
        annotations = self.annotate_samples(self.data, self.markers, return_nonzero_annotations=False)

        self.assertEqual(type(annotations), Table)
        self.assertEqual(len(annotations), len(self.data))
        self.assertEqual(len(annotations[0]), 2)  # two types in the data
        self.assertGreater(np.nansum(annotations.X), 0)
        self.assertLessEqual(np.nanmax(annotations.X), 1)
        self.assertGreaterEqual(np.nanmin(annotations.X), 0)

    @dense_sparse
    def test_select_attributes(self, array):
        self.data.X = array(self.data.X)
        z = self.annotator.mann_whitney_test(self.data)

        self.assertEqual(z.X.shape, self.data.X.shape)
        self.assertGreaterEqual(z.X[0, 0], 1)
        self.assertGreaterEqual(z.X[0, 1], 1)
        self.assertGreaterEqual(z.X[0, 3], 1)

    @dense_sparse
    def test_assign_annotations(self, array):
        z = np.array(
            [
                [1.1, 1.1, 1.1, 1.1, 0, 0, 0, 0],
                [1.1, 1.1, 0, 0, 1.1, 0, 0, 0],
                [1.1, 0, 0, 0, 1.1, 1.1, 0, 0],
                [0, 0, 0, 0, 1.1, 1.1, 1.1, 1.1],
            ]
        )
        z_table = Table(self.domain, z)
        attrs = [
            {"111", "112", "113", "114"},
            {"111", "112", "211"},
            {"211", "212", "111"},
            {"211", "212", "213", "214"},
        ]
        exp_ann = np.array([[1, 0], [1 / 2, 1 / 4], [1 / 4, 1 / 2], [0, 1]])
        self.data.X = array(self.data.X)
        annotations, fdrs = self.annotator.assign_annotations(z_table, self.markers, self.data[:4])

        self.assertEqual(len(attrs), len(annotations))
        self.assertEqual(len(attrs), len(fdrs))
        self.assertEqual(2, annotations.X.shape[1])  # only two types in markers
        self.assertEqual(2, fdrs.X.shape[1])
        np.testing.assert_array_almost_equal(exp_ann, annotations)

        exp_fdrs_smaller = np.array([[0.05, 2], [2, 2], [2, 2], [2, 0.05]])

        np.testing.assert_array_less(fdrs, exp_fdrs_smaller)

    @dense_sparse
    def test_scoring(self, array):
        # scoring SCORING_EXP_RATIO
        # last two cases with explicit zero are skipped due to pandas issue
        # https://github.com/pandas-dev/pandas/issues/28992
        if isinstance(array, functools.partial):
            return
        self.data.X = array(self.data.X)
        annotations = self.annotate_samples(self.data, self.markers, scoring=SCORING_EXP_RATIO)

        self.assertEqual(type(annotations), Table)
        self.assertEqual(len(annotations), len(self.data))
        self.assertEqual(len(annotations[0]), 2)  # two types in the data
        self.assertGreater(np.nansum(annotations.X), 0)
        self.assertLessEqual(np.nanmax(annotations.X), 1)
        self.assertGreaterEqual(np.nanmin(annotations.X), 0)

        # scoring SCORING_MARKERS_SUM
        annotations = self.annotate_samples(self.data, self.markers, scoring=SCORING_MARKERS_SUM)

        self.assertEqual(type(annotations), Table)
        self.assertEqual(len(annotations), len(self.data))
        self.assertEqual(len(annotations[0]), 2)  # two types in the data

        # based on provided data it should match
        # the third row is skipped, since it is special
        self.assertEqual(annotations[0, 0].value, self.data.X[0].sum())
        self.assertEqual(annotations[5, 1].value, self.data.X[5].sum())

        # scoring SCORING_LOG_FDR
        annotations = self.annotate_samples(self.data, self.markers, scoring=SCORING_LOG_FDR)

        self.assertEqual(type(annotations), Table)
        self.assertEqual(len(annotations), len(self.data))
        self.assertEqual(len(annotations[0]), 2)  # two types in the data

        # scoring SCORING_LOG_PVALUE
        annotations = self.annotate_samples(self.data, self.markers, scoring=SCORING_LOG_PVALUE)

        self.assertEqual(type(annotations), Table)
        self.assertEqual(len(annotations), len(self.data))
        self.assertEqual(len(annotations[0]), 2)  # two types in the data

    @dense_sparse
    def test_entrez_id_not_string(self, array):
        """
        It seems that some datasets (e.g. AML dataset) have Entrez ID as int
        although they should be strings. Here we add the test for those cases.
        """
        # change Entrez IDs to int
        for i, att in enumerate(self.data.domain.attributes):
            att.attributes["Entrez ID"] = int(att.attributes["Entrez ID"])

        self.data.X = array(self.data.X)
        annotations = self.annotate_samples(self.data, self.markers)

        self.assertEqual(type(annotations), Table)
        self.assertEqual(len(annotations), len(self.data))
        self.assertEqual(len(annotations[0]), 2)  # two types in the data
        self.assertGreater(np.nansum(annotations.X), 0)
        self.assertLessEqual(np.nanmax(annotations.X), 1)
        self.assertGreaterEqual(np.nanmin(annotations.X), 0)

    @dense_sparse
    def test_celltypes_as_feature(self, array):
        """
        Test case when cell type or genes appears as a feature
        """
        ct_values = self.markers.get_column_view("Cell Type")[0].tolist()
        entrezid = self.markers.get_column_view("Entrez ID")[0].tolist()
        ct_variable = DiscreteVariable("Cell Type", values=list(set(ct_values)))
        entrez_variable = DiscreteVariable("Entrez ID", values=list(set(entrezid)))

        markers = Table(
            Domain([ct_variable, entrez_variable]),
            list(
                zip(
                    [ct_variable.values.index(x) for x in ct_values],
                    [entrez_variable.values.index(x) for x in entrezid],
                )
            ),
        )
        self.data.X = array(self.data.X)
        annotations = self.annotate_samples(self.data, markers)

        self.assertEqual(type(annotations), Table)
        self.assertEqual(len(annotations), len(self.data))
        self.assertEqual(len(annotations[0]), 2)  # two types in the data
        self.assertGreater(np.nansum(annotations.X), 0)
        self.assertLessEqual(np.nanmax(annotations.X), 1)
        self.assertGreaterEqual(np.nanmin(annotations.X), 0)

    @dense_sparse
    def test_celltypes_as_feature_and_meta(self, array):
        """
        Test case when cell type is a feature and Entrez id is in metas
        """
        ct_values = self.markers.get_column_view("Cell Type")[0].tolist()
        entrezid = self.markers.get_column_view("Entrez ID")[0].tolist()
        ct_variable = DiscreteVariable("Cell Type", values=list(set(ct_values)))
        entrez_variable = DiscreteVariable("Entrez ID", values=list(set(entrezid)))

        markers = Table(
            Domain([ct_variable], metas=[entrez_variable]),
            list(
                zip(
                    [ct_variable.values.index(x) for x in ct_values],
                    [entrez_variable.values.index(x) for x in entrezid],
                )
            ),
        )
        self.data.X = array(self.data.X)
        annotations = self.annotate_samples(self.data, markers)

        self.assertEqual(type(annotations), Table)
        self.assertEqual(len(annotations), len(self.data))
        self.assertEqual(len(annotations[0]), 2)  # two types in the data
        self.assertGreater(np.nansum(annotations.X), 0)
        self.assertLessEqual(np.nanmax(annotations.X), 1)
        self.assertGreaterEqual(np.nanmin(annotations.X), 0)
