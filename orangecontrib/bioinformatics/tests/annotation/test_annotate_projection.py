import unittest

import numpy as np

from Orange.data import Table, Domain, DiscreteVariable, ContinuousVariable
from Orange.clustering import DBSCAN, KMeans

from orangecontrib.bioinformatics.annotation.annotate_projection import (
    _angle,
    get_epsilon,
    assign_labels,
    _filter_clusters,
    labels_locations,
    annotate_projection,
    compute_concave_hulls,
    cluster_additional_points,
)


class TestAnnotateProjection(unittest.TestCase):
    def setUp(self):
        domain = Domain([ContinuousVariable("x"), ContinuousVariable("y")])
        data = np.array(
            [
                [1.0, 1.0],
                [1.0, 1.5],
                [1.5, 1.0],
                [1.2, 1.2],
                [1.7, 1.8],
                [10.0, 10.0],
                [10.0, 10.5],
                [10.5, 11.0],
                [10.5, 10.7],
                [10.5, 10.6],
                [10.5, 10.7],
                [10.5, 10.7],
            ]
        )
        self.data = Table.from_numpy(domain, data)

        domain_ann = Domain([ContinuousVariable("a"), ContinuousVariable("b"), ContinuousVariable("c")])
        data_ann = np.array(
            [
                [0.9, 1, 0.4],
                [0.9, 0.4, 0.2],
                [0.9, 0.4, 0.2],
                [0.9, 1, 0.2],
                [0.9, 0.4, 0.2],
                [0.2, 0.4, 0.7],
                [0.1, 0.4, 0.6],
                [0.2, 0.9, 0.2],
                [0.1, 0.4, 0.6],
                [0.1, 0.8, 0.6],
                [0.1, 0.8, 0.6],
                [0.1, 0.8, 0.6],
            ]
        )
        self.annotations = Table.from_numpy(domain_ann, data_ann)

    def test_annotate_projection(self):
        clusters, clusters_meta, eps = annotate_projection(
            self.annotations, self.data, clustering_algorithm=DBSCAN, labels_per_cluster=2, eps=3.0
        )

        self.assertListEqual(
            list(map(clusters.domain.attributes[0].repr_val, clusters.X[:, 0])), ["C2"] * 5 + ["C1"] * 7
        )

        self.assertEqual(len(clusters_meta), 2)
        self.assertEqual(len(clusters_meta["C1"]), 3)
        self.assertEqual(len(clusters_meta["C2"]), 3)
        self.assertEqual(clusters_meta["C1"][0][0][0], 'b')
        self.assertAlmostEqual(clusters_meta["C1"][0][0][1], 4 / 7, 5)
        self.assertEqual(clusters_meta["C1"][0][1][0], 'c')
        self.assertAlmostEqual(clusters_meta["C1"][0][1][1], 3 / 7, 5)

        self.assertEqual(2, len(clusters_meta["C1"][1]))
        self.assertEqual(2, clusters_meta["C1"][2].shape[1])
        self.assertEqual(2, clusters_meta["C1"][2].shape[1])
        self.assertEqual(float, type(eps))

    def test_example_not_clustered(self):
        self.data[-1] = [23, 23]
        clusters, clusters_meta, eps = annotate_projection(
            self.annotations, self.data, clustering_algorithm=DBSCAN, eps=2, min_samples=3
        )

        self.assertListEqual(
            list(map(clusters.domain.attributes[0].repr_val, clusters.X[:, 0])), ["C2"] * 5 + ["C1"] * 6 + ["?"]
        )

        self.assertEqual(len(clusters_meta), 2)
        self.assertEqual(len(clusters_meta["C1"]), 3)
        self.assertEqual(len(clusters_meta["C2"]), 3)
        self.assertEqual(clusters_meta["C1"][0][0][0], 'c')
        self.assertAlmostEqual(clusters_meta["C1"][0][0][1], 0.5, 5)
        self.assertEqual(clusters_meta["C1"][0][1][0], 'b')
        self.assertAlmostEqual(clusters_meta["C1"][0][1][1], 0.5, 5)

    def test_other_clustering(self):
        clusters, clusters_meta, locs = annotate_projection(
            self.annotations, self.data, clustering_algorithm=KMeans, n_clusters=2, labels_per_cluster=2
        )

        self.assertEqual(2, len(set(clusters.X[:, 0].flatten())))
        self.assertEqual(3, len(clusters_meta["C1"]))
        self.assertEqual(3, len(clusters_meta["C2"]))
        self.assertEqual(2, len(clusters_meta["C1"][0]))
        self.assertEqual(2, len(clusters_meta["C2"][0]))

    def test_labels_location(self):
        clusters = Table.from_list(Domain([DiscreteVariable("cl", values=["1", "2"])]), [[0]] * 6 + [[1]] * 6)
        locs = labels_locations(self.data, clusters)

        self.assertEqual(dict, type(locs))
        self.assertEqual(tuple, type(locs["1"]))
        self.assertEqual(tuple, type(locs["2"]))

    def test_get_epsilon(self):
        data = Table.from_file("iris")
        eps = get_epsilon(data)

        self.assertGreaterEqual(eps, 0.9)

    def test_compute_concave_hulls(self):
        data = Table.from_file("iris")[:, 2:4]
        clusters = Table.from_list(
            Domain([DiscreteVariable("cl", values=["1", "2", "3"])]), [[0]] * 50 + [[1]] * 50 + [[2]] * 50
        )
        hulls = compute_concave_hulls(data, clusters, epsilon=0.5)

        self.assertEqual(3, len(hulls))
        self.assertEqual(2, hulls["1"].shape[1])  # hull have x and y
        self.assertEqual(2, hulls["2"].shape[1])  # hull have x and y
        self.assertEqual(2, hulls["3"].shape[1])  # hull have x and y

    def test_compute_concave_hulls_subsampling(self):
        """
        When more than 1000 points passed they are sub-sampled in order to
        compute a concave hull
        """
        iris = Table.from_file("iris")
        data = Table.from_numpy(  # here we pass 1500 points
            Domain(iris.domain.attributes[2:4]), np.repeat(iris.X[:, 2:4], 10, axis=0)
        )
        clusters = Table.from_list(
            Domain([DiscreteVariable("cl", values=["1", "2", "3"])]),
            [[0]] * 50 * 10 + [[1]] * 50 * 10 + [[2]] * 50 * 10,
        )
        hulls = compute_concave_hulls(data, clusters, epsilon=0.5)

        self.assertEqual(3, len(hulls))
        self.assertEqual(2, hulls["1"].shape[1])  # hull have x and y
        self.assertEqual(2, hulls["2"].shape[1])  # hull have x and y
        self.assertEqual(2, hulls["3"].shape[1])  # hull have x and y

    def test_compute_concave_hulls_triangle(self):
        """
        Concave hull must also work for tree points - it is a special case
        """
        data = Table.from_numpy(
            Domain([ContinuousVariable("x"), ContinuousVariable("y")]), np.array([[1, 1], [1, 2], [2, 1]])
        )
        clusters = Table.from_list(Domain([DiscreteVariable("cl", values=["1"])]), [[0]] * 3)
        hulls = compute_concave_hulls(data, clusters, epsilon=0.5)

        self.assertEqual(1, len(hulls))
        self.assertEqual(2, hulls["1"].shape[1])  # hull have x and y

    def test_empty_annotations(self):
        ann = Table.from_numpy(Domain([]), np.empty((len(self.data), 0)))
        clusters, clusters_meta, eps = annotate_projection(ann, self.data)

        self.assertGreaterEqual(len(clusters), 0)
        self.assertEqual(len(clusters_meta), 0)

    def test_one_label(self):
        """
        Test whether having only one label works fine, in this case one cluster
        will not have a label assigned - this must work fine as well.
        """
        domain_ann = Domain([ContinuousVariable("a")])
        data_ann = np.array([[0.9], [0.9], [0.9], [0.9], [0.9], [0], [0], [0], [0], [0], [0], [0]])
        ann = Table.from_numpy(domain_ann, data_ann)
        clusters, clusters_meta, eps = annotate_projection(ann, self.data, eps=2)
        self.assertGreater(len(clusters), 0)
        self.assertGreater(len(clusters_meta), 0)
        self.assertGreater(len(clusters_meta["C1"][0]), 0)
        self.assertGreater(len(clusters_meta["C1"][1]), 0)
        self.assertGreater(len(clusters_meta["C1"][2]), 0)

    def test_socres_with_nan(self):
        self.annotations.X[0, 0] = np.nan
        self.annotations.X[1, 1] = np.nan
        self.annotations.X[1, 2] = np.nan

        self.annotations.X[2, 0] = np.nan
        self.annotations.X[2, 1] = np.nan
        self.annotations.X[2, 2] = np.nan

        clusters = Table.from_numpy(
            Domain([DiscreteVariable("Cluster", values=["C1", "C3"])]), np.array([[0] * 5 + [1] * 7]).reshape(-1, 1)
        )

        labels_dict, labels_items = assign_labels(clusters, self.annotations, 3)

        transformed_labels = list(map(labels_items.domain.attributes[0].repr_val, labels_items.X[:, 0]))
        self.assertListEqual(['b', 'a', '?', 'b', 'a', 'c', 'c', 'b', 'c', 'b', 'b', 'b'], transformed_labels)
        self.assertTupleEqual((0.4, 0.4), list(zip(*labels_dict["C1"]))[1])

        # check all nans
        self.annotations.X[:, :] = np.nan
        labels_dict, labels_items = assign_labels(clusters, self.annotations, 3)
        transformed_labels = list(map(labels_items.domain.attributes[0].repr_val, labels_items.X[:, 0]))
        self.assertListEqual(['?'] * len(self.annotations), transformed_labels)
        self.assertFalse("C1" in labels_dict)
        self.assertFalse("C3" in labels_dict)

    def test_cluster_additional_points_simple(self):
        """
        Test with the simple convex hull
        """
        hull = {"C1": np.array([[0.0, 0.0], [0.0, 1.0], [1.0, 1.0], [1.0, 0.0]])}
        x, y = np.meshgrid(np.arange(-1.0, 2.0, 0.2), np.arange(-1, 2.0, 0.2))
        points = Table.from_numpy(
            Domain([ContinuousVariable("x"), ContinuousVariable("y")]),
            np.concatenate((x.reshape(-1, 1), y.reshape(-1, 1)), axis=1),
        )

        clusters = cluster_additional_points(points, hull)
        expected = [0.0 if 0 < x < 1 and 0 < y < 1 else np.nan for x, y in points.X]
        np.testing.assert_array_equal(clusters.X.flatten(), expected)

    def test_cluster_additional_points_concave(self):
        """
        Test with the concave hull
        """
        hull = {"C1": np.array([[0.0, 0.0], [0.3, 0.5], [0.0, 1.0], [1.0, 1.0], [0.7, 0.5], [1.0, 0.0]])}
        p = [
            [-1.0, -1],  # point that is totally out polygon
            [0.2, 0.5],  # point outside the polygon, but close
            [0.8, 0.5],  # point outside the polygon, but close
            [0.5, 0.5],  # point in
            [0.5, 0.7],  # point in
            [0.5, 0.3],  # point in
            [0.2, 0.9],  # point in
            [0.2, 0.1],  # point in
        ]
        points = Table.from_list(Domain([ContinuousVariable("x"), ContinuousVariable("y")]), p)

        clusters = cluster_additional_points(points, hull)
        expected = [np.nan] * 3 + [0.0] * 5
        np.testing.assert_array_equal(clusters.X.flatten(), expected)

    points = [
        [-1.0, -1],  # point that is totally out  of any polygon
        [0.2, 0.5],  # point outside, but in concave part of C1
        [0.8, 0.5],  # point outside, but in concave part of C1
        [0.5, 0.5],  # point in C1
        [0.5, 0.7],  # point in C1
        [0.5, 0.3],  # point in C1
        [0.2, 0.9],  # point in C1
        [0.2, 0.1],  # point in C1
        [-0.5, 0.1],  # point in C2
        [-0.8, 0.1],  # point in C2
        [-0.5, 0.8],  # point out, but in concave part of C2
        [-0.6, 0.9],  # point out, but in concave part of C2
    ]

    hull = {
        "C1": np.array([[0.0, 0.0], [0.3, 0.5], [0.0, 1.0], [1.0, 1.0], [0.7, 0.5], [1.0, 0.0]]),
        "C2": np.array([[0.0, 0.0], [0.0, 1.0], [-0.5, 0.7], [-1.0, 1.0], [-1.0, 0.0]]),
    }

    def test_cluster_additional_points_more_clusters(self):
        """
        Test with more concave hulls
        """

        points_table = Table.from_list(Domain([ContinuousVariable("x"), ContinuousVariable("y")]), self.points)

        clusters = cluster_additional_points(points_table, self.hull)
        clusters = list(map(clusters.domain.attributes[0].repr_val, clusters.X[:, 0]))
        expected = ["?"] * 3 + ["C1"] * 5 + ["C2"] * 2 + ["?"] * 2
        np.testing.assert_array_equal(clusters, expected)

    def test_cluster_additional_points_with_existing_domain(self):
        """
        Test with the concave hull
        """
        points = Table.from_list(Domain([ContinuousVariable("x"), ContinuousVariable("y")]), self.points)
        attr = DiscreteVariable("Clusters", values=["C2", "C1"])

        clusters = cluster_additional_points(points, self.hull, cluster_attribute=attr)
        cluster_map = list(map(clusters.domain.attributes[0].repr_val, clusters.X[:, 0]))
        expected = ["?"] * 3 + ["C1"] * 5 + ["C2"] * 2 + ["?"] * 2
        np.testing.assert_array_equal(cluster_map, expected)
        self.assertEqual(attr, clusters.domain.attributes[0])
        self.assertSequenceEqual(attr.values, clusters.domain.attributes[0].values)

    def test_angle(self):
        # test cases with first vector from x, y = [1, 0] to x, y = [0, 0]
        # and second from x, y = [0, 0] to points in test_cases
        v1 = ([1, 0], [0, 0])
        v2_1 = [0, 0]
        test_cases = [[1, 0], [1, -1], [0, -1], [-1, -1], [-1, 0], [-1, 1], [0, 1], [1, 1]]
        real_angles = [np.pi * i / 4 for i in range(8)]
        computed_angles = [_angle(v1, (v2_1, x)) for x in test_cases]
        np.testing.assert_array_almost_equal(computed_angles, real_angles)

        # test cases with first vector from x, y = [2, 2] to x, y = [1, 2]
        # and second from x, y = [1, 2] to points in test_cases
        v1 = ([2, 2], [1, 2])
        v2_1 = np.array([1, 2])
        test_cases = np.array([[1, 0], [1, -1], [0, -1], [-1, -1], [-1, 0], [-1, 1], [0, 1], [1, 1]]) + v2_1
        real_angles = [np.pi * i / 4 for i in range(8)]
        computed_angles = [_angle(v1, (v2_1, x)) for x in test_cases]
        np.testing.assert_array_almost_equal(computed_angles, real_angles)

        # some more complex cases
        self.assertAlmostEqual(_angle(([-2, 1], [0, 0]), ([0, 0], [1, 2])), np.pi / 2)
        self.assertAlmostEqual(_angle(([-4, 1], [-2, 0]), ([-2, 0], [-1, 2])), np.pi / 2)
        self.assertAlmostEqual(_angle(([-4, 3], [-2, 2]), ([-2, 2], [-1, 4])), np.pi / 2)

    def test_cluster_filtering(self):
        """
        Cluster filtering was introduced in order to remove clusters that
        does not have any labels.
        """

        def assert_dict_same(d1, d2):
            self.assertTrue(len(d1) == len(d2) and sorted(d1) == sorted(d2))

        clusters = Table.from_list(Domain([DiscreteVariable("Clusters", values=["C1", "C2"])]), [[0]] * 6 + [[1]] * 5)
        metas = {"C1": "test1", "C2": "test 2"}

        new_clusters, new_metas = _filter_clusters(clusters, metas)

        # check old table unchanged
        self.assertEqual(11, len(clusters))
        self.assertSequenceEqual(["C1", "C2"], clusters.domain["Clusters"].values)
        np.testing.assert_array_equal(np.array([[0]] * 6 + [[1]] * 5), clusters.X)
        assert_dict_same(metas, {"C1": "test1", "C2": "test 2"})

        # new clusters should be same in this case
        self.assertEqual(11, len(new_clusters))
        self.assertSequenceEqual(["C1", "C2"], new_clusters.domain["Clusters"].values)
        np.testing.assert_array_equal(np.array([[0]] * 6 + [[1]] * 5), new_clusters.X)
        assert_dict_same(new_metas, {"C1": "test1", "C2": "test 2"})

        clusters = Table.from_list(Domain([DiscreteVariable("Clusters", values=["C1", "C2"])]), [[0]] * 6 + [[1]] * 5)
        metas = {"C1": "test1"}

        new_clusters, new_metas = _filter_clusters(clusters, metas)

        # check old table unchanged
        self.assertEqual(11, len(clusters))
        self.assertSequenceEqual(["C1", "C2"], clusters.domain["Clusters"].values)
        np.testing.assert_array_equal(np.array([[0]] * 6 + [[1]] * 5), clusters.X)
        assert_dict_same(metas, {"C1": "test1"})

        # new clusters should be same in this case
        self.assertEqual(11, len(new_clusters))
        self.assertSequenceEqual(["C1"], new_clusters.domain["Clusters"].values)
        np.testing.assert_array_equal(np.array([[0]] * 6 + [[np.nan]] * 5), new_clusters.X)
        assert_dict_same(new_metas, {"C1": "test1"})
