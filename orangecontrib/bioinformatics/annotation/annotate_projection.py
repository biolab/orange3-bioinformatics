"""
This module cluster the projection of data (usually it is 2D projection) with
one of the standard algorithms and attach a certain number of labels per
cluster.

Example:

>>> from Orange.projection import TSNE
>>> from Orange.data import Table
>>> from orangecontrib.bioinformatics.utils import serverfiles
>>> from orangecontrib.bioinformatics.annotation.annotate_projection import \
...     annotate_projection
>>> from orangecontrib.bioinformatics.annotation.annotate_samples import \
...     AnnotateSamples
>>>
>>> # load data
>>> data = Table("https://datasets.orange.biolab.si/sc/aml-1k.tab.gz")
>>> marker_p = serverfiles.localpath_download(
...     'marker_genes','panglao_gene_markers.tab')
>>> markers = Table(marker_p)
>>>
>>> # annotate data with labels
>>> annotator = AnnotateSamples()
>>> annotations = annotator.annotate_samples(data, markers)
>>>
>>> # project data in 2D
>>> tsne = TSNE(n_components=2)
>>> tsne_model = tsne(data)
>>> embedding = tsne_model(data)
>>>
>>> # get clusters and annotations for clusters
>>> clusters, clusters_meta, eps = annotate_projection(annotations, embedding,
...     clustering_algorithm=DBSCAN, eps=1.2)

In case when user uses a DBSCAN algorithm and do not provide eps to the
`annotate_projection` function it is computed automatically with knn method.

>>> clusters, clusters_meta, eps = annotate_projection(annotations, embedding,
...     clustering_algorithm=DBSCAN)

"""
from math import sqrt
from bisect import bisect_left, bisect_right
from collections import Counter

import numpy as np
import pyclipper
from scipy.spatial import Delaunay, distance

from Orange.data import Table, Domain, DiscreteVariable
from Orange.clustering import DBSCAN


def cluster_data(coordinates, clustering_algorithm=DBSCAN, **kwargs):
    """
    This function receives data and cluster them.

    Parameters
    ----------
    coordinates : Orange.data.Table
        Visualisation coordinates - embeddings
    clustering_algorithm : callable
        Algorithm used for clustering.

    Returns
    -------
    Orange.data.Table
        List of cluster indices.
    """
    learner = clustering_algorithm(**kwargs)
    clusters = learner(coordinates)
    if not isinstance(clusters, np.ndarray):  # old clustering method
        clusters = clusters(coordinates)
        clusters = np.array(list(map(int, map(clusters.domain.attributes[0].repr_val, clusters.X[:, 0])))).flatten()

    # sort classes in descending order base on number of cases in the cluster
    sorted_clust_idx = [v for v, _ in Counter(clusters).most_common() if v != -1]

    # re-indexed array
    new_clustering = np.empty(len(clusters))
    new_clustering[:] = np.nan  # nan for not clustered
    # reindex based on descending cluster size
    for i, v in enumerate(sorted_clust_idx):
        new_clustering[clusters == v] = i

    # create the table
    new_domain = Domain(
        [DiscreteVariable("Clusters", values=["C{}".format(i) for i in range(1, len(sorted_clust_idx) + 1)])]
    )
    return Table(new_domain, new_clustering.reshape((-1, 1)))


def assign_labels(clusters, annotations, labels_per_cluster):
    """
    This function assigns a certain number of labels per cluster. Each cluster
    gets `labels_per_cluster` number of most common labels in cluster assigned.

    Parameters
    ----------
    clusters : Orange.data.Table
        Cluster indices for each item.
    annotations : Orange.data.Table
        Table with annotations and their probabilities.
    labels_per_cluster : int
        Number of labels that need to be assigned to each cluster.

    Returns
    -------
    dict
        Dictionary with cluster index as a key and list of annotations as a
        value. Each list include tuples with the annotation name and their
        proportion in the cluster.
    Orange.data.Table
        The array with the annotation assigned to the item.
    """
    clusters_unique = set(clusters.domain[0].values)

    if len(annotations.domain) == 0:
        return {}, Table(Domain([DiscreteVariable("Annotation", values=[])]), np.ones((len(clusters), 1)) * np.nan)

    labels = np.array(list(map(str, annotations.domain.attributes)))

    # remove rows with all nans
    nan_mask = np.isnan(annotations.X).all(axis=1)
    ann_not_nan = annotations.X[~nan_mask]

    # find indices and labels
    annotation_best_idx = np.nanargmax(ann_not_nan, axis=1)
    annotation_best = labels[annotation_best_idx]

    # join back together
    items_annotations = np.empty(annotations.X.shape[0], dtype=labels.dtype)
    items_annotations[~nan_mask] = annotation_best

    annotations_clusters = {}
    for cl in clusters_unique:
        mask = np.array(list(map(clusters.domain.attributes[0].repr_val, clusters.X[:, 0]))).flatten() == cl
        labels_cl = items_annotations[mask]
        # remove nans from labels
        labels_cl_filtered = labels_cl[~(labels_cl == "")]

        counts = Counter(labels_cl_filtered)
        common_labels = counts.most_common(labels_per_cluster)

        if len(common_labels) > 0:
            annotations_clusters[cl] = [(l, c / len(labels_cl)) for l, c in common_labels]

    # pack item annotations to Table
    nan_mask = items_annotations == ""
    values, indices = np.unique(items_annotations[~nan_mask], return_inverse=True)
    corrected_idx = np.ones(items_annotations.shape) * np.nan
    corrected_idx[~nan_mask] = indices
    domain = Domain([DiscreteVariable("Annotation", values=values)])
    item_annotations = Table(domain, corrected_idx.reshape((-1, 1)))

    return annotations_clusters, item_annotations


def labels_locations(coordinates, clusters):
    """
    Function computes the location of the label for each cluster.
    The location is compute as a center point.

    Parameters
    ----------
    coordinates : Orange.data.Table
        Visualisation coordinates - embeddings
    clusters : Orange.data.Table
        Cluster indices for each item.

    Returns
    -------
    dict
        The coordinates for locating the label. Dictionary with cluster index
        as a key and tuple (x, y) as a value.
    """
    clusters_unique = set(clusters.domain[0].values) - {"-1"}  # -1 is not clustered
    locations = {}
    for cl in clusters_unique:
        mask = np.array(list(map(clusters.domain.attributes[0].repr_val, clusters.X[:, 0]))).flatten() == cl
        cl_coordinates = coordinates.X[mask, :]
        x, y = 1 / 2 * (np.min(cl_coordinates, axis=0) + np.max(cl_coordinates, axis=0))
        locations[cl] = (x, y)
    return locations


def get_epsilon(coordinates, k=10, skip=0.1):
    """
    The function computes the epsilon parameter for DBSCAN through method
    proposed in the paper.

    Parameters
    ----------
    coordinates : Orange.data.Table
        Visualisation coordinates - embeddings
    k : int
        Number kth observed neighbour
    skip : float
        Percentage of skipped neighborus.

    Returns
    -------
    float
        Epsilon parameter for DBSCAN
    """
    x = coordinates.X
    if len(x) > 1000:  # subsampling is required
        i = len(x) // 1000
        x = x[::i]

    d = distance.squareform(distance.pdist(x))
    k = min(k + 1, len(coordinates) - 1)
    kth_point = np.argpartition(d, k, axis=1)[:, k]
    # k+1 since first one is item itself
    kth_dist = np.sort(d[np.arange(0, len(kth_point)), kth_point])

    # currently mark proportion equal to skip as a noise
    return kth_dist[-int(np.round(len(kth_dist) * skip))]


def _angle(v1, v2):
    """
    Compute clockwise angles between v1 and v2. Both vectors are
    given as a tuple with two points ([x1, y1], [x2, y2]) such that
    [x2, y2] of v1 == [x1, y1] of v2.
    The angle is given in degrees between 0 and 2 * pi
    """
    v1, v2 = np.array(v1), np.array(v2)
    x1, y1 = v1[0] - v1[1]
    x2, y2 = v2[1] - v2[0]
    dot = np.sum(x1 * x2 + y1 * y2)  # dot product
    det = x1 * y2 - y1 * x1  # determinant
    angle = np.arctan2(det, dot)  # atan2(y, x) or atan2(sin, cos)
    return -angle if angle <= 0 else np.pi + np.pi - angle


def _find_hull(edges_list, points_list, starting_edge):
    """
    This function return a single hull which starts and ends in
    the starting_edge.

    Parameters
    edges_list : list
        List of edges. Each edge is presented as a tuple of two indices which
        tell the starting and ending note. The index correspond to location
        of point in points_list. This list must be sorted in the ascending
        order.
    points_list : list
        List of points location. Each point has x and y location.
    starting_edge : int
        The index of the list where hull starts.

    Returns
    -------
    np.ndarray
        The array with the hull/polygon points
    list
        List of booleans that indicates whether each point was used or not in
        the resulting polygon.
    """
    firsts = [x[0] for x in edges_list]  # first elements for bisection
    used_edges = [False] * len(edges_list)
    used_edges[starting_edge] = True

    # remember start and next point
    # start point is required for stop condition
    start, next_point = edges_list[starting_edge]

    # remember current polygon around points
    poly = [points_list[start], points_list[next_point]]

    # we count number of steps to stop iteration in case it is locked
    # in some dead cycle. It can be a result of some unexpected cases.
    count = 0
    while start != next_point and count < len(edges_list):
        # find the index where the first value equal to next_point
        # appear
        ind_left = bisect_left(firsts, next_point)
        # find the index next to the last value
        ind_right = bisect_right(firsts, next_point)

        # check if there are more edges available from the same point
        if ind_right - ind_left > 1:
            # select the most distant one in clockwise direction
            # it is probably the point on the outer hull
            # with this we prevent a hull to discover cycles inside
            # a polygon
            ang = -1
            for i in range(ind_left, ind_right):
                cur_ang = _angle((poly[-2], poly[-1]), (poly[-1], points_list[edges_list[i][1]]))
                if cur_ang > ang:
                    ang = cur_ang
                    ind_left = i
        # save a next point of the polgon
        used_edges[ind_left] = True
        next_point = edges_list[ind_left][1]
        poly.append(points_list[next_point])
        count += 1
    return np.array(poly), used_edges


def edges_to_polygon(edges_list, points_list):
    """
    This function connect edges in polygons. It computes all possible hulls -
    yes some clusters have more of them when they have a hole in the middle.
    It then selects one that is outer hull.

    Parameters
    ----------
    edges_list : list
        List of edges. Each edge is presented as a tuple of two indices which
        tell the starting and ending note. The index correspond to location
        of point in points_list
    points_list : list
        List of points location. Each point has x and y location.

    Returns
    -------
    np.ndarray
        The array with the hull/polygon points.
    """
    # sort based on first element of tuple to enable bisection search
    edges_list = sorted(edges_list, key=lambda x: x[0])
    # need to use all edges
    used = [False] * len(edges_list)

    # it is possible that we will find more separate hulls -
    # it happen in cases when a polygon has inner cycles
    polygons = []
    while not all(used):
        i = used.index(False)
        poly, new_used = _find_hull(edges_list, points_list, i)
        polygons.append(poly)
        used = [u1 or u2 for u1, u2 in zip(used, new_used)]

    # select polygon that is outside - the widest and the highest is the
    # most outside
    height_width = [np.sum(p.max(axis=0) - p.min(axis=0)) for p in polygons]
    i = height_width.index(max(height_width))

    return polygons[i]


def compute_concave_hulls(coordinates, clusters, epsilon):
    """
    Function computes the points of the concave hull around points.

    Parameters
    ----------
    coordinates : Orange.data.Table
        Visualisation coordinates - embeddings
    clusters : Orange.data.Table
       Cluster indices for each item.
    epsilon : float
        Epsilon used by DBSCAN to cluster the data

    Returns
    -------
    dict
       The points of the concave hull. Dictionary with cluster index
       as a key and np.ndaray of points as a value -
       [[x1, y1], [x2, y2], [x3, y3], ...]
    """

    def get_shape(pts, eps):
        """
        Compute the shape (concave hull) of a set of a cluster.
        """

        def add_edge(edges_list, i, j):
            """
            Add a line between the i-th and j-th points,
            if not in the list already. Remove the lines that are not outer
            edges - when (j, i) already in list means that two triangles
            has same edge what means it is not an outer edge
            """
            if (j, i) in edges_list:
                # if both neighboring triangles are in shape, it's not a
                # boundary edge
                edges_list.remove((j, i))
                return
            edges_list.add((i, j))

        if len(pts) < 4:
            # When you have a triangle, there is no sense in computing the hull
            rng = list(range(3))
            return edges_to_polygon(zip(rng, rng[1:] + rng[:1]), pts)

        tri = Delaunay(pts)
        edges = set()
        # loop over triangles:
        # ia, ib, ic = indices of corner points of the triangle
        for ia, ib, ic in tri.vertices:
            pa = pts[ia]
            pb = pts[ib]
            pc = pts[ic]

            # Lengths of sides of triangle
            a = sqrt((pa[0] - pb[0]) ** 2 + (pa[1] - pb[1]) ** 2)
            b = sqrt((pb[0] - pc[0]) ** 2 + (pb[1] - pc[1]) ** 2)
            c = sqrt((pc[0] - pa[0]) ** 2 + (pc[1] - pa[1]) ** 2)

            # filter - longest edge of triangle smaller than epsilon
            if max(a, b, c) <= eps:
                add_edge(edges, ia, ib)
                add_edge(edges, ib, ic)
                add_edge(edges, ic, ia)

        polygon = edges_to_polygon(edges, pts)
        return polygon

    hulls = {}
    clusters_array = np.array(list(map(clusters.domain.attributes[0].repr_val, clusters.X[:, 0])))
    for cl in set(clusters_array) - {"None", "?"}:
        points = coordinates.X[clusters_array == cl]

        # subsample when more than 1000 points
        # it keeps time finding hull under 0.3 s on my computer
        if points.shape[0] > 1000:
            points = points[np.random.randint(points.shape[0], size=1000), :]

        # compute the concave hul
        hull = get_shape(points, eps=epsilon * 2)

        # epsilon seems to be good parameter for lines to be smooth enough
        # selecting epsilon for the distance
        # shows approximately what is DBSCAN neighbourhood
        # those line buffer the line for 3 * epsilon and move it back for
        # 2 * epsilon. The it will make a hull more smooth.

        # pyclipper work with integer so points need to be scaled first to keep
        # precision. It is undo with scal_from_clipper
        scaling_factor = 1000
        scaled_hull = pyclipper.scale_to_clipper(hull, scaling_factor)

        # buffer the hull fro epsilon * 3
        pco1 = pyclipper.PyclipperOffset()
        pco1.AddPath(scaled_hull, pyclipper.JT_ROUND, pyclipper.ET_CLOSEDPOLYGON)
        im_solution = pco1.Execute(epsilon * scaling_factor * 3)
        # buffer the hull fro epsilon * -2
        pco2 = pyclipper.PyclipperOffset()
        pco2.AddPath(im_solution[0], pyclipper.JT_ROUND, pyclipper.ET_CLOSEDPOLYGON)
        solution = pyclipper.scale_from_clipper(pco2.Execute(epsilon * scaling_factor * (-2)), scaling_factor)

        hulls[cl] = np.array(solution).reshape(-1, 2)
    return hulls


def _filter_clusters(clusters, clusters_meta):
    """
    Function removes clusters that does not has any labels
    """
    clust_map = {c: i for i, c in enumerate(c for c in clusters.domain["Clusters"].values if c in clusters_meta)}

    # change cluster indices
    clust_idx = np.zeros(clusters.X.shape) * np.nan
    for i, cl in enumerate(map(clusters.domain.attributes[0].repr_val, clusters.X[:, 0])):
        clust_idx[i, 0] = clust_map.get(cl, np.nan)
    new_clusters = Table(
        Domain([DiscreteVariable("Clusters", values=clusters.domain["Clusters"].values[: len(clust_map)])]), clust_idx
    )
    # change cluster names in metas
    new_clusters_meta = {"C{}".format(clust_map[cl] + 1): v for cl, v in clusters_meta.items()}
    return new_clusters, new_clusters_meta


def annotate_projection(annotations, coordinates, clustering_algorithm=DBSCAN, labels_per_cluster=3, **kwargs):
    """
    Function cluster the data based on coordinates, and assigns a certain number
    of labels per cluster. Each cluster gets `labels_per_cluster` number of most
    common labels in cluster assigned.

    Parameters
    ----------
    annotations : Orange.data.Table
        Table with annotations and their probabilities.
    coordinates : Orange.data.Table
        Visualisation coordinates - embeddings
    clustering_algorithm : callable, optional (default = DBSCAN)
        Algorithm used in clustering.
    labels_per_cluster : int, optional (default = 3)
        Number of labels that need to be assigned to each cluster.

    Returns
    -------
    Orange.data.Table
        List of cluster indices.
    dict
        Dictionary with cluster index as a key and list of annotations as a
        value. Each list include tuples with the annotation name and their
        proportion in the cluster.
    dict
        The coordinates for locating the label. Dictionary with cluster index
        as a key and tuple (x, y) as a value.
    """
    assert len(annotations) == len(coordinates), "Number of coordinates does not match to number of annotations"
    # sklearn clustering want to have one example
    assert len(coordinates) > 0, "At least one data point need to be provided"
    assert len(coordinates.domain) > 0, "Coordinates need to have at least one attribute"

    eps = kwargs.get("eps", get_epsilon(coordinates))
    if clustering_algorithm == DBSCAN:
        kwargs["eps"] = eps

    # get clusters
    clusters = cluster_data(coordinates, clustering_algorithm, **kwargs)

    # assign top n labels to group
    annotations_cl, item_annotations = assign_labels(clusters, annotations, labels_per_cluster)

    labels_loc = labels_locations(coordinates, clusters)

    concave_hull = compute_concave_hulls(coordinates, clusters, eps)

    # crate the dictionary with annotations, labels locations, and hulls for
    # each cluster
    clusters_meta = {}
    # note clusters without labels are not added to the dictionary
    for cl in annotations_cl.keys():
        clusters_meta[cl] = (annotations_cl[cl], labels_loc[cl], concave_hull[cl])

    # we do not want to color clusters that does not have labels so we
    # remove clusters without labels from the table
    clusters, clusters_meta = _filter_clusters(clusters, clusters_meta)

    # add the labels to the cluster table
    clusters_ann = Table(
        Domain(clusters.domain.attributes + item_annotations.domain.attributes),
        np.concatenate((clusters.X, item_annotations.X), axis=1),
    )

    return clusters_ann, clusters_meta, eps


def cluster_additional_points(coordinates, hulls, cluster_attribute=None):
    """
    This function receives additional points and assign them current existing
    clusters based on current concave hull.

    Parameters
    ----------
    coordinates : Orange.data.Table
        Visualisation coordinates - embeddings
    hulls : dict
        Concave hull for each cluster
    cluster_attribute : Orange.data.DiscreteVariable (optional)
        A variable for clusters. If cluster_attribute is provided it will be
        used in the creation of the resulting Table.

    Returns
    -------
    Orange.data.Table
        Cluster label for each point
    """

    def point_in_polygon_test(test_point, polygon_points):
        """
        This function uses the horizontal ray casting to find out if the point
        is in the hull/polygon. For each point, it tests how many times the
        horizontal ray from test_point to infinity crosses the polygon edge. If
        it happens odd many times the point is in the polygon.
        https://stackoverflow.com/a/2922778/3551700
        """
        test_x = test_point[0]
        test_y = test_point[1]
        # flipping bool from True to False is similar to counting odd numbers
        # of intersections. If it will be True at the end odd number of
        # intersections happened
        is_inside = False

        for (x1, y1), (x2, y2) in zip(
            polygon_points, np.concatenate((polygon_points[1:], polygon_points[:1]), axis=0)
        ):
            # ray crosses the edge if test_y between both y from an edge
            # and if intersection on the right of the test_x
            if (y1 > test_y) != (y2 > test_y):
                # compute the intersection between the horizontal ray and
                # polygon edge
                intersection_x = (x2 - x1) * (test_y - y1) / (y2 - y1) + x1
                if test_x < intersection_x:
                    is_inside = not is_inside
        return is_inside

    clusters = [None] * len(coordinates)
    for cluster, hull in hulls.items():
        for i, c in enumerate(coordinates.X):
            if point_in_polygon_test(c, hull):
                clusters[i] = cluster

    if cluster_attribute is not None:
        assert all(
            i in cluster_attribute.values for i in set(clusters) - {None}
        ), "cluster_attribute does not have all required values."
    # create the table
    new_domain = Domain(
        [
            DiscreteVariable("Clusters", values=sorted(list(hulls.keys())))
            if cluster_attribute is None
            else cluster_attribute
        ]
    )
    return Table(new_domain, np.array(list(map(new_domain[0].to_val, clusters))).reshape(-1, 1))


if __name__ == "__main__":
    # run hull creation at Iris data
    data = Table("iris")[:, 2:4]
    clustered_data = Table(
        Domain([DiscreteVariable("cl", values=["1", "2", "3"])]), [[0]] * 50 + [[1]] * 50 + [[2]] * 50
    )
    compute_concave_hulls(data, clustered_data, epsilon=0.5)
