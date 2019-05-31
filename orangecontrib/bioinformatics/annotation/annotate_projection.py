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
>>> annotator = AnnotateSamples(p_value_th=0.05)
>>> annotations = annotator.annotate_samples(data, markers)
>>>
>>> # project data in 2D
>>> tsne = TSNE(n_components=2)
>>> tsne_model = tsne(data)
>>> embedding = tsne_model(data)
>>>
>>> # get clusters and annotations for clusters
>>> clusters, annotations_cl = annotate_projection(annotations, embedding,
...     clustering_algorithm=DBSCAN(eps=2))
"""


from collections import Counter

from Orange.clustering import DBSCAN
import numpy as np
from scipy.spatial import distance


def cluster_data(coordinates, clustering_algorithm, **kwargs):
    """
    This function receives data and cluster them.

    Parameters
    ----------
    coordinates : Orange.data.Table
        Data to be clustered
    clustering_algorithm : callable
        Algorithm used for clustering.

    Returns
    -------
    Orange.data.Table
        List of cluster indices.
    """
    learner = clustering_algorithm(**kwargs)
    model = learner(coordinates)
    return model(coordinates)
    # TODO: this need to be changed when clustering in orange is changed


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
    """
    labels = np.array(list(map(str, annotations.domain.attributes)))
    annotation_best_idx = np.argmax(annotations.X, axis=1)
    annotation_best = labels[annotation_best_idx]

    clusters_unique = set(
        clusters.domain[0].values) - {"-1"}  # -1 is not clustered
    annotations_clusters = {}
    for cl in clusters_unique:
        mask = np.array(list(
            map(clusters.domain.attributes[0].repr_val,
                clusters.X[:, 0]))).flatten() == cl
        labels_cl = annotation_best[mask]
        counts = Counter(labels_cl)
        annotations_clusters[cl] = [
            (l, c / len(labels_cl))
            for l, c in counts.most_common(labels_per_cluster)]

    return annotations_clusters


def labels_locations(coordinates, clusters):
    """
    Function computes the location of the label for each cluster.
    The location is compute as a center point.

    Parameters
    ----------
    coordinates : Orange.data.Table
        Data points
    clusters : Orange.data.Table
        Cluster indices for each item.

    Returns
    -------
    dict
        The coordinates for locating the label. Dictionary with cluster index
        as a key and tuple (x, y) as a value.
    """
    clusters_unique = set(
        clusters.domain[0].values) - {"-1"}  # -1 is not clustered
    locations = {}
    for cl in clusters_unique:
        mask = np.array(list(
            map(clusters.domain.attributes[0].repr_val,
                clusters.X[:, 0]))).flatten() == cl
        cl_coordinates = coordinates.X[mask, :]
        x, y = np.mean(cl_coordinates, axis=0)
        locations[cl] = (x, y)
    return locations


def get_epsilon(data, k=10, skip=0.1):
    """
    The function computes the epsilon parameter for DBSCAN through method
    proposed in the paper.

    Parameters
    ----------
    data : Orange.data.Table
        Input data which wil be clustered
    k : int
        Number kth observed neighbour
    skip : float
        Percentage of skipped neighborus.

    Returns
    -------
    float
        Epsilon parameter for DBSCAN
    """
    x = data.X
    if len(x) > 1000:  # subsampling is required
        i = len(x) // 1000
        x = x[::i]
    print(len(x))
    d = distance.squareform(distance.pdist(x))
    kth_point = np.argpartition(d, k+1, axis=1)[:, k+1]
    # k+1 since first one is item itself
    kth_dist = np.sort(d[np.arange(0, len(kth_point)), kth_point])

    # currently mark proportion equal to skip as a noise
    return kth_dist[-int(np.round(len(kth_dist) * skip))]

def annotate_projection(annotations, coordinates,
                        clustering_algorithm=DBSCAN,
                        labels_per_cluster=3, **kwargs):
    """
    Function cluster the data based on coordinates, and assigns a certain number
    of labels per cluster. Each cluster gets `labels_per_cluster` number of most
    common labels in cluster assigned.

    Parameters
    ----------
    annotations : Orange.data.Table
        Table with annotations and their probabilities.
    coordinates : Orange.data.Table
        Data to be clustered
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
    assert len(annotations) == len(coordinates)
    assert len(coordinates) > 0  # sklearn clustering want to have one example
    assert len(annotations.domain) > 0
    assert len(coordinates.domain) > 0
    # get clusters
    clusters = cluster_data(coordinates, clustering_algorithm, **kwargs)

    # assign top n labels to group
    annotations_cl = assign_labels(clusters, annotations, labels_per_cluster)

    labels_loc = labels_locations(coordinates, clusters)

    return clusters, annotations_cl, labels_loc
