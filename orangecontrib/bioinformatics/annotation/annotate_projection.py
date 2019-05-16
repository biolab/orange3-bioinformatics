"""
This module cluster the projection of data (usually it is 2D projection) with
one of the standard algorithms and attach a certain number of labels per
cluster.

Example:

>>> from Orange.projection import TSNE
>>> from Orange.data import Table
>>> from orangecontrib.bioinformatics.utils import serverfiles
>>> from orangecontrib.bioinformatics.annotation.annotate_projection import annotate_projection
>>> from orangecontrib.bioinformatics.annotation.annotate_samples import AnnotateSamples
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


def _cluster_data(coordinates, clustering_algorithm):
    """
    This function receives data and cluster them.

    Parameters
    ----------
    coordinates : Orange.data.Table
        Data to be clustered
    clustering_algorithm : callable
        Algorithm used in clustering.

    Returns
    -------
    ndarray
        List of cluster indices.
    """
    model = clustering_algorithm(coordinates)
    return model(coordinates).X.flatten()
    # TODO: this need to be changed when clustering in orange is changed


def _assign_labels(clusters, annotations, labels_per_cluster):
    """
    This function assigns a certain number of labels per cluster. Each cluster
    gets `labels_per_cluster` number of most common labels in cluster assigned.

    Parameters
    ----------
    clusters : ndarray
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

    clusters_unique = set(clusters) - {-1}  # -1 means that item not clustered
    annotations_clusters = {}
    for cl in clusters_unique:
        labels_cl = annotation_best[clusters == cl]
        counts = Counter(labels_cl)
        annotations_clusters[cl] = [
            (l, c / len(labels_cl))
            for l, c in counts.most_common(labels_per_cluster)]

    return annotations_clusters


def annotate_projection(annotations, coordinates,
                        clustering_algorithm=DBSCAN(),
                        labels_per_cluster=3):
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
    ndarray
        List of cluster indices.
    dict
        Dictionary with cluster index as a key and list of annotations as a
        value. Each list include tuples with the annotation name and their
        proportion in the cluster.
    """
    assert len(annotations) == len(coordinates)
    assert len(coordinates) > 0  # sklearn clustering want to have one example
    assert len(annotations.domain) > 0
    assert len(coordinates.domain) > 0
    # get clusters
    clusters = _cluster_data(coordinates, clustering_algorithm)

    # assign top n labels to group
    annotations_cl = _assign_labels(clusters, annotations, labels_per_cluster)
    return clusters, annotations_cl
