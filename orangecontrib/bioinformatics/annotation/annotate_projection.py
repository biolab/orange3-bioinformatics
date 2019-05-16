from collections import Counter

from Orange.clustering import DBSCAN, KMeans
import numpy as np


def cluster_data(coordinates, clustering_algorithm):
    model = clustering_algorithm.fit(coordinates.X)  # popravi
    return model(coordinates.X).flatten()


def assign_labels(clusters, annotations, labels_per_cluster):
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
    # get clusters
    clusters = cluster_data(coordinates, clustering_algorithm)

    # assign top n labels to group
    annotations_cl = assign_labels(clusters, annotations, labels_per_cluster)
    return clusters, annotations_cl



