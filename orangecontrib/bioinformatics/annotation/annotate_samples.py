from collections import defaultdict

import numpy as np
from Orange.data import Domain, ContinuousVariable, Table

from orangecontrib.bioinformatics.utils import statistics
from scipy.stats.mstats import rankdata
from scipy.stats import hypergeom


class AnnotateSamples:
    """
    AnnotateSamples is class used for the annotation of data items with the
    labels Mann-Whitney U test for selecting important values and
    the Hyper-geometric for assigning the labels.

    Example of use:

    >>> from Orange.data import Table
    >>> from orangecontrib.bioinformatics.utils import serverfiles
    >>> from orangecontrib.bioinformatics.annotation.annotate_samples import AnnotateSamples
    >>>
    >>> data = Table("https://datasets.orange.biolab.si/sc/aml-1k.tab.gz")
    >>> markers_path = serverfiles.localpath_download(
    >>>     'marker_genes','panglao_gene_markers.tab')
    >>> marker = Table(markers_path)
    >>> annotator = AnnotateSamples(p_value_th=0.05)
    >>> annotations = annotator.annotate_samples(data, markers)

    Attributes
    ----------
    p_value_th : float
        A threshold for accepting the annotations. Annotations that has FDR
        value bellow this threshold are used.
    p_value_fun : callable, optional (defaults: statistics.Binomial().p_value)
        A function that calculates p-value. It can be either
        statistics.Binomial().p_value or hypergeom.sf.
    """

    def __init__(self, p_value_th, p_value_fun=statistics.Binomial().p_value):
        if p_value_fun == hypergeom.sf:  # because sf accept x-1 instead of x
            self.p_value_fun = lambda x, n, m, k: p_value_fun(x-1, n, m, k)
        else:
            self.p_value_fun = p_value_fun

        self.p_threshold = p_value_th

    @staticmethod
    def _select_genes(data):
        """
        Function selects "over"-expressed attributes for items with Mann-Whitney
        U test.

        Parameters
        ----------
        data : Orange.data.Table
            Tabular data

        Returns
        -------
        :obj:`list`
            Sets of selected attributes for each cell
        """
        if len(data.X) <= 1:
            return [], []
        # rank data
        data_ge_ranked = rankdata(data.X, axis=0)

        # compute U, mu, sigma
        n = data_ge_ranked.shape[0]
        n2 = n - 1
        u = data_ge_ranked - 1
        mu = n2 / 2
        sigma = np.zeros(data_ge_ranked.shape[1])
        for i in range(data_ge_ranked.shape[1]):
            _, counts = np.unique(data_ge_ranked[:, i], return_counts=True)
            sigma[i] = np.sqrt(
                1 * n2 / 12 * ((n + 1) - np.sum((counts ** 3 - counts)) /
                               (n * (n - 1))))

        # compute z
        z = (u - mu) / (sigma + 1e-16)

        # gene selection
        attributes_np = np.array([
            a.attributes.get("Entrez ID") for a in data.domain.attributes])
        attributes_sets = [set(attributes_np[row > 1]) - {None} for row in z]
        return attributes_sets, z

    @staticmethod
    def _group_marker_attributes(markers):
        """
        Function transforms annotations table to dictionary with format
        {annotation1: [attributes], annotation2: [attributes], ...}
        """
        types_dict = defaultdict(set)
        for m in markers:
            types_dict[str(m["Cell Type"])].add(int(m["Entrez ID"].value))
        return types_dict

    def _assign_annotations(self, items_sets, available_annotations):
        """
        Function get set of attributes (e.g. genes) that represents each item
        and attributes for each annotation. It returns the annotations most
        significant for each cell.

        Parameters
        ----------
        items_sets : list of sets
            Set of most important attributes for each item.
        available_annotations : Orange.data.Table
            Available annotations (e.g. cell types)

        Returns
        -------
        ndarray
            Annotation probabilities
        list
            Annotation list
        """
        N = 20365  # number of all genes for human

        grouped_annotations = self._group_marker_attributes(
            available_annotations)
        grouped_annotations_items = list(grouped_annotations.items())

        def hg_cell(item_attributes):
            p_values = []
            prob = np.zeros(len(grouped_annotations_items))
            for i, (ct, attributes) in enumerate(grouped_annotations_items):
                x = len(item_attributes & attributes)
                k = len(item_attributes)  # drawn balls - expressed for item
                m = len(attributes)  # marked balls - items for a process

                p_value = self.p_value_fun(x, N, m, k)

                p_values.append(p_value)
                prob[i] = x / (m + 1e-16)

            fdrs = statistics.FDR(p_values)
            prob[np.array(fdrs) > self.p_threshold] = 0
            return prob

        return np.array([hg_cell(cg) for cg in items_sets]), \
            [x[0] for x in grouped_annotations_items]

    def annotate_samples(self, data, available_annotations):
        """
        Function marks the data with annotations that are provided provided.

        Parameters
        ----------
        data : Orange.data.Table
            Tabular data
        available_annotations : Orange.data.Table
            Available annotations (e.g. cell types)

        Returns
        -------
        Orange.data.Table
            Cell type most important for each cell.
        """
        assert len(data) > 1, "At least two data items are required for " \
                              "method to work."

        selected_attributes, z = self._select_genes(data)
        annotations_scores, annotations = self._assign_annotations(
            selected_attributes, available_annotations)

        domain = Domain([ContinuousVariable(ct) for ct in annotations])
        if len(annotations_scores) <= 0:
            annotations_scores = np.empty((0, len(domain)))
        scores_table = Table(domain, annotations_scores)
        return scores_table
