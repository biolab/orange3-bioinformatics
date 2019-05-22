from collections import defaultdict

import numpy as np
from Orange.data import Domain, ContinuousVariable, Table

from orangecontrib.bioinformatics.utils import statistics
from scipy.stats.mstats import rankdata
from scipy.stats import hypergeom


class AnnotateSamples:
    """
    This class annotate data items with the labels using Hyper-geometric test

    Attributes
    ----------
    p_value_th : float
        A threshold for the FDR. FDR values bellow this value are accepted
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
    def select_genes(data):
        """
        This function selects "over"-expressed genes for cells
        with Mann-Whitney U test.

        Parameters
        ----------
        data : Orange.data.Table
            Gene expressions

        Returns
        -------
        :obj:`list`
            Sets of selected genes for each cell
        """
        if len(data.X) <= 0:
            return [], []
        # rank data
        data_ge_ranked = rankdata(data.X, axis=0)

        # compute U, mu, sigma
        n2 = data_ge_ranked.shape[0]
        n = n2 + 1
        u = data_ge_ranked - 1
        mu = n2 / 2
        sigma = np.zeros(data_ge_ranked.shape[1])
        for i in range(data_ge_ranked.shape[1]):
            _, counts = np.unique(data_ge_ranked[:, i], return_counts=True)
            sigma[i] = np.sqrt(
                1 * n2 / 12 * ((n + 1) - np.sum((counts ** 3 - counts)) /
                               (n * (n - 1))))

        # compute z
        z = (u - mu) / sigma

        # gene selection
        genes_np = np.array([
            a.attributes.get("Entrez ID") for a in data.domain.attributes])
        ge_expressed_sets = [set(genes_np[row > 1]) - {None} for row in z]
        return ge_expressed_sets, z

    @staticmethod
    def group_marker_genes(markers):
        types_genes_dict = defaultdict(set)
        for m in markers:
            types_genes_dict[str(m["Cell Type"])].add(int(m["Entrez ID"].value))
        return types_genes_dict

    def assign_annotations(self, gene_sets, cell_types_markers):
        """
        Function get set of genes that represents each cell and marker genes for
        each cell type. It returns the cell type most significant for each cell.

        Parameters
        ----------
        gene_sets : list of sets
            Set of most expressed genes for each cell.
        cell_types_markers : dict
            marker genes for each cell type

        Returns
        -------
            Cell type most important for each cell.
        """
        N = 20365  # number of all genes for human

        marker_genes = self.group_marker_genes(cell_types_markers)
        marker_genes_items = list(marker_genes.items())

        def hg_cell(cell_genes):
            p_values = []
            prob = np.zeros(len(marker_genes_items))
            for i, (ct, genes) in enumerate(marker_genes_items):
                x = len(cell_genes & genes)
                k = len(cell_genes)  # drawn balls - features expressed for item
                m = len(genes)  # number of marked balls - items for a process

                p_value = self.p_value_fun(x, N, m, k)

                p_values.append(p_value)
                prob[i] = x / (m + 1e-16)

            fdrs = statistics.FDR(p_values)
            prob[np.array(fdrs) > self.p_threshold] = 0
            return prob

        return np.array([hg_cell(cg) for cg in gene_sets]), \
            [x[0] for x in marker_genes_items]

    def annotate_samples(self, data, marker_genes):
        selected_genes, z = self.select_genes(data)
        cell_types_scores, cell_types = self.assign_annotations(
            selected_genes, marker_genes)

        domain = Domain([ContinuousVariable(ct) for ct in cell_types])
        if len(cell_types_scores) <= 0:
            cell_types_scores = np.empty((0, len(domain)))
        cell_types_scores_table = Table(domain, cell_types_scores)
        return cell_types_scores_table
