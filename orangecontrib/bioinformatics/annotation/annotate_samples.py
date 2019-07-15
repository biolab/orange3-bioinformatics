import numpy as np
from Orange.data import Domain, ContinuousVariable, Table

from orangecontrib.bioinformatics.ncbi.gene import GeneInfo
from orangecontrib.bioinformatics.utils import statistics
from scipy.stats.mstats import rankdata
from scipy.stats import hypergeom, binom

from orangecontrib.bioinformatics.widgets.utils.data import TAX_ID

SCORING_EXP_RATIO = "scoring_exp_ratio"
SCORING_MARKERS_SUM = "scoring_sum_of_expressed_markers"
SCORING_LOG_FDR = "scoring_log_fdr"
SCORING_LOG_PVALUE = "scoring_log_p_value"

PFUN_BINOMIAL = "binomial_p_function"
PFUN_HYPERGEOMETRIC = "HYPERGEOMETRIC_p_function"


class ScoringNotImplemented(Exception):
    pass


class AnnotateSamples:
    """
    AnnotateSamples is class used for the annotation of data items with the
    labels Mann-Whitney U test for selecting important values and
    the Hyper-geometric for assigning the labels.

    Example for full annotation:

    >>> from Orange.data import Table
    >>> from orangecontrib.bioinformatics.utils import serverfiles
    >>> from orangecontrib.bioinformatics.annotation.annotate_samples import \
    ...     AnnotateSamples
    >>> from orangecontrib.bioinformatics.widgets.utils.data import TAX_ID
    >>>
    >>> data = Table("https://datasets.orange.biolab.si/sc/aml-1k.tab.gz")
    >>> data.attributes[TAX_ID] = "9606"  # table needs to have an organism ID
    >>> markers_path = serverfiles.localpath_download(
    ...     'marker_genes','panglao_gene_markers.tab')
    >>> marker = Table(markers_path)
    >>>
    >>> # filter only human markers
    >>> from Orange.data.filter import FilterString, Values
    >>> f = FilterString("Organism", FilterString.Equal, "Human")
    >>> markers = Values([f])(marker)
    >>>
    >>> annotator = AnnotateSamples()
    >>> annotations = annotator.annotate_samples(
    ...     data, markers, p_threshold=0.05)

    Example for full manual annotation. Here annotation is split in three
    phases. We assume that data are already loaded.

    >>> annotator = AnnotateSamples()
    >>> z = annotator.mann_whitney_test(data)
    >>> scores, p_val = AnnotateSamples.assign_annotations(z, markers, data)
    >>> scores = annotator.filter_annotations(scores, p_val, p_threshold=0.05)

    Attributes
    ----------

    """
    @staticmethod
    def log_cpm(data):
        """
        Function normalizes data with CPM methods and normalize them.

        Parameters
        ----------
        data : Orange.data.Table
            Tabular data with gene expressions

        Returns
        -------
        Orange.data.Table
            Normalized gene expression data
        """
        norm_data = np.log(1 + AnnotateSamples._cpm(data.X))
        return Table(data.domain, norm_data)

    @staticmethod
    def _cpm(data):
        """
        This function normalizes data with CPM methods.

        Parameters
        ----------
        data : array_like
            Numpy array with data. Columns are genes, rows are cells.
        """
        return data / np.sum(data, axis=1)[:, None] * 1e6

    @staticmethod
    def _ranks(data):
        """
        This function computes ranks for data in the table along axis=0.
        """
        x_len = data.shape[0]
        x_mask = data.sum(axis=0) > 0

        # create a matrix of ranges - init with average rank
        # for columns without nonzero expressions
        data_ge_ranked = np.ones(data.shape) * (1 + data.shape[0]) / 2

        # compute ranks only for nonzero columns
        for i in np.where(x_mask)[0]:
            mask = data[:, i] > 0
            col = np.ones(x_len) * (1 + (x_len - mask.sum())) / 2
            col[mask] = rankdata(data[mask, i]) + (x_len - mask.sum())
            data_ge_ranked[:, i] = col
        return data_ge_ranked

    @staticmethod
    def mann_whitney_test(data):
        """
        Compute z values with test Mann-Whitney U test.

        Parameters
        ----------
        data : Orange.data.Table
            Tabular data with gene expressions

        Returns
        -------
        Orange.data.Table
            Z-value for each item.
        """
        if len(data.X) <= 1:
            return [], []

        # rank data
        data_ge_ranked = AnnotateSamples._ranks(data.X)

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

        # pack z values to data table
        # take only attributes in new domain
        domain = Domain([x for x in data.domain.attributes])
        z_table = Table(domain, z)

        return z_table

    @staticmethod
    def _reorder_matrix(matrix, genes_order):
        """
        Function reorder the columns of the array to fit to the genes_order

        Parameters
        ----------
        matrix : Orange.data.Table
            Tabular data tha needs to be reordered
        genes_order : list
            Desired genes order

        Returns
        ------
        np.ndarray
            Reordered array.
        """
        current_order = np.array(
            [x.attributes.get("Entrez ID")
             for x in matrix.domain.attributes])
        values = matrix.X

        # filter out genes without entrez ID
        has_entrez_id = current_order != None
        # just in case if Entrez ID not strings
        current_order = np.array([str(x) for x in current_order[has_entrez_id]])
        values = values[:, has_entrez_id]

        genes_order = np.array(genes_order)

        xsorted = np.argsort(genes_order)
        ypos = np.searchsorted(genes_order[xsorted], current_order)
        indices = xsorted[ypos]  # index which tell where should be the column

        reordered_values = np.zeros((values.shape[0], len(genes_order)))
        for i_curr, i_dest in enumerate(indices):
            reordered_values[:, i_dest] = values[:, i_curr]

        return reordered_values

    @staticmethod
    def _select_attributes(z, genes_order, z_threshold=1):
        """
        Function selects "over"-expressed attributes for items based on z
        values. It also reorder the matrix columns.

        Parameters
        ----------
        z : Orange.data.Table
            Tabular data z values for each item in the table
        genes_order : list
            Desired genes order
        z_threshold : float
            The threshold for selecting the attribute. For each item the
            attributes with z-value above this value are selected.

        Returns
        -------
        np.ndarray
            Reordered and thresholded z-values/
        """
        reordered_z = AnnotateSamples._reorder_matrix(z, genes_order)

        return reordered_z > z_threshold

    @staticmethod
    def _group_marker_attributes(markers, genes_order):
        """
        Function transforms annotations table to dictionary with format
        {annotation1: [attributes], annotation2: [attributes], ...}
        """
        types = sorted(list(set(markers[:, "Cell Type"].metas.flatten())))
        genes_celltypes = np.zeros((len(genes_order), len(types)))

        for m in markers:
            g = str(m["Entrez ID"].value)
            m = m["Cell Type"].value
            if g is not None and not g == "?":
                genes_celltypes[genes_order.index(g), types.index(m)] = 1

        return genes_celltypes, types

    @staticmethod
    def _score(scoring_type, p_values, fdrs, data, M, x, m, genes_order):
        if scoring_type == SCORING_MARKERS_SUM:
            return AnnotateSamples._reorder_matrix(data, genes_order).dot(M)
        if scoring_type == SCORING_EXP_RATIO:
            return x / m
        if scoring_type == SCORING_LOG_FDR:
            return -np.log(fdrs)
        if scoring_type == SCORING_LOG_PVALUE:
            return -np.log(p_values)
        else:
            raise ScoringNotImplemented()

    @staticmethod
    def assign_annotations(z_values, available_annotations, data, z_threshold=1,
                           p_value_fun=PFUN_BINOMIAL,
                           scoring=SCORING_EXP_RATIO):
        """
        The function gets a set of attributes (e.g. genes) for each cell and
        attributes for each annotation. It returns the annotations significant
        for each cell.

        Parameters
        ----------
        z_values : Orange.data.Table
            Table which show z values for each item
        available_annotations : Orange.data.Table
            Available annotations (e.g. cell types)
        z_threshold : float
            The threshold for selecting the attribute. For each item the
            attributes with z-value above this value are selected.
        p_value_fun : str, optional (defaults: TEST_BINOMIAL)
            A function that calculates p-value. It can be either
            PFUN_BINOMIAL that uses statistics.Binomial().p_value or
            PFUN_HYPERGEOMETRIC that uses hypergeom.sf.
        data : Orange.data.Table
            Tabular data with gene expressions - we need that to compute scores.
        scoring : str, optional (default=SCORING_EXP_RATIO)
            Type of scoring

        Returns
        -------
        Orange.data.Table
            Annotation probabilities
        Orange.data.Table
            Annotation fdrs
        """
        assert TAX_ID in data.attributes, "The input table needs to have a " \
                                          "tax_id attribute"
        assert any(
            "Entrez ID" in x.attributes for x in data.domain.attributes),\
            "Input data do not contain gene expression data."
        tax_id = data.attributes[TAX_ID]

        # select function for p-value
        if p_value_fun == PFUN_HYPERGEOMETRIC:
            p_fun = lambda x, N, m, k: hypergeom.sf(x, N, m, k)
        else:
            p_fun = lambda x, N, m, k: binom.sf(x, k, m / N)

        N = len(GeneInfo(tax_id))  # number of genes for organism

        # make an attributes order
        genes_data = [str(x.attributes["Entrez ID"])
                      for x in z_values.domain.attributes
                      if "Entrez ID" in x.attributes]
        genes_celltypes = [
            str(x) for x in available_annotations[:, "Entrez ID"].metas.flatten()
            if x is not None and not x == "?"]
        genes_order = list(set(genes_data) | set(genes_celltypes))

        # get marker genes matrix M
        M, annotations = AnnotateSamples._group_marker_attributes(
            available_annotations, genes_order)

        Z = AnnotateSamples._select_attributes(
            z_values, genes_order, z_threshold)

        x = Z.dot(M)
        k = np.repeat(Z.sum(axis=1).reshape(-1, 1), x.shape[1], axis=1)
        m = np.repeat(M.sum(axis=0).reshape(1, -1), x.shape[0], axis=0)

        p_values = p_fun(x - 1, N, m, k)

        fdrs = np.zeros(p_values.shape)
        for i, row in enumerate(p_values):
            fdrs[i] = np.array(statistics.FDR(row.tolist()))

        scores = AnnotateSamples._score(
            scoring, p_values, fdrs, data, M, x, m, genes_order)

        domain = Domain(
            [ContinuousVariable(ct) for ct in annotations])
        scores_table = Table(domain, scores)
        fdrs_table = Table(domain, fdrs)

        return scores_table, fdrs_table

    @staticmethod
    def filter_annotations(scores, p_values, return_nonzero_annotations=True,
                           p_threshold=0.05):
        """
        This function filters the probabilities on places that do not reach the
        threshold for p-value and filter zero columns
        return_nonzero_annotations is True.

        Parameters
        ----------
        scores : Orange.data.Table
            Scores for each annotations for each cell
        p_values : Orange.data.Table
            p-value scores for annotations for each cell
        return_nonzero_annotations : bool
            Flag that enables filtering the non-zero columns.
        p_threshold : float
            A threshold for accepting the annotations. Annotations that has FDR
            value bellow this threshold are used.


        Returns
        -------
        Orange.data.Table
            Filtered scores for each annotations for each cell
        """
        scores_x = np.copy(scores.X)  # do not want to edit values inplace
        scores_x[p_values.X > p_threshold] = np.nan
        scores = Table(scores.domain, scores_x)

        if return_nonzero_annotations:
            col_not_empty = ~np.isnan(scores).all(axis=0)
            new_domain = Domain(
                np.array(scores.domain.attributes)[col_not_empty])
            scores = Table(new_domain, scores)
        return scores

    @staticmethod
    def annotate_samples(data, available_annotations,
                         return_nonzero_annotations=True, p_threshold=0.05,
                         p_value_fun=PFUN_BINOMIAL, z_threshold=1,
                         scoring=SCORING_EXP_RATIO, normalize=False):
        """
        Function marks the data with annotations that are provided. This
        function implements the complete functionality. First select genes,
        then annotate them and filter them.

        Parameters
        ----------
        data : Orange.data.Table
            Tabular data
        available_annotations : Orange.data.Table
            Available annotations (e.g. cell types)
        return_nonzero_annotations : bool, optional (default=True)
            If true return scores for only annotations present in at least one
            sample.
        p_threshold : float
            A threshold for accepting the annotations. Annotations that has FDR
            value bellow this threshold are used.
        p_value_fun : str, optional (defaults: TEST_BINOMIAL)
            A function that calculates p-value. It can be either
            PFUN_BINOMIAL that uses statistics.Binomial().p_value or
            PFUN_HYPERGEOMETRIC that uses hypergeom.sf.
        z_threshold : float
            The threshold for selecting the attribute. For each item the
            attributes with z-value above this value are selected.
        scoring : str, optional (default = SCORING_EXP_RATIO)
            Type of scoring
        normalize : bool, optional (default = False)
            This variable tells whether to normalize data or not.

        Returns
        -------
        Orange.data.Table
            Cell type most important for each cell.
        """
        assert len(data) > 1, "At least two data items are required for " \
                              "method to work."

        if normalize:
            data = AnnotateSamples.log_cpm(data)

        z = AnnotateSamples.mann_whitney_test(
            data)

        annotation_probs, annotation_fdrs = AnnotateSamples.assign_annotations(
            z, available_annotations, data, z_threshold=z_threshold,
            p_value_fun=p_value_fun, scoring=scoring)

        annotation_probs = AnnotateSamples.filter_annotations(
            annotation_probs, annotation_fdrs, return_nonzero_annotations,
            p_threshold
        )

        return annotation_probs
