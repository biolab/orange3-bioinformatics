from collections import defaultdict

import numpy as np
from Orange.data import Domain, ContinuousVariable, Table

from orangecontrib.bioinformatics.ncbi.gene import GeneInfo
from orangecontrib.bioinformatics.utils import statistics
from scipy.stats.mstats import rankdata
from scipy.stats import hypergeom

from orangecontrib.bioinformatics.widgets.utils.data import TAX_ID

SCORING_EXP_RATIO = "scoring_exp_ratio"
SCORING_MARKERS_SUM = "scoring_sum_of_expressed_markers"
SCORING_LOG_FDR = "scoring_log_fdr"
SCORING_LOG_PVALUE = "scoring_log_p_value"

PFUN_BINOMIAL = "binomial_p_function"
PFUN_HYPERGEOMETRIC = "HYPERGEOMETRIC_p_function"


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
    >>> selected_attributes, z = annotator.select_attributes(data)
    >>> scores, p_val = annotator.assign_annotations(
    ...     selected_attributes, markers, data.attributes[TAX_ID])
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
    def select_attributes(data, z_threshold=1):
        """
        Function selects "over"-expressed attributes for items with Mann-Whitney
        U test.

        Parameters
        ----------
        data : Orange.data.Table
            Tabular data with gene expressions
        z_threshold : float
            The threshold for selecting the attribute. For each item the
            attributes with z-value above this value are selected.

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
        attributes_sets = [
            set(map(str, set(attributes_np[row > z_threshold]) - {None}))
            for row in z]
        # map to string was added since there seems to be no guarantee that
        # Entrez ID is a string.

        # pack z values to data table
        # take only attributes in new domain
        domain = Domain([x for x in data.domain.attributes])
        z_table = Table(domain, z)

        return attributes_sets, z_table

    @staticmethod
    def _group_marker_attributes(markers, genes_order):
        """
        Function transforms annotations table to dictionary with format
        {annotation1: [attributes], annotation2: [attributes], ...}
        """
        # dictionary with structure {celltype: [gene1, gene2, ...], ...}
        types_dict = defaultdict(set)
        for m in markers:
            if m["Entrez ID"].value is not None and \
                    not m["Entrez ID"].value == "?":
                types_dict[str(m["Cell Type"])].add(m["Entrez ID"].value)
        # dictionary as list of items
        types_list = list(types_dict.items())

        # create numpy matrix with cell type - gene affiliation
        genes_celltypes = np.zeros(
            (len(genes_order), len(types_list)))
        for i, (_, genes) in enumerate(types_list):
            for g in genes:
                if g in genes_order:
                    genes_celltypes[genes_order.index(g), i] = 1
        return types_list, genes_celltypes

    @staticmethod
    def _scores_markers_sum(data, genes_types):
        return data.X.dot(genes_types)

    @staticmethod
    def _scores_fdr(fdrs):
        return -np.log(np.array(fdrs))

    @staticmethod
    def assign_annotations(items_sets, available_annotations, data,
                           p_value_fun=PFUN_BINOMIAL,
                           scoring=SCORING_EXP_RATIO):
        """
        The function gets a set of attributes (e.g. genes) for each cell and
        attributes for each annotation. It returns the annotations significant
        for each cell.

        Parameters
        ----------
        items_sets : list of sets
            Set of most important attributes for each item.
        available_annotations : Orange.data.Table
            Available annotations (e.g. cell types)
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
        tax_id = data.attributes[TAX_ID]

        # select function for p-value
        if p_value_fun == PFUN_HYPERGEOMETRIC:  # sf accept x-1 instead of x
            p_fun = lambda x, n, m, k: hypergeom.sf(x-1, n, m, k)
        else:
            p_fun = statistics.Binomial().p_value

        # retrieve number of genes for organism
        N = len(GeneInfo(tax_id))

        grouped_annotations_items, genes_celltypes = \
            AnnotateSamples._group_marker_attributes(
                available_annotations,
                [d.attributes.get("Entrez ID")
                 for d in data.domain.attributes])

        def hg_cell(item_attributes):
            p_values = []
            scores = []
            for i, (ct, attributes) in enumerate(grouped_annotations_items):
                intersect = item_attributes & attributes
                x = len(intersect)
                k = len(item_attributes)  # drawn balls - expressed for item
                m = len(attributes)  # marked balls - items for a process

                if x > 2:  # avoid the heavy computation when intersect small
                    p_value = p_fun(x, N, m, k)
                else:
                    p_value = 1
                p_values.append(p_value)

                if scoring == SCORING_EXP_RATIO:
                    scores.append(x / (m + 1e-16))

            fdrs = statistics.FDR(p_values)
            if scoring == SCORING_LOG_FDR or scoring == SCORING_LOG_PVALUE:
                scores = AnnotateSamples._scores_fdr(
                    fdrs if scoring == SCORING_LOG_FDR else p_values)

            return scores, fdrs

        prob_fdrs = [hg_cell(its) for its in items_sets]
        probs, fdrs = zip(*prob_fdrs)

        if scoring == SCORING_MARKERS_SUM:
            probs = AnnotateSamples._scores_markers_sum(data, genes_celltypes)

        domain = Domain(
            [ContinuousVariable(ct[0]) for ct in grouped_annotations_items])
        probs_table = Table(domain, np.array(probs))
        fdrs_table = Table(domain, np.array(fdrs))

        return probs_table, fdrs_table

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
        scores_x[p_values.X > p_threshold] = 0
        probs = Table(scores.domain, scores_x)

        if return_nonzero_annotations:
            col_nonzero = np.sum(probs, axis=0) > 0
            new_domain = Domain(
                np.array(scores.domain.attributes)[col_nonzero])
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

        selected_attributes, z = AnnotateSamples.select_attributes(
            data, z_threshold=z_threshold)
        annotation_probs, annotation_fdrs = AnnotateSamples.assign_annotations(
            selected_attributes, available_annotations, data,
            p_value_fun=p_value_fun, scoring=scoring)

        annotation_probs = AnnotateSamples.filter_annotations(
            annotation_probs, annotation_fdrs, return_nonzero_annotations,
            p_threshold
        )

        return annotation_probs
