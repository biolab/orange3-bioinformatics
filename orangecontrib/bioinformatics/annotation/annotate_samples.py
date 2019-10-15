import numpy as np
import pandas as pd
from scipy.sparse import issparse
from pointannotator.annotate_samples import (
    PFUN_BINOMIAL,
    SCORING_LOG_FDR,
    SCORING_EXP_RATIO,
    PFUN_HYPERGEOMETRIC,
    SCORING_MARKERS_SUM,
    AnnotateSamples,
)

from Orange.data import Table, Domain, ContinuousVariable

from orangecontrib.bioinformatics.ncbi.gene import GeneInfo
from orangecontrib.bioinformatics.widgets.utils.data import TAX_ID

__all__ = [
    "PFUN_BINOMIAL",
    "PFUN_HYPERGEOMETRIC",
    "SCORING_EXP_RATIO",
    "SCORING_LOG_FDR",
    "SCORING_MARKERS_SUM",
    "AnnotateSamplesMeta",
]


class ScoringNotImplemented(Exception):
    pass


class AnnotateSamplesMeta:
    """
    AnnotateSamples is a meta class used for the annotation of data items with
    the labels Mann-Whitney U test for selecting important values and
    the Hyper-geometric for assigning the labels. This class is a proxy to the
    external library - point-annotator.

    Example for full annotation:

    >>> from Orange.data import Table
    >>> from orangecontrib.bioinformatics.utils import serverfiles
    >>> from orangecontrib.bioinformatics.annotation.annotate_samples import \
    ...     AnnotateSamples
    >>> from orangecontrib.bioinformatics.widgets.utils.data import TAX_ID
    >>> from Orange.data.filter import FilterString, Values
    >>>
    >>> data = Table("https://datasets.orange.biolab.si/sc/aml-1k.tab.gz")
    >>> data.attributes[TAX_ID] = "9606"  # table needs to have an organism ID
    >>> markers_path = serverfiles.localpath_download(
    ...     'marker_genes','panglao_gene_markers.tab')
    >>> marker = Table(markers_path)
    >>>
    >>> # filter only human markers
    >>> f = FilterString("Organism", FilterString.Equal, "Human")
    >>> markers = Values([f])(marker)
    >>>
    >>> annotator = AnnotateSamples
    >>> z = annotator.mann_whitney_test(data)
    >>> scores, p_val = AnnotateSamples.assign_annotations(z, markers, data)
    >>> scores = annotator.filter_annotations(scores, p_val, p_threshold=0.05)
    """

    @staticmethod
    def _to_pandas(data, use_entrez_id=False):
        """
        Transform Orange table to pandas dataframe.

        Parameters
        ----------
        data : Orange.data.Table
            Orange table with gene data. We suppose that genes have `Entrez ID`
            assigned.
        use_entrez_id : bool
            This parameter tells whether to use Entrez ID as column names.
            When it is False variables names will be used

        Returns
        -------
        pd.DataFrame
            Pandas DataFrame with EntrezID
        list
            List of booleans with True on location of included columns in
            dataframe.
        """
        if use_entrez_id:
            entrez_ids = [x.attributes.get("Entrez ID") for x in data.domain.attributes]
            has_entrez_id = np.array([x is not None for x in entrez_ids])
            columns = list(map(str, np.array(entrez_ids)))
            columns_subset = has_entrez_id.tolist()
        else:
            columns = list(map(str, data.domain.attributes))
            columns_subset = None

        if issparse(data.X):
            df = pd.DataFrame.sparse.from_spmatrix(data.X, columns=columns)
        else:
            df = pd.DataFrame(data.X, columns=columns)
        return (df if columns_subset is None else df.loc[:, columns_subset], columns_subset)

    @staticmethod
    def mann_whitney_test(data):
        """
        Compute z values with Mann-Whitney U test.

        Parameters
        ----------
        data : Orange.data.Table
            Tabular data with gene expressions

        Returns
        -------
        Orange.data.Table
            Z-value for each item.
        """
        df, columns = AnnotateSamplesMeta._to_pandas(data, use_entrez_id=False)
        df_z_values = AnnotateSamples.mann_whitney_test(df)

        z_table = Table(data.domain, df_z_values.values, data.Y, metas=data.metas)

        return z_table

    @staticmethod
    def _prepare_annotations(ann):
        """
        This function prepares annotation to be acceptable by the method.
        """
        # transform to string if discrete else keep string
        ct_values = (
            list(map(ann.domain["Cell Type"].repr_val, ann.get_column_view("Cell Type")[0]))
            if ann.domain["Cell Type"].is_discrete
            else ann.get_column_view("Cell Type")[0]
        )
        entrez_values = (
            list(map(ann.domain["Entrez ID"].repr_val, ann.get_column_view("Entrez ID")[0]))
            if ann.domain["Entrez ID"].is_discrete
            else ann.get_column_view("Entrez ID")[0]
        )

        # the framework recognizes Gene instead of Entrez ID
        df_available_annotations = pd.DataFrame(list(zip(ct_values, entrez_values)), columns=["Cell Type", "Gene"])
        return df_available_annotations

    @staticmethod
    def assign_annotations(
        z_values, available_annotations, data, z_threshold=1, p_value_fun=PFUN_BINOMIAL, scoring=SCORING_EXP_RATIO
    ):
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
        # checks that assures that data are ok
        assert TAX_ID in data.attributes, "The input table needs to have a " "tax_id attribute"
        assert any(
            "Entrez ID" in x.attributes for x in data.domain.attributes
        ), "Input data do not contain gene expression data."

        # retrieve number of genes
        tax_id = data.attributes[TAX_ID]
        n = len(GeneInfo(tax_id))  # number of genes for organism

        # transform data to pandas dataframe
        df_z_values, _ = AnnotateSamplesMeta._to_pandas(z_values, use_entrez_id=True)
        df_data, _ = AnnotateSamplesMeta._to_pandas(data, use_entrez_id=True)
        df_available_annotations = AnnotateSamplesMeta._prepare_annotations(available_annotations)
        df_available_annotations = df_available_annotations[df_available_annotations["Gene"] != "?"]

        # call the method
        scores, fdrs = AnnotateSamples.assign_annotations(
            df_z_values,
            df_available_annotations,
            df_data,
            num_all_attributes=n,
            attributes_col="Gene",
            annotations_col="Cell Type",
            z_threshold=z_threshold,
            p_value_fun=p_value_fun,
            scoring=scoring,
        )
        # create orange tables
        domain = Domain([ContinuousVariable(ct) for ct in scores.columns.values])
        scores_table = Table(domain, scores.values)
        fdrs_table = Table(domain, fdrs.values)

        return scores_table, fdrs_table

    @staticmethod
    def filter_annotations(scores, p_values, return_nonzero_annotations=True, p_threshold=0.05):
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
        df_scores, _ = AnnotateSamplesMeta._to_pandas(scores)
        df_p_values, _ = AnnotateSamplesMeta._to_pandas(p_values)
        scores = AnnotateSamples.filter_annotations(df_scores, df_p_values, return_nonzero_annotations, p_threshold)

        domain = Domain([ContinuousVariable(ct) for ct in scores.columns.values])
        scores_table = Table(domain, scores.values)
        return scores_table
