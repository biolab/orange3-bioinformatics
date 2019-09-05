""" Cluster analysis module """
import threading
import concurrent.futures
from typing import Union
from operator import attrgetter
from functools import partial

import numpy as np

from AnyQt.QtCore import Qt, Slot, QThread, QVariant, QAbstractListModel

from Orange.widgets.gui import ProgressBar
from Orange.widgets.utils.concurrent import FutureWatcher, ThreadExecutor, methodinvoke

from orangecontrib.bioinformatics.geneset import GeneSet
from orangecontrib.bioinformatics.ncbi.gene import Gene
from orangecontrib.bioinformatics.utils.statistics import FDR, ALT_GREATER
from orangecontrib.bioinformatics.widgets.utils.gui import gene_scoring_method

DISPLAY_GENE_COUNT = 20
DISPLAY_GENE_SETS_COUNT = 5


class ClusterGene(Gene):
    __slots__ = ['score', 'p_val', 'fdr']

    def __init__(self, gene_name, gene_id):
        super().__init__()
        self.input_identifier = gene_name
        self.gene_id = gene_id

        # default values
        self.fdr = 1
        self.p_val = 1
        self.score = 0


class ClusterGeneSet(GeneSet):
    __slots__ = ['p_val', 'fdr', 'count']

    def __init__(self):
        super().__init__()


class Cluster:

    CLUSTER_VS_REST = False
    CLUSTER_VS_CLUSTER = True

    def __init__(self, name, index):
        # type: (str, int) -> None
        self.name = name
        self.index = index

        # gene enrichment
        self.genes = []
        self.filtered_genes = []

        # set enrichment
        self.gene_sets = []
        self.filtered_gene_sets = []

        # Statistical method used when performing analysis
        self.method_used = None

    def set_genes(self, gene_names, gene_ids):
        self.genes = []

        for name, entrez in zip(gene_names, gene_ids):
            self.genes.append(ClusterGene(name, entrez))

    @staticmethod
    def apply_filter(obj, p_val, fdr):
        filter_status = []

        if p_val is not None:
            filter_status.append(obj.p_val < p_val)

        if fdr is not None:
            filter_status.append(obj.fdr < fdr)

        return all(filter_status)

    def filter_enriched_genes(self, p_val, fdr, max_gene_count=None):
        all_args_none = all(arg is None for arg in [p_val, fdr])
        filter_function = partial(self.apply_filter, p_val=p_val, fdr=fdr)
        sorted_list = sorted(self.genes, key=attrgetter('p_val' if all_args_none else 'fdr'))
        filtered_list = list(filter(filter_function, sorted_list))

        if max_gene_count is not None:
            self.filtered_genes = filtered_list[:max_gene_count]
        else:
            self.filtered_genes = filtered_list

    def filter_gene_sets(self, p_val, fdr, max_set_count=None):
        all_args_none = all(arg is None for arg in [p_val, fdr])
        filter_function = partial(self.apply_filter, p_val=p_val, fdr=fdr)
        sorted_list = sorted(self.gene_sets, key=attrgetter('p_val' if all_args_none else 'fdr'))
        filtered_list = list(filter(filter_function, sorted_list))

        if max_set_count is not None:
            self.filtered_gene_sets = filtered_list[:max_set_count]
        else:
            self.filtered_gene_sets = filtered_list

    def gene_set_enrichment(self, gene_sets, selected_sets, genes, ref_genes):
        self.gene_sets = []

        if not genes:
            return

        # calculate gene set enrichment
        for gene_set in gene_sets:

            if gene_set.hierarchy not in selected_sets:
                continue

            gs = ClusterGeneSet()
            enrichment_result = gene_set.set_enrichment(ref_genes, genes.intersection(genes))
            gs.count = len(enrichment_result.query)
            gs.p_val = enrichment_result.p_value
            gs.name = gene_set.name
            gs.gs_id = gene_set.gs_id
            self.gene_sets.append(gs)

        # calculate FDR
        p_vals = [gs.p_val for gs in self.gene_sets]
        fdr_values = FDR(p_vals)

        for gs, fdr in zip(self.gene_sets, fdr_values):
            gs.fdr = fdr

    def __update_gene_objects(self, scores, p_vals, fdr_vals):
        # type: (Union[np.ndarray, list], Union[np.ndarray, list], Union[np.ndarray, list]) ->  None
        """ update gene objects with computed results

        :param p_vals:   Computed p-values
        :param fdr_vals: Computed fdr-values

        """
        for score, p_val, fdr, gene in zip(scores, p_vals, fdr_vals, self.genes):
            gene.score = score
            gene.p_val = p_val
            gene.fdr = fdr

    def cluster_scores(self, table_x, rows_by_cluster, method, design, **kwargs):
        # type: (np.ndarray, np.ndarray, gene_scoring_method, str) -> None
        """
        General scoring of genes in the cluster.
        Ways to score genes are determined by design.
        If a batch variable index is defined, it is accounted for in gene scoring.

        :param table_x:
        :param rows_by_cluster:
        :param method:
        :param design:
        :param kwargs:
        :param aggregation:
        :param rows_by_batch:
        :return:
        """
        aggregation = kwargs.get('aggregation', 'max')
        alternative = kwargs.get('alternative', ALT_GREATER)
        rows_by_batch = kwargs.get('rows_by_batch', None)
        if not isinstance(rows_by_batch, np.ndarray):
            rows_by_batch = np.zeros((len(table_x),))
        uniq_batches = set(rows_by_batch)

        # Determine clusters
        self.method_used = method.name
        uniq_clusters = set(rows_by_cluster) - {self.index}
        this_cluster = self.index
        if design == self.CLUSTER_VS_REST:
            rows_by_cluster = rows_by_cluster == self.index
            uniq_clusters = {False}
            this_cluster = True

        calculated_p_values = np.ones(
            (table_x.shape[1], len(uniq_clusters), len(uniq_batches))  # genes  # other clusters
        )  # batches

        calculated_scores = np.ones(
            (table_x.shape[1], len(uniq_clusters), len(uniq_batches))  # genes  # other clusters
        )  # batches

        for bi, b in enumerate(uniq_batches):
            for ci, c in enumerate(uniq_clusters):
                cluster = table_x[np.logical_and(rows_by_cluster == this_cluster, rows_by_batch == b)]
                rest = table_x[np.logical_and(rows_by_cluster == c, rows_by_batch == b)]
                if cluster.any() and rest.any():
                    scores, p_values = method.score_function(cluster, rest, alternative=alternative)
                    scores[np.isnan(p_values)] = 0
                    calculated_scores[:, ci, bi] = scores
                    p_values[np.isnan(p_values)] = 1
                    calculated_p_values[:, ci, bi] = p_values

        if aggregation == 'max':
            max_p_values = np.max(calculated_p_values, axis=(1, 2))
            max_p_indexes = np.where(np.max(calculated_p_values, axis=(1, 2), keepdims=True) == calculated_p_values)
            # this holds true only if max_p_indexes.ndim == 3
            scores = calculated_scores[max_p_indexes[0], max_p_indexes[1], max_p_indexes[2]]
            fdr_values = FDR(max_p_values.tolist())
            self.__update_gene_objects(scores, max_p_values, fdr_values)
        else:
            raise NotImplementedError("Aggregation %s is not implemented" % aggregation)
        return

    def to_html(self):
        gene_sets = '(no enriched gene sets)'
        if self.filtered_gene_sets:
            sets_to_display = self.filtered_gene_sets[:DISPLAY_GENE_SETS_COUNT]
            gene_sets = '<br>'.join(
                [
                    '<b>{}</b> (FDR={:0.2e}, n={})'.format(g_set.name, g_set.fdr, g_set.count)
                    for g_set in sets_to_display
                ]
            )

            if len(self.filtered_gene_sets) > len(sets_to_display):
                gene_sets += '<br> ... ({} more gene sets)'.format(
                    len(self.filtered_gene_sets) - DISPLAY_GENE_SETS_COUNT
                )

        genes = '(all genes are filtered out)'
        if self.filtered_genes:
            genes_to_display = [gene.input_identifier for gene in self.filtered_genes[:DISPLAY_GENE_COUNT]]

            genes = ', '.join(genes_to_display)
            if len(self.filtered_genes) > len(genes_to_display):
                genes += ', ... ({} more genes)'.format(len(self.filtered_genes) - DISPLAY_GENE_COUNT)

        html_string = """
        <html>
        <body>
        <table >
            <tr>
            </tr>
            <tr>
                <td width=\"20%\" heigth=\"100%\" align=left valign=top><p>{}</p></td>
                <td width=\"80%\" heigth=\"100%\" align=left valign=top>
                    {}
                </td>
            </tr>
            <tr> </tr>
            <tr>
                <td width=\"20%\" ></td>
                <td width=\"80%\" align=left valign=top>{}</td>
            </tr>
            <tr>
            </tr>
        </table>
        </body>
        </html>
        """.format(
            self.name, gene_sets, genes
        )
        return html_string


class Task:
    future = None
    watcher = None
    cancelled = False

    def cancel(self):
        self.cancelled = True
        self.future.cancel()
        concurrent.futures.wait([self.future])


class ClusterModel(QAbstractListModel):
    def __init__(self, parent=None):
        QAbstractListModel.__init__(self)
        self.__items = []
        self.parent = parent

        self._task = None  # type: Union[Task, None]
        self._executor = ThreadExecutor()

    def add_rows(self, rows):
        self.__items = rows

    def get_rows(self):
        return self.__items

    def rowCount(self, *args, **kwargs):
        return len(self.__items)

    def data(self, model_index, role=None):
        # check if data is set
        if not self.__items:
            return QVariant()

        # return empty QVariant if model index is unknown
        if not model_index.isValid() or not (0 <= model_index.row() < len(self.__items)):
            return QVariant()

        row_obj = self.__items[model_index.row()]

        if role == Qt.DisplayRole:
            return row_obj

    @Slot(concurrent.futures.Future)
    def _end_task(self, f):
        assert self.thread() is QThread.currentThread()
        assert threading.current_thread() == threading.main_thread()
        assert self._task is not None
        assert self._task.future is f
        assert f.done()

        self._task = None
        self.parent.progressBarFinished()
        self.parent.filter_genes()

        try:
            f.result()
        except Exception as ex:
            raise ex

    def _score_genes(self, callback, **kwargs):
        for item in self.get_rows():
            item.cluster_scores(**kwargs)
            callback()

    @Slot(bool)
    def progress_advance(self, finish):
        # GUI should be updated in main thread. That's why wex are calling advance method here
        if self.parent.progress_bar:
            if finish:
                self.parent.progressBarFinished()
            else:
                self.parent.progress_bar.advance()

    def cancel(self):
        """
        Cancel the current task (if any).
        """
        if self._task is not None:
            self._task.cancel()
            assert self._task.future.done()
            # disconnect the `_task_finished` slot
            self._task.watcher.done.disconnect(self._end_task)
            self._task = None

    def score_genes(self, **kwargs):
        """ Run gene enrichment.

        :param design:
        :param data_x:
        :param rows_by_cluster:
        :param method:


        Note:
            We do not apply filter nor notify view that data is changed. This is done after filters
        """
        if self._task is not None:
            # First make sure any pending tasks are cancelled.
            self.cancel()
        assert self._task is None

        progress_advance = methodinvoke(self, "progress_advance", (bool,))

        def callback():
            if self._task.cancelled:
                raise KeyboardInterrupt()
            progress_advance(self._task.cancelled)

        self.parent.progress_bar = ProgressBar(self.parent, iterations=len(self.get_rows()))
        f = partial(self._score_genes, callback=callback, **kwargs)
        self._task = Task()
        self._task.future = self._executor.submit(f)

        self._task.watcher = FutureWatcher(self._task.future)
        self._task.watcher.done.connect(self._end_task)

    def gene_sets_enrichment(self, gs_object, gene_sets, reference_genes):
        """ Run gene sets enrichment.

        :param gs_object:
        :param gene_sets:
        :param reference_genes:

        Note:
            We do not apply filter nor notify view that data is changed. This is done after filters

        """

        for item in self.get_rows():
            genes = [gene.gene_id for gene in item.filtered_genes]
            item.gene_set_enrichment(gs_object, gene_sets, set(genes), reference_genes)

    def apply_gene_filters(self, p_val=None, fdr=None, count=None):
        [item.filter_enriched_genes(p_val, fdr, max_gene_count=count) for item in self.get_rows()]
        self.dataChanged.emit(self.createIndex(0, 0), self.createIndex(self.rowCount(0), 0))

    def apply_gene_sets_filters(self, p_val=None, fdr=None, count=None):
        [item.filter_gene_sets(p_val, fdr, max_set_count=count) for item in self.get_rows()]
        self.dataChanged.emit(self.createIndex(0, 0), self.createIndex(self.rowCount(0), 0))
