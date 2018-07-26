""" Cluster analysis module """
import numpy as np
import threading
import concurrent.futures

from typing import Union, Optional
from operator import attrgetter
from functools import partial


from AnyQt.QtCore import (
    Qt, QAbstractListModel, QVariant, QThreadPool, Slot, QThread
)
from Orange.widgets.utils.concurrent import ThreadExecutor, FutureWatcher, methodinvoke
from Orange.widgets.gui import ProgressBar

from orangecontrib.bioinformatics.geneset import GeneSet
from orangecontrib.bioinformatics.widgets.utils.gui import gene_scoring_method
from orangecontrib.bioinformatics.utils.statistics import FDR
from orangecontrib.bioinformatics.ncbi.gene import Gene

GENE_COUNT = 20
GENE_SETS_COUNT = 5


class ClusterGene(Gene):
    __slots__ = ['score', 'p_val', 'fdr']

    def __init__(self, gene_name, gene_id):
        super().__init__()
        self.input_name = gene_name
        self.ncbi_id = gene_id

        # default values
        self.fdr = 1
        self.p_val = 1


class ClusterGeneSet(GeneSet):
    __slots__ = ['p_val', 'fdr', 'count']

    def __init__(self):
        super().__init__()


class Cluster:

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

        # misc
        self.use_gene_count = True
        self.delta = 0

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

    def filter_enriched_genes(self, count, p_val, fdr):
        filter_function = partial(self.apply_filter, p_val=p_val, fdr=fdr)
        filtered_list = list(sorted(filter(filter_function, self.genes), key=attrgetter('fdr')))
        self.delta = len(filtered_list)
        self.filtered_genes = filtered_list[:count]

    def filter_gene_sets(self, count, p_val, fdr):
        filter_function = partial(self.apply_filter, p_val=p_val, fdr=fdr)
        self.filtered_gene_sets = list(sorted(filter(filter_function, self.gene_sets), key=attrgetter('fdr'))[:count])

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
            self.gene_sets.append(gs)

        # calculate FDR
        p_vals = [gs.p_val for gs in self.gene_sets]
        fdr_values = FDR(p_vals)

        for gs, fdr in zip(self.gene_sets, fdr_values):
            gs.fdr = fdr

    def __update_gene_objects(self, p_vals, fdr_vals):
        # type: (Union[np.ndarray, list], Union[np.ndarray, list]) ->  None
        """ update gene objects with computed results

        :param p_vals:   Computed p-values
        :param fdr_vals: Computed fdr-values

        """
        for p_val, fdr, gene in zip(p_vals, fdr_vals, self.genes):
            gene.p_val = p_val
            gene.fdr = fdr

    def cluster_vs_rest(self, table_x, rows_by_cluster, method):
        # type: (np.ndarray, np.ndarray, gene_scoring_method) -> None

        cluster = table_x[rows_by_cluster == self.index]
        rest = table_x[rows_by_cluster != self.index]
        if cluster.any():
            scores, p_values = method.score_function(cluster, rest)
            fdr_values = FDR(p_values)
            self.__update_gene_objects(p_values, fdr_values)

    def cluster_vs_cluster(self, table_x, rows_by_cluster, method, **kwargs):
        # type: (np.ndarray, np.ndarray, gene_scoring_method) -> None

        aggregation = kwargs.get('aggregation', 'max')

        calculated_p_values = []
        for cluster_index in set(rows_by_cluster):

            # ignore current cluster index
            if cluster_index != self.index:
                cluster = table_x[rows_by_cluster == self.index]
                rest = table_x[rows_by_cluster == cluster_index]

                _, p_values = method.score_function(cluster, rest)
                calculated_p_values.append(p_values)

        if aggregation == 'max':
            max_p_values = np.max(np.array(calculated_p_values), axis=0)
            fdr_values = FDR(max_p_values)

            self.__update_gene_objects(max_p_values, fdr_values)

    def to_html(self):
        gene_sets = '(no enriched gene sets)'
        if self.filtered_gene_sets:
            gene_sets = '<br>'.join(['<b>{}</b> (FDR={:0.2e}, n={})'.format(g_set.name, g_set.fdr, g_set.count)
                                    for g_set in self.filtered_gene_sets])

        genes = '(all genes are filtered out)'
        if self.filtered_genes:
            genes = ', '.join([gene.input_name for gene in self.filtered_genes])
            self.delta = GENE_COUNT if self.delta > GENE_COUNT else self.delta

            if (self.delta - len(self.filtered_genes)) > 0 and self.use_gene_count:
                genes += ', ... ({} more genes)'.format(self.delta - len(self.filtered_genes))

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
        </table>
        </body>
        </html>
        """.format(self.name, gene_sets, genes)
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
        self.parent.progress_bar.finish()

        try:
            f.result()
        except Exception as ex:
            print(ex)

    def _score_genes(self, design, data_x, rows_by_cluster, method, callback):

        for item in self.get_rows():
            if design:
                item.cluster_vs_cluster(data_x, rows_by_cluster, method)
            else:
                item.cluster_vs_rest(data_x, rows_by_cluster, method)

            callback()

    @Slot()
    def progress_advance(self):
        # GUI should be updated in main thread. That's why we are calling advance method here
        if self.parent.progress_bar:
            self.parent.progress_bar.advance()

    def cancel(self):
        """
        Cancel the current task (if any).
        """
        if self._task is not None:
            self._task.cancel()
            assert self._task.future.done()
            # disconnect the `_task_finished` slot
            self._task.watcher.done.disconnect(self._score_genes)
            self._task = None

    def score_genes(self, design, data_x, rows_by_cluster, method):
        """ Run gene enrichment.

        :param design:
        :param data_x:
        :param rows_by_cluster:
        :param method:


        Note:
            We do not apply filter nor notify view that data is changed. This is done after filters
        """

        self._task = Task()
        progress_advance = methodinvoke(self, "progress_advance")

        def callback():
            if self._task.cancelled:
                raise KeyboardInterrupt()
            if self.parent.progress_bar:
                progress_advance()

        self.parent.progress_bar = ProgressBar(self.parent, iterations=len(self.get_rows()))
        f = partial(self._score_genes, design, data_x, rows_by_cluster, method, callback=callback)
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
            genes = [gene.ncbi_id for gene in item.filtered_genes]
            item.gene_set_enrichment(gs_object,
                                     gene_sets,
                                     set(genes),
                                     reference_genes)

    def apply_gene_filters(self, count=GENE_COUNT, p_val=None, fdr=None):
        [item.filter_enriched_genes(count, p_val, fdr) for item in self.get_rows()]
        self.dataChanged.emit(self.createIndex(0, 0), self.createIndex(self.rowCount(0), 0))

    def apply_gene_sets_filters(self, count=GENE_SETS_COUNT, p_val=None, fdr=None):
        [item.filter_gene_sets(count, p_val, fdr) for item in self.get_rows()]
        self.dataChanged.emit(self.createIndex(0, 0), self.createIndex(self.rowCount(0), 0))




