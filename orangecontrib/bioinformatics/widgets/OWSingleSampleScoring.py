from typing import Any, Dict, Optional

import numpy as np
import pandas as pd
from single_sample_gsea import ss_gsea

from Orange.data import Table, Domain, ContinuousVariable
from Orange.widgets import gui
from Orange.widgets.widget import Msg, Input, Output, OWWidget
from Orange.widgets.settings import Setting
from Orange.widgets.utils.concurrent import TaskState, ConcurrentWidgetMixin
from Orange.widgets.utils.widgetpreview import WidgetPreview

from orangecontrib.bioinformatics.geneset import GeneSets
from orangecontrib.bioinformatics.widgets.utils.data import (
    TableAnnotation,
    check_table_annotation,
)


def ss_mean(df, gene_sets, callback=None):
    scores = {}

    for gs_name, gs_genes in gene_sets.items():
        genes = [gene for gene in gs_genes if gene in df.index]
        scores[gs_name] = df.loc[genes].mean()

        if callback is not None:
            callback()

    df_results = pd.DataFrame(scores)
    df_results.index = df.columns
    return df_results


scoring_methods = [('mean', ss_mean), ('ssGSEA', ss_gsea)]


def worker(data: Table, gene_sets, selected_method: int, state: TaskState):
    progress_steps = iter(np.linspace(0, 100, len(gene_sets)))

    def callback():
        state.set_progress_value(next(progress_steps))

    # Prepare data
    gene_id_attribute = data.attributes[TableAnnotation.gene_id_attribute]
    genes = [attr.attributes[gene_id_attribute] for attr in data.domain.attributes]
    data_x = np.transpose(data.X)
    df = pd.DataFrame(data_x)
    df.index = genes

    _, method = scoring_methods[selected_method]
    return method(df, gene_sets, callback=callback)


class OWSingleSampleScoring(OWWidget, ConcurrentWidgetMixin):
    name = 'Single sample scoring'
    description = 'Scoring gene sets by single sample.'
    icon = 'icons/OWSingleSampleScoring.svg'
    priority = 150
    keywords = ['single sample scoring', 'gene sets', 'ssGSEA']
    want_main_area = False

    scoring_method: int = Setting(1)
    auto_commit: bool = Setting(True, schema_only=True)

    class Inputs:
        data = Input('Data', Table)
        gene_sets = Input('Gene Sets', (Table, GeneSets))

    class Outputs:
        data = Output('Data', Table)

    class Warning(OWWidget.Warning):
        wrong_data_format = Msg(
            'Please transpose the data. The widget expects '
            'samples as rows and genes as columns.'
        )

    def __init__(self):
        OWWidget.__init__(self)
        ConcurrentWidgetMixin.__init__(self)

        self.data: Optional[Table] = None
        self.gene_sets: Optional[Dict[str, set]] = None

        box = gui.widgetBox(self.controlArea, 'Method', margin=0)
        gui.comboBox(
            box,
            self,
            'scoring_method',
            items=[label for label, _ in scoring_methods],
            callback=self.commit.deferred,
        )

        gui.rubber(self.controlArea)
        self.commit_button = gui.auto_commit(
            self.buttonsArea, self, 'auto_commit', '&Commit', box=False
        )

    @Inputs.data
    @check_table_annotation
    def set_data(self, data: Table):
        self.Warning.wrong_data_format.clear()

        if not data:
            return

        gene_as_attr_name = data.attributes[TableAnnotation.gene_as_attr_name]
        if not gene_as_attr_name:
            self.Warning.wrong_data_format()
            self.sample_column_box.hide()
            return

        self.data = data
        self.start_worker()

    @Inputs.gene_sets
    def set_gene_sets(self, gene_sets):
        self.gene_sets = None
        if not gene_sets:
            return
        self.gene_sets = {gs.name: gs.genes for gs in gene_sets}

        self.start_worker()

    @gui.deferred
    def commit(self):
        self.start_worker()

    def start_worker(self):
        if not self.data or not self.gene_sets:
            return

        self.start(worker, self.data, self.gene_sets, self.scoring_method)

    def on_done(self, worker_result: pd.DataFrame):
        table = None

        if worker_result is not None:
            scores = []
            for column in worker_result.columns:
                var = ContinuousVariable(column)
                var.attributes['Score'] = 'Score'
                scores.append(var)

            domain = Domain(
                self.data.domain.attributes,
                metas=self.data.domain.metas + tuple(scores),
                class_vars=self.data.domain.class_vars,
            )

            table = self.data.transform(domain)
            with table.unlocked():
                table[:, scores] = worker_result.values

        self.Outputs.data.send(table)

    def on_partial_result(self, result: Any) -> None:
        pass

    def on_exception(self, ex):
        raise ex


if __name__ == '__main__':
    previewer = WidgetPreview(OWSingleSampleScoring)
    previewer.run(Table('iris'))
