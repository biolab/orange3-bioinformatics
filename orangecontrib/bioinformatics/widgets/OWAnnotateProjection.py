# pylint: disable=too-many-ancestors
from enum import IntEnum
from typing import Tuple, Optional
from types import SimpleNamespace as namespace
from itertools import chain
import numpy as np

from AnyQt.QtCore import Qt, QRectF, QObject
from AnyQt.QtGui import QColor

import pyqtgraph as pg

from Orange.data import Table, Domain
from Orange.data.filter import FilterString, Values
from Orange.projection import Projector, PCA
from Orange.widgets import gui, report
from Orange.widgets.settings import Setting, SettingProvider
from Orange.widgets.unsupervised.owtsne import pca_preprocessing, \
    prepare_tsne_obj
from Orange.widgets.utils.colorpalette import ColorPaletteGenerator
from Orange.widgets.utils.concurrent import TaskState, ConcurrentWidgetMixin
from Orange.widgets.utils.widgetpreview import WidgetPreview
from Orange.widgets.visualize.owscatterplotgraph import OWScatterPlotBase
from Orange.widgets.visualize.utils.widget import OWDataProjectionWidget
from Orange.widgets.widget import Input, Msg

from orangecontrib.bioinformatics.annotation.annotate_projection import \
    annotate_projection
from orangecontrib.bioinformatics.annotation.annotate_samples import \
    AnnotateSamples, SCORING_EXP_RATIO, SCORING_MARKERS_SUM, \
    SCORING_LOG_FDR, PFUN_BINOMIAL, PFUN_HYPERGEOMETRIC
from orangecontrib.bioinformatics.widgets.utils.data import TAX_ID


class Result(namespace):
    scores = None  # type: Optional[Table]
    embedding = None  # type: Optional[Table]
    clusters = None  # type: Clusters


class Runner:
    @staticmethod
    def compute_scores(data: Table, genes: Table, p_threshold: float,
                       p_value_fun: str, scoring: str, start: float,
                       end: float, result: Result, state: TaskState):
        if not data or not genes:
            result.scores.z_vals = None
            result.scores.annotations = None
            result.scores.p_vals = None
            result.scores.table = None
        else:
            state.set_status("Computing scores...")
            weights = np.array([15, 75, 10]) * (end - start) / 100

            if not result.scores.z_vals:
                result.scores.z_vals = AnnotateSamples.mann_whitney_test(data)
                state.set_partial_result(("scores", result))
            state.set_progress_value(weights[0])
            if state.is_interruption_requested():
                return

            if not result.scores.annotations or not result.scores.p_vals:
                annot, p_vals = AnnotateSamples.assign_annotations(
                    result.scores.z_vals, genes, data,
                    p_value_fun=p_value_fun, scoring=scoring)
                result.scores.annotations = annot
                result.scores.p_vals = p_vals
                state.set_partial_result(("scores", result))
            state.set_progress_value(weights[1])
            if state.is_interruption_requested():
                return

            result.scores.table = AnnotateSamples.filter_annotations(
                result.scores.annotations, result.scores.p_vals,
                p_threshold=p_threshold)

        state.set_partial_result(("scores", result))

    @staticmethod
    def compute_tsne_embedding(data: Table, start: float, end: float,
                               state: TaskState):
        count, weights = (0, np.array([7, 1, 1, 5, 1, 23, 61, 1]) *
                          (end - start) / 100)
        tsne = prepare_tsne_obj(data, 30, False, 1)

        def set_progress():
            nonlocal start
            nonlocal count
            start += weights[count]
            count += 1
            state.set_progress_value(start)
            return 0 if state.is_interruption_requested() else 1

        pca_data = pca_preprocessing(data, 20, True)
        if not set_progress():
            return None

        pp_data = tsne.preprocess(pca_data)
        if not set_progress():
            return None

        init = tsne.compute_initialization(pp_data.X)
        if not set_progress():
            return None

        affinities = tsne.compute_affinities(pp_data.X)
        if not set_progress():
            return None

        tsne_embedding = tsne.prepare_embedding(affinities, init)
        if not set_progress():
            return None

        def tsne_optimization(n_iter, **kwargs):
            nonlocal tsne_embedding
            step_size = 25
            n_iter_done = 0
            while n_iter_done < n_iter:
                n_iter_done += step_size
                tsne_embedding = tsne_embedding.optimize(
                    n_iter=step_size, inplace=True, **kwargs)
                step = n_iter_done / n_iter * weights[count]
                state.set_progress_value(start + step)
                if state.is_interruption_requested():
                    return None

        kwargs = {"exaggeration": tsne.early_exaggeration, "momentum": 0.5}
        tsne_optimization(tsne.early_exaggeration_iter, **kwargs)
        if not set_progress():
            return None

        kwargs = {"exaggeration": tsne.exaggeration, "momentum": 0.8}
        tsne_optimization(tsne.n_iter, **kwargs)
        if not set_progress():
            return None

        tsne_embedding = tsne.convert_embedding_to_model(
            pp_data, tsne_embedding)
        set_progress()
        return tsne_embedding.embedding

    @staticmethod
    def compute_embedding(data: Table, projector: Optional[Projector],
                          start: float, end: float,
                          result: Result, state: TaskState):
        if not data:
            result.embedding = None
        else:
            state.set_status("Computing 2D projection...")
            if projector is None:
                # Use t-SNE by default
                embedding = Runner.compute_tsne_embedding(
                    data, start, end, state)
            elif isinstance(projector, PCA):
                embedding = projector(data)(data)[:, :2]
            else:
                raise NotImplementedError
            result.embedding = embedding
        state.set_partial_result(("embedding", result))

    @staticmethod
    def compute_clusters(result: Result, state: TaskState):
        if not result.scores.table or not result.embedding:
            result.clusters.table = None
            result.clusters.groups = None
        else:
            state.set_status("Finding clusters...")
            kwargs = {}
            if result.clusters.epsilon is not None:
                kwargs["eps"] = result.clusters.epsilon
            clusters = annotate_projection(
                result.scores.table, result.embedding, **kwargs)
            result.clusters.table = clusters[0]
            result.clusters.groups = clusters[1]
            result.clusters.epsilon = clusters[2]
        state.set_partial_result(("clusters", result))

    @classmethod
    def run(cls, data: Table, genes: Table, projector: Optional[Projector],
            p_threshold: float, p_value_fun: str, scoring: str,
            result: Result, state: TaskState):

        start, step, weights = 0, 0, np.array([25, 70, 5])

        def set_progress():
            nonlocal start
            nonlocal step
            start = int(start + weights[step])
            step += 1
            state.set_progress_value(start)
            return 0 if state.is_interruption_requested() else 1

        if not result.scores.table:
            end = start + weights[step]
            cls.compute_scores(data, genes, p_threshold, p_value_fun,
                               scoring, start, end, result, state)
        if not set_progress():
            return result

        if not result.embedding:
            end = start + weights[step]
            cls.compute_embedding(data, projector, start, end, result, state)
        if not set_progress():
            return result

        if not result.clusters.table:
            cls.compute_clusters(result, state)
        set_progress()
        return result


class CenteredTextItem(pg.TextItem):
    def __init__(self, view_box, x, y, text, tooltip, rgb):
        bg_color = QColor(*rgb)  # TODO - remove rgb, or use it
        bg_color = QColor(Qt.white)
        bg_color.setAlpha(200)
        color = QColor(Qt.black)
        super().__init__(text, pg.mkColor(color), fill=pg.mkBrush(bg_color))
        self._x = x
        self._y = y
        self._view_box = view_box
        self._view_box.sigStateChanged.connect(self.center)
        self.textItem.setToolTip(tooltip)
        self.setPos(x, y)
        self.center()

    def center(self):
        br = self.boundingRect()
        dx = br.width() / 2 * self._view_box.viewPixelSize()[0]
        dy = br.height() / 2 * self._view_box.viewPixelSize()[1]
        self.setPos(self._x - dx, self._y + dy)


class EventDelegate(QObject):
    def eventFilter(self, *_):
        return False


class OWAnnotateProjectionGraph(OWScatterPlotBase):
    show_cluster_hull = Setting(True)
    n_cluster_labels = Setting(1)

    def __init__(self, scatter_widget, parent):
        super().__init__(scatter_widget, parent)
        self.cluster_hulls_items = []
        self.cluster_labels_items = []
        self._tooltip_delegate = EventDelegate()  # remove points tooltip

    def clear(self):
        super().clear()
        self.cluster_hulls_items.clear()
        self.cluster_labels_items.clear()

    def reset_view(self):
        x, y = self.get_coordinates()
        hulls = self.master.get_cluster_hulls()
        if not hulls or x is None or y is None:
            return
        x = np.hstack([hull[:, 0] for hull, _ in hulls] + [x])
        min_x, max_x = np.min(x), np.max(x)
        y = np.hstack([hull[:, 1] for hull, _ in hulls] + [y])
        min_y, max_y = np.min(y), np.max(y)
        rect = QRectF(min_x, min_y, max_x - min_x or 1, max_y - min_y or 1)
        self.view_box.setRange(rect, padding=0.025)

    def update_coordinates(self):
        super().update_coordinates()
        if self.scatterplot_item is not None:
            self.update_cluster_hull()
            self.update_cluster_labels()
            self.view_box.setAspectLocked(True, 1)

    def update_cluster_hull(self):
        for item in self.cluster_hulls_items:
            self.plot_widget.removeItem(item)
        if not self.show_cluster_hull:
            return
        hulls = self.master.get_cluster_hulls()
        if hulls is None:
            return
        for hull, color in hulls:
            pen = pg.mkPen(color=QColor(*color), style=Qt.DashLine, width=3)
            item = pg.PlotCurveItem(
                x=hull[:, 0], y=hull[:, 1], pen=pen, antialias=True)
            self.plot_widget.addItem(item)
            self.cluster_hulls_items.append(item)

    def update_cluster_labels(self):
        for item in self.cluster_labels_items:
            self.plot_widget.removeItem(item)
        if not self.n_cluster_labels:
            return
        labels = self.master.get_cluster_labels()
        if labels is None:
            return
        for label_per, (xp, yp), color in labels:
            text = "\n".join([l for l, _ in label_per[:self.n_cluster_labels]])
            ttip = "\n".join([f"{round(p * 100)}%  {l}" for l, p in label_per])
            item = CenteredTextItem(self.view_box, xp, yp, text, ttip, color)
            self.plot_widget.addItem(item)
            self.cluster_labels_items.append(item)


class ScoringMethod(IntEnum):
    ExpRatio, MarkersSum, LogFdr = range(3)

    @staticmethod
    def values():
        return [SCORING_LOG_FDR, SCORING_MARKERS_SUM, SCORING_EXP_RATIO]

    @staticmethod
    def items():
        return ["-Log(FDR)", "Marker expression", "Marker expression (%)"]


class StatisticalTest(IntEnum):
    Binomial, Hypergeometric = range(2)

    @staticmethod
    def values():
        return [PFUN_BINOMIAL, PFUN_HYPERGEOMETRIC]

    @staticmethod
    def items():
        return ["Binomial", "Hypergeometric"]


class OWAnnotateProjection(OWDataProjectionWidget, ConcurrentWidgetMixin):
    name = "Annotator"
    description = "Annotates projection clusters."
    icon = "icons/OWAnnotateProjection.svg"
    priority = 3050
    keywords = ["annotate", "projection", "tsne", "pca"]

    GRAPH_CLASS = OWAnnotateProjectionGraph
    graph = SettingProvider(OWAnnotateProjectionGraph)
    left_side_scrolling = True  # remove this line when https://github.com/biolab/orange3/pull/3863 is merged

    scoring_method = Setting(ScoringMethod.ExpRatio)
    statistical_test = Setting(StatisticalTest.Binomial)
    p_threshold = Setting(0.05)
    use_user_epsilon = Setting(False)
    epsilon = Setting(0)
    color_by_cluster = Setting(False)

    class Scores(namespace):
        z_vals = None  # type: Optional[Table]
        annotations = None  # type: Optional[Table]
        p_vals = None  # type: Optional[Table]
        table = None  # type: Optional[Table]

    class Clusters(namespace):
        table = None  # type: Optional[Table]
        groups = None  # type: Optional[Dict[str, Tuple]]
        epsilon = None  # type: Optional[float]

    class Inputs(OWDataProjectionWidget.Inputs):
        genes = Input("Genes", Table)
        projector = Input("Projector", Projector)

    class Information(OWDataProjectionWidget.Information):
        modified = Msg("The parameter settings have been changed. Press "
                       "\"Start\" to rerun with the new settings.")

    class Warning(OWDataProjectionWidget.Warning):
        no_genes = Msg("Missing genes table on input.")

    class Error(OWDataProjectionWidget.Error):
        proj_error = Msg("An error occurred while annotating data.\n{}")
        missing_tax_id = Msg(f"'{TAX_ID}' is missing in data table. "
                             f"Try using 'Genes' widget.")
        missing_entrez_id = Msg("'Entred ID' is missing in genes table.")
        missing_cell_type = Msg("'Cell Type' is missing in genes table.")

    def __init__(self):
        OWDataProjectionWidget.__init__(self)
        ConcurrentWidgetMixin.__init__(self)
        # inputs
        self.projector = None  # type: Optional[Projector]
        self.genes = None  # type: Optional[Table]
        # annotations
        self.scores = None  # type: Scores
        self.embedding = None  # type: Optional[Table]
        self.clusters = None  # type: Clusters
        self.__invalidate_scores()
        self.__invalidate_clusters()

    # GUI
    def _add_controls(self):
        self.__add_annotation_controls()
        super()._add_controls()
        self.gui.add_control(
            self._effects_box, gui.hSlider, "Cluster labels:",
            master=self.graph, value="n_cluster_labels", minValue=0,
            maxValue=3, step=1,
            createLabel=False, callback=self.graph.update_cluster_labels
        )
        self._plot_box.children()[1].hide()  # Hide 'Show color regions'
        gui.checkBox(
            self._plot_box, self.graph, "show_cluster_hull",
            "Show cluster hull", callback=self.graph.update_cluster_hull)
        gui.checkBox(
            self._plot_box, self, "color_by_cluster",
            "Color points by cluster",
            callback=self.__color_by_cluster_changed)

    def __add_annotation_controls(self):
        box = gui.vBox(self.controlArea, True)
        gui.comboBox(
            box, self, "scoring_method", label="Scoring method:",
            items=ScoringMethod.items(), orientation=Qt.Horizontal,
            contentsLength=13, labelWidth=100,
            callback=self.__scoring_combo_changed)
        gui.comboBox(
            box, self, "statistical_test", label="Statistical test:",
            items=StatisticalTest.items(), orientation=Qt.Horizontal,
            labelWidth=100, callback=self.__scoring_combo_changed)
        gui.doubleSpin(
            box, self, "p_threshold", 0, 1, 0.01, label="FDR threshold:",
            callback=self.__p_threshold_changed)
        hbox = gui.hBox(box)
        gui.checkBox(
            hbox, self, "use_user_epsilon", "ε for DBSCAN:",
            callback=self.__epsilon_check_changed)
        self.epsilon_spin = gui.doubleSpin(
            hbox, self, "epsilon", 0, 10, 0.1, callback=self.__epsilon_changed)
        self.run_button = gui.button(box, self, "Start", self._toggle_run)

    def __color_by_cluster_changed(self):
        self.controls.attr_color.setEnabled(not self.color_by_cluster)
        self.graph.update_colors()

    def __scoring_combo_changed(self):
        self.__invalidate_scores_annotations()
        self.__parameter_changed()

    def __p_threshold_changed(self):
        self.__invalidate_scores_table()
        self.__parameter_changed()

    def __epsilon_changed(self):
        self.__parameter_changed()

    def __parameter_changed(self):
        self.__invalidate_clusters()
        self.Information.modified()
        self.Error.proj_error.clear()

    def __epsilon_check_changed(self):
        self.enable_epsilon_spin()
        if not self.use_user_epsilon:
            self.__epsilon_changed()

    def __invalidate_scores(self):
        self.scores = self.Scores(z_vals=None, annotations=None,
                                  p_vals=None, scores=None)

    def __invalidate_scores_table(self):
        self.scores.table = None

    def __invalidate_scores_annotations(self):
        self.scores.annotations = None
        self.scores.p_vals = None
        self.__invalidate_scores_table()

    def __invalidate_clusters(self):
        epsilon = self.epsilon if self.use_user_epsilon else None
        self.clusters = self.Clusters(table=None, groups=None, epsilon=epsilon)

    def _toggle_run(self):
        if self.task is not None:
            self.cancel()
            self.run_button.setText("Resume")
            self.commit()
        else:
            self._run()

    def _run(self):
        self.Information.modified.clear()
        if not self.data:
            return

        # Remove cluster hulls and labels
        self.graph.update_cluster_hull()
        self.graph.update_cluster_labels()

        self.run_button.setText("Stop")
        result = Result(scores=self.scores,
                        embedding=self.embedding,
                        clusters=self.clusters)
        self.start(Runner.run, self.data, self.genes, self.projector,
                   self.p_threshold,
                   StatisticalTest.values()[self.statistical_test],
                   ScoringMethod.values()[self.scoring_method], result)

    def on_partial_result(self, value: Tuple[str, Result]):
        which, result = value
        if which == "scores":
            self.scores = result.scores
        elif which == "embedding":
            self.embedding = result.embedding
            self.setup_plot()
        elif which == "clusters":
            self.clusters = result.clusters
            if result.clusters.epsilon is not None:
                self.epsilon = result.clusters.epsilon
            self.graph.update_cluster_hull()
            self.graph.update_cluster_labels()
            self.graph.update_colors()
            self.graph.reset_view()

    def on_done(self, result: Result):
        self.scores = result.scores
        self.embedding = result.embedding
        self.clusters = result.clusters
        if result.clusters.epsilon is not None:
            self.epsilon = result.clusters.epsilon
        self.run_button.setText("Start")
        self.commit()

    def on_exception(self, ex: Exception):
        self.Error.proj_error(ex)
        self.run_button.setText("Start")

    def check_data(self):
        super().check_data()
        if self.data:
            if TAX_ID not in self.data.attributes:
                self.Error.missing_tax_id()
                self.data = None
        if self.data and not self.genes:
            self.Warning.no_genes()

    def check_genes(self):
        if self.genes:
            metas = [m.name for m in self.genes.domain.metas]
            if "Entrez ID" not in metas:
                self.Error.missing_entrez_id()
                self.genes = None
            elif "Cell Type" not in metas:
                self.Error.missing_cell_type()
                self.genes = None

    def enable_controls(self):
        super().enable_controls()
        self.enable_epsilon_spin()
        self.controls.attr_color.setEnabled(not self.color_by_cluster)

    def enable_epsilon_spin(self):
        self.epsilon_spin.setEnabled(self.use_user_epsilon)

    @Inputs.genes
    def set_genes(self, genes: Optional[Table]):
        self.genes = genes
        self.__invalidate_scores_annotations()
        self.__invalidate_clusters()
        self.check_data()
        self.check_genes()

    @Inputs.projector
    def set_projector(self, projector: Optional[Projector]):
        self.cancel()
        self.embedding = None
        self._invalidated = True
        self.__invalidate_clusters()
        self.projector = projector

    def handleNewSignals(self):
        self._run()
        super().handleNewSignals()

    def get_coordinates_data(self):
        coordinates = super().get_coordinates_data()
        if coordinates[0] is not None and len(coordinates) == 1:
            return np.vstack((coordinates, np.zeros_like(coordinates)))
        return coordinates

    def get_embedding(self):
        if not self.embedding:
            self.valid_data = None
            return None

        self.valid_data = np.all(np.isfinite(self.embedding.X), 1)
        return self.embedding.X

    def get_color_data(self):
        if not self.color_by_cluster or not self.clusters.table:
            return super().get_color_data()

        attr = self.clusters.table.domain["Clusters"]
        all_data = self.clusters.table.get_column_view(attr)[0]
        if all_data.dtype == object and attr.is_primitive():
            all_data = all_data.astype(float)
        if self.valid_data is not None:
            all_data = all_data[self.valid_data]
        return all_data

    def get_color_labels(self):
        if not self.color_by_cluster or not self.clusters.table:
            return super().get_color_labels()
        return self.clusters.table.domain["Clusters"].values

    def is_continuous_color(self):
        if not self.color_by_cluster or not self.clusters.table:
            return super().is_continuous_color()
        return False

    def get_palette(self):
        if not self.color_by_cluster or not self.clusters.table:
            return super().get_palette()

        colors = self.clusters.table.domain["Clusters"].colors
        return ColorPaletteGenerator(number_of_colors=len(colors),
                                     rgb_colors=colors)

    def get_cluster_hulls(self):
        if not self.clusters.groups:
            return None
        var = self.clusters.table.domain["Clusters"]
        return [(hull, var.colors[var.values.index(key)])
                for key, (_, _, hull) in self.clusters.groups.items()]

    def get_cluster_labels(self):
        if not self.clusters.groups:
            return None
        var = self.clusters.table.domain["Clusters"]
        return [(ann, lab, var.colors[var.values.index(key)])
                for key, (ann, lab, _) in self.clusters.groups.items()]

    def _get_projection_data(self):
        if not self.scores.table or not self.embedding or not self.clusters.table:
            return self.data
        n_sco = self.scores.table.X.shape[1]
        n_emb = self.embedding.X.shape[1]
        n_clu = self.clusters.table.X.shape[1]
        domain = self.data.domain
        metas = chain(self.scores.table.domain.attributes,
                      self.embedding.domain.attributes,
                      self.clusters.table.domain.attributes,
                      domain.metas)
        domain = Domain(domain.attributes, domain.class_vars, metas)

        augmented_data = self.data.transform(domain)
        augmented_data.metas[:, :n_sco] = self.scores.table.X
        augmented_data.metas[:, n_sco:n_sco + n_emb] = self.embedding.X
        ind = n_sco + n_emb + n_clu
        augmented_data.metas[:, n_sco + n_emb: ind] = self.clusters.table.X
        augmented_data.metas[:, ind:] = self.data.metas
        return augmented_data

    def _get_send_report_caption(self):
        color_vr_name = "Clusters" if self.color_by_cluster else \
            self._get_caption_var_name(self.attr_color)
        return report.render_items_vert((
            ("Color", color_vr_name),
            ("Label", self._get_caption_var_name(self.attr_label)),
            ("Shape", self._get_caption_var_name(self.attr_shape)),
            ("Size", self._get_caption_var_name(self.attr_size)),
            ("Jittering", self.graph.jitter_size != 0 and
             "{} %".format(self.graph.jitter_size))))

    def clear(self):
        super().clear()
        self.cancel()
        self.embedding = None
        self.__invalidate_scores()
        self.__invalidate_clusters()

    def onDeleteWidget(self):
        self.shutdown()
        super().onDeleteWidget()


if __name__ == "__main__":
    from orangecontrib.bioinformatics.utils import serverfiles

    data_path = "https://datasets.orange.biolab.si/sc/aml-1k.tab.gz"
    table_data = Table(data_path)
    table_data.attributes[TAX_ID] = "9606"
    genes_path = serverfiles.localpath_download(
        "marker_genes", "panglao_gene_markers.tab")
    f = FilterString("Organism", FilterString.Equal, "Human")
    table_genes = Values([f])(Table(genes_path))

    WidgetPreview(OWAnnotateProjection).run(
        set_data=table_data,
        set_subset_data=table_data[::10],
        set_genes=table_genes,
        # set_projector=PCA()
    )
