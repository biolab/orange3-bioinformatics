# pylint: disable=too-many-ancestors
from typing import Tuple
from types import SimpleNamespace as namespace
from itertools import chain
import numpy as np

from AnyQt.QtCore import Qt, QRectF, QObject
from AnyQt.QtGui import QColor

import pyqtgraph as pg

from Orange.data import Table, Domain
from Orange.data.filter import FilterString, Values
from Orange.projection import Projector, TSNE
from Orange.widgets import gui
from Orange.widgets.settings import Setting, SettingProvider
from Orange.widgets.utils.concurrent import TaskState, ConcurrentWidgetMixin
from Orange.widgets.utils.widgetpreview import WidgetPreview
from Orange.widgets.visualize.owscatterplotgraph import OWScatterPlotBase
from Orange.widgets.visualize.utils.widget import OWDataProjectionWidget
from Orange.widgets.widget import Input, Msg

from orangecontrib.bioinformatics.annotation.annotate_projection import \
    annotate_projection
from orangecontrib.bioinformatics.annotation.annotate_samples import \
    AnnotateSamples
from orangecontrib.bioinformatics.widgets.utils.data import TAX_ID


class Result(namespace):
    scores = None  # type: Optional[Table]
    embedding = None  # type: Optional[Table]
    clusters = None  # type: Clusters


class Runner:
    @staticmethod
    def compute_scores(data: Table, genes: Table, p_value_th: float,
                       result: Result, state: TaskState):
        if not data or not genes:
            result.scores = None
        else:
            state.set_status("Computing scores...")
            result.scores = AnnotateSamples.annotate_samples(
                data, genes, p_threshold=p_value_th)
        state.set_partial_result(("scores", result))

    @staticmethod
    def compute_embedding(data: Table, projector: Projector,
                          result: Result, state: TaskState):
        if not data:
            result.embedding = None
        else:
            state.set_status("Computing 2D projection...")
            if isinstance(projector, TSNE):
                embedding = projector(data).embedding
            else:
                # Only OWPCA has Projector output
                embedding = projector(data)(data)[:, :2]
            result.embedding = embedding
        state.set_partial_result(("embedding", result))

    @staticmethod
    def compute_clusters(result: Result, state: TaskState):
        if not result.scores or not result.embedding:
            result.scores = None
        else:
            state.set_status("Finding clusters...")
            kwargs = {}
            if result.clusters.epsilon is not None:
                kwargs["eps"] = result.clusters.epsilon
            clusters = annotate_projection(
                result.scores, result.embedding, **kwargs)
            result.clusters.table = clusters[0]
            result.clusters.groups = clusters[1]
            result.clusters.epsilon = clusters[2]
        state.set_partial_result(("clusters", result))

    @classmethod
    def run(cls, data: Table, genes: Table, projector: Projector,
            p_value_th: float, result: Result, state: TaskState):

        state.set_progress_value(0)
        if state.is_interruption_requested():
            return result

        if result.scores is None:
            cls.compute_scores(data, genes, p_value_th, result, state)
        state.set_progress_value(20)
        if state.is_interruption_requested():
            return result

        if result.embedding is None:
            cls.compute_embedding(data, projector, result, state)
        state.set_progress_value(80)
        if state.is_interruption_requested():
            return result

        if result.clusters is None or result.clusters.table is None:
            cls.compute_clusters(result, state)
        state.set_progress_value(100)

        return result


class CenteredTextItem(pg.TextItem):
    def __init__(self, view_box, x, y, text, color, tooltip):
        super().__init__(text, color)
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
        hulls = self.master.get_cluster_hulls()
        if not hulls:
            return
        x = np.hstack([hull[:, 0] for hull in hulls])
        min_x, max_x = np.min(x), np.max(x)
        y = np.hstack([hull[:, 1] for hull in hulls])
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
        for hull in hulls:
            item = pg.PlotCurveItem(
                x=hull[:, 0], y=hull[:, 1],
                pen=pg.mkPen(QColor(Qt.black), width=1), antialias=True)
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
        for label_per, (xp, yp) in labels:
            text = "\n".join([l for l, _ in label_per[:self.n_cluster_labels]])
            ttip = "\n".join([f"{l}  {round(p * 100)}%" for l, p in label_per])
            item = CenteredTextItem(self.view_box, xp, yp, text,
                                    pg.mkColor(0, 0, 0), ttip)
            self.plot_widget.addItem(item)
            self.cluster_labels_items.append(item)


class OWAnnotateProjection(OWDataProjectionWidget, ConcurrentWidgetMixin):
    name = "Annotator"
    description = "Annotates projection clusters."
    icon = "icons/OWAnnotateProjection.svg"
    priority = 3050
    keywords = ["annotate", "projection", "tsne", "pca"]

    GRAPH_CLASS = OWAnnotateProjectionGraph
    graph = SettingProvider(OWAnnotateProjectionGraph)
    left_side_scrolling = True  # remove this line when https://github.com/biolab/orange3/pull/3863 is merged

    p_value_th = Setting(0.05)
    use_user_epsilon = Setting(False)
    epsilon = Setting(0)

    PROJECTOR = TSNE(n_components=2)

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
        missing_tax_id = Msg(f"'{TAX_ID}' is missing in data table.")
        missing_entrez_id = Msg("'Entred ID' is missing in genes table.")
        missing_cell_type = Msg("'Cell Type' is missing in genes table.")

    def __init__(self):
        OWDataProjectionWidget.__init__(self)
        ConcurrentWidgetMixin.__init__(self)
        # inputs
        self.projector = self.PROJECTOR  # type: Projector
        self.genes = None  # type: Optional[Table]
        # annotations
        self.scores = None  # type: Optional[Table]
        self.embedding = None  # type: Optional[Table]
        self.clusters = None  # type: Clusters
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
        gui.checkBox(
            self._plot_box, self.graph, "show_cluster_hull",
            "Show cluster hull", callback=self.graph.update_cluster_hull)

    def __add_annotation_controls(self):
        box = gui.vBox(self.controlArea, True)
        gui.doubleSpin(
            box, self, "p_value_th", 0, 1, 0.01, label="p-value threshold:",
            callback=self.__p_value_th_changed)
        hbox = gui.hBox(box)
        gui.checkBox(
            hbox, self, "use_user_epsilon", "ε for DBSCAN:",
            callback=self.__epsilon_check_changed)
        self.epsilon_spin = gui.doubleSpin(
            hbox, self, "epsilon", 0, 10, 0.1, callback=self.__epsilon_changed)
        self.run_button = gui.button(box, self, "Start", self._toggle_run)

    def __p_value_th_changed(self):
        self.scores = None
        self.__invalidate_clusters()
        self.Information.modified()

    def __epsilon_changed(self):
        self.__invalidate_clusters()
        self.Information.modified()

    def __epsilon_check_changed(self):
        self.enable_epsilon_spin()
        if not self.use_user_epsilon:
            self.__epsilon_changed()

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
        if self.data is None:
            return
        self.run_button.setText("Stop")
        result = Result(scores=self.scores,
                        embedding=self.embedding,
                        clusters=self.clusters)
        self.start(Runner.run, self.data, self.genes, self.projector,
                   self.p_value_th, result)

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

    def enable_epsilon_spin(self):
        self.epsilon_spin.setEnabled(self.use_user_epsilon)

    @Inputs.genes
    def set_genes(self, genes: Table):
        self.genes = genes
        self.scores = None
        self.__invalidate_clusters()
        self.check_data()
        self.check_genes()

    @Inputs.projector
    def set_projector(self, projector: Projector):
        self.clear()
        self._invalidated = True
        self.__invalidate_clusters()
        self.projector = projector or self.PROJECTOR

    def handleNewSignals(self):
        self._run()
        super().handleNewSignals()

    def get_coordinates_data(self):
        coordinates = super().get_coordinates_data()
        if coordinates[0] is not None and len(coordinates) == 1:
            return np.vstack((coordinates, np.zeros_like(coordinates)))
        return coordinates

    def get_embedding(self):
        if self.embedding is None:
            self.valid_data = None
            return None

        self.valid_data = np.all(np.isfinite(self.embedding.X), 1)
        return self.embedding.X

    def get_cluster_hulls(self):
        if self.clusters is None or self.clusters.groups is None:
            return None
        return [hull for _, _, hull in self.clusters.groups.values()]

    def get_cluster_labels(self):
        if self.clusters is None or self.clusters.groups is None:
            return None
        return [(ann, lab) for ann, lab, _ in self.clusters.groups.values()]

    def _get_projection_data(self):
        if not self.scores or not self.embedding or not self.clusters.table:
            return self.data
        n_sco = self.scores.X.shape[1]
        n_emb = self.embedding.X.shape[1]
        n_clu = self.clusters.table.X.shape[1]
        domain = self.data.domain
        metas = chain(self.scores.domain.attributes,
                      self.embedding.domain.attributes,
                      self.clusters.table.domain.attributes,
                      domain.metas)
        domain = Domain(domain.attributes, domain.class_vars, metas)

        augmented_data = self.data.transform(domain)
        augmented_data.metas[:, :n_sco] = self.scores.X
        augmented_data.metas[:, n_sco:n_sco + n_emb] = self.embedding.X
        ind = n_sco + n_emb + n_clu
        augmented_data.metas[:, n_sco + n_emb: ind] = self.clusters.table.X
        augmented_data.metas[:, ind:] = self.data.metas
        return augmented_data

    def clear(self):
        super().clear()
        self.cancel()
        self.embedding = None

    def onDeleteWidget(self):
        self.shutdown()
        super().onDeleteWidget()


if __name__ == "__main__":
    from Orange.projection import PCA
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
        set_projector=PCA()
    )
