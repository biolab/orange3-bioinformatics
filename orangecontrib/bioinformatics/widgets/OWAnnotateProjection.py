# pylint: disable=too-many-ancestors
from enum import IntEnum
from types import SimpleNamespace
from typing import Dict, Tuple, Optional
from itertools import chain

import numpy as np
import pyqtgraph as pg

from AnyQt.QtGui import QColor
from AnyQt.QtCore import Qt, QRectF, QObject

from Orange.data import Table, Domain, DiscreteVariable, ContinuousVariable
from Orange.widgets import gui, report
from Orange.data.util import array_equal
from Orange.data.filter import Values, FilterString
from Orange.widgets.widget import Msg, Input
from Orange.widgets.settings import Setting, ContextSetting, SettingProvider
from Orange.widgets.utils.sql import check_sql_input
from Orange.widgets.utils.concurrent import TaskState, ConcurrentWidgetMixin
from Orange.widgets.utils.itemmodels import DomainModel
from Orange.widgets.utils.colorpalette import ColorPaletteGenerator
from Orange.widgets.utils.widgetpreview import WidgetPreview
from Orange.widgets.visualize.utils.widget import OWDataProjectionWidget
from Orange.widgets.visualize.owscatterplotgraph import OWScatterPlotBase

from orangecontrib.bioinformatics.widgets.utils.data import TAX_ID
from orangecontrib.bioinformatics.annotation.annotate_samples import (
    PFUN_BINOMIAL,
    SCORING_LOG_FDR,
    SCORING_EXP_RATIO,
    PFUN_HYPERGEOMETRIC,
    SCORING_MARKERS_SUM,
    AnnotateSamplesMeta,
)
from orangecontrib.bioinformatics.annotation.annotate_projection import annotate_projection, cluster_additional_points

CELL_TYPE = "Cell Type"
ENTREZ_ID = "Entrez ID"


class Result(SimpleNamespace):
    scores = None  # type: Optional[Table]
    clusters = None  # type: OWAnnotateProjection.Clusters


class Runner:
    @staticmethod
    def compute_scores(
        data: Table,
        genes: Table,
        p_threshold: float,
        p_value_fun: str,
        scoring: str,
        start: float,
        end: float,
        result: Result,
        state: TaskState,
    ):
        if not data or not genes:
            result.scores.z_vals = None
            result.scores.annotations = None
            result.scores.p_vals = None
            result.scores.table = None
        else:
            state.set_status("Computing scores...")
            weights = np.array([15, 75, 10]) * (end - start) / 100

            if not result.scores.z_vals:
                result.scores.z_vals = AnnotateSamplesMeta.mann_whitney_test(data)
                state.set_partial_result(("scores", result))
            state.set_progress_value(weights[0])
            if state.is_interruption_requested():
                return

            if not result.scores.annotations or not result.scores.p_vals:
                annot, p_vals = AnnotateSamplesMeta.assign_annotations(
                    result.scores.z_vals, genes, data, p_value_fun=p_value_fun, scoring=scoring
                )
                result.scores.annotations = annot
                result.scores.p_vals = p_vals
                state.set_partial_result(("scores", result))
            state.set_progress_value(weights[1])
            if state.is_interruption_requested():
                return

            result.scores.table = AnnotateSamplesMeta.filter_annotations(
                result.scores.annotations, result.scores.p_vals, p_threshold=p_threshold
            )

        state.set_partial_result(("scores", result))

    @staticmethod
    def compute_clusters(embedding: Table, result: Result, state: TaskState):
        if not result.scores.table or not embedding:
            result.clusters.table = None
            result.clusters.groups = None
        else:
            state.set_status("Finding clusters...")
            kwargs = {}
            if result.clusters.epsilon is not None:
                kwargs["eps"] = result.clusters.epsilon
            clusters = annotate_projection(result.scores.table, embedding, **kwargs)
            result.clusters.table = clusters[0]
            result.clusters.groups = clusters[1]
            result.clusters.epsilon = clusters[2]
        state.set_partial_result(("clusters", result))

    @staticmethod
    def compute_secondary_clusters(embedding: Table, result: Result, state: TaskState):
        if not result.clusters.groups or not embedding:
            result.clusters.secondary_table = None
        else:
            state.set_status("Finding secondary clusters...")
            hulls = {k: v[2] for k, v in result.clusters.groups.items()}
            clusters = result.clusters.table
            domain = clusters and clusters.domain["Clusters"]
            table = cluster_additional_points(embedding, hulls, domain)
            result.clusters.secondary_table = table
        state.set_partial_result(("secondary_clusters", result))

    @classmethod
    def run(
        cls,
        data: Table,
        secondary_data: Table,
        attr_x: ContinuousVariable,
        attr_y: ContinuousVariable,
        genes: Table,
        p_threshold: float,
        p_value_fun: str,
        scoring: str,
        result: Result,
        state: TaskState,
    ):

        start, step, weights = 0, 0, np.array([80, 15, 5])

        def set_progress():
            nonlocal start
            nonlocal step
            start = int(start + weights[step])
            step += 1
            state.set_progress_value(start)
            return 0 if state.is_interruption_requested() else 1

        if not result.scores.table:
            end = start + weights[step]
            cls.compute_scores(data, genes, p_threshold, p_value_fun, scoring, start, end, result, state)
        if not set_progress():
            return result

        if not result.clusters.table:
            embedding = data.transform(Domain([attr_x, attr_y]))
            cls.compute_clusters(embedding, result, state)
        if not set_progress():
            return result

        if not result.clusters.secondary_table and secondary_data:
            embedding = secondary_data.transform(Domain([attr_x, attr_y]))
            cls.compute_secondary_clusters(embedding, result, state)
        set_progress()
        return result


class CenteredTextItem(pg.TextItem):
    def __init__(self, view_box, x, y, text, tooltip):
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
    show_ref_data = Setting(True)

    def __init__(self, scatter_widget, parent):
        super().__init__(scatter_widget, parent)
        self.cluster_hulls_items = []
        self.cluster_labels_items = []
        self.ref_scatterplot_item = None
        self._tooltip_delegate = EventDelegate()  # remove points tooltip

    def clear(self):
        super().clear()
        self.cluster_hulls_items.clear()
        self.cluster_labels_items.clear()
        self.ref_scatterplot_item = None

    def reset_view(self):
        x, y = [self.get_coordinates()[0]], [self.get_coordinates()[1]]
        if x[0] is None or y[0] is None:
            return

        hulls = self.master.get_cluster_hulls()
        if hulls is not None:
            x.extend([hull[:, 0] for hull, _ in hulls])
            y.extend([hull[:, 1] for hull, _ in hulls])

        points = self.master.get_coordinates_reference_data()
        if points is not None:
            x.append(points[0])
            y.append(points[1])

        x, y = np.hstack(x), np.hstack(y)
        min_x, max_x, min_y, max_y = np.min(x), np.max(x), np.min(y), np.max(y)
        rect = QRectF(min_x, min_y, max_x - min_x or 1, max_y - min_y or 1)
        self.view_box.setRange(rect, padding=0.025)

    def update_coordinates(self):
        self.update_reference_coordinates()
        super().update_coordinates()
        self.update_clusters()
        self.view_box.setAspectLocked(True, 1)
        self.update_reference_item()
        self.reset_view()

    def update_clusters(self):
        self._update_cluster_hull()
        self._update_cluster_labels()

    def _update_cluster_hull(self):
        for item in self.cluster_hulls_items:
            self.plot_widget.removeItem(item)
        if not self.show_cluster_hull:
            return
        hulls = self.master.get_cluster_hulls()
        if hulls is None:
            return
        for hull, color in hulls:
            pen = pg.mkPen(color=QColor(*color), style=Qt.DashLine, width=3)
            item = pg.PlotCurveItem(x=hull[:, 0], y=hull[:, 1], pen=pen, antialias=True)
            self.plot_widget.addItem(item)
            self.cluster_hulls_items.append(item)

    def _update_cluster_labels(self):
        for item in self.cluster_labels_items:
            self.plot_widget.removeItem(item)
        if not self.n_cluster_labels:
            return
        labels = self.master.get_cluster_labels()
        if labels is None:
            return
        for label_per, (x, y), _ in labels:
            text = "\n".join([l for l, _ in label_per[: self.n_cluster_labels]])
            ttip = "\n".join([f"{round(p * 100)}%  {l}" for l, p in label_per])
            item = CenteredTextItem(self.view_box, x, y, text, ttip)
            self.plot_widget.addItem(item)
            self.cluster_labels_items.append(item)

    def update_reference_coordinates(self):
        points = self.master.get_coordinates_reference_data()
        if points is None:
            return
        if self.ref_scatterplot_item is None:
            color = pg.mkColor(200, 200, 200)
            pen, brush = pg.mkPen(color=color), pg.mkBrush(color=color)
            size = OWScatterPlotBase.MinShapeSize + 3
            self.ref_scatterplot_item = pg.ScatterPlotItem(x=points[0], y=points[1], pen=pen, brush=brush, size=size)
            self.plot_widget.addItem(self.ref_scatterplot_item)
        else:
            self.ref_scatterplot_item.setData(x=points[0], y=points[1])

    def update_reference_item(self):
        if self.ref_scatterplot_item is not None:
            self.ref_scatterplot_item.setVisible(self.show_ref_data)
        elif self.scatterplot_item is not None:
            self.scatterplot_item.setVisible(self.show_ref_data)


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
    priority = 140
    keywords = ["annotate"]

    GRAPH_CLASS = OWAnnotateProjectionGraph
    graph = SettingProvider(OWAnnotateProjectionGraph)
    # remove this line when https://github.com/biolab/orange3/pull/3863 is released
    left_side_scrolling = True

    attr_x = ContextSetting(None)
    attr_y = ContextSetting(None)
    scoring_method = Setting(ScoringMethod.ExpRatio)
    statistical_test = Setting(StatisticalTest.Binomial)
    p_threshold = Setting(0.05)
    use_user_epsilon = Setting(False)
    epsilon = Setting(0)
    color_by_cluster = Setting(False)

    class Scores(SimpleNamespace):
        z_vals = None  # type: Optional[Table]
        annotations = None  # type: Optional[Table]
        p_vals = None  # type: Optional[Table]
        table = None  # type: Optional[Table]

    class Clusters(SimpleNamespace):
        table = None  # type: Optional[Table]
        groups = None  # type: Optional[Dict[str, Tuple]]
        epsilon = None  # type: Optional[float]
        secondary_table = None  # type: Optional[Table]

    class Inputs:
        data = Input("Reference Data", Table, default=True)
        secondary_data = Input("Secondary Data", Table)
        genes = Input("Genes", Table)

    class Information(OWDataProjectionWidget.Information):
        modified = Msg("The parameter settings have been changed. Press " "\"Start\" to rerun with the new settings.")

    class Warning(OWDataProjectionWidget.Warning):
        no_genes = Msg("Missing Genes table on input.")
        missing_compute_value = Msg(
            "{} or {} is missing in Secondary data and" " cannot be computed from Reference data."
        )
        missing_tax_id = Msg(f"'{TAX_ID}' is missing in Reference data. " f"Try using 'Genes' widget.")
        missing_entrez_id = Msg("'Entred ID' is missing in Genes table.")
        missing_cell_type = Msg("'Cell Type' is missing in Genes table.")
        missing_tax_id_genes = Msg(f"'{TAX_ID}' is missing in Genes table. " f"Try using 'Marker Genes' widget.")
        different_tax_id = Msg(f"Data and Genes appear to describe different " f"organisms (mismatching {TAX_ID}).")
        same_axis_features = Msg(f"Selected features for Axis x and " f"Axis y should differ.")

    class Error(OWDataProjectionWidget.Error):
        no_reference_data = Msg("Missing Reference data on input.")
        no_continuous_vars = Msg("Data has no continuous variables.")
        not_enough_inst = Msg("Not enough instances in Reference data.")
        proj_error = Msg("An error occurred while annotating data.\n{}")
        sec_proj_error = Msg("An error occurred while projecting " "Secondary data.\n{}")

    def __init__(self):
        OWDataProjectionWidget.__init__(self)
        ConcurrentWidgetMixin.__init__(self)
        # inputs
        self.reference_data = None  # type: Optional[Table]
        self.secondary_data = None  # type: Optional[Table]
        self.genes = None  # type: Optional[Table]
        # annotations
        self.scores = None  # type: OWAnnotateProjection.Scores
        self.clusters = None  # type: OWAnnotateProjection.Clusters
        self.__invalidate_scores()
        self.__invalidate_clusters()

    # GUI
    def _add_controls(self):
        self.__add_annotation_controls()
        super()._add_controls()
        self.gui.add_control(
            self._effects_box,
            gui.hSlider,
            "Cluster labels:",
            master=self.graph,
            value="n_cluster_labels",
            minValue=0,
            maxValue=3,
            step=1,
            createLabel=False,
            callback=self.graph.update_clusters,
        )
        self._plot_box.children()[1].hide()  # Hide 'Show color regions'
        gui.checkBox(
            self._plot_box, self.graph, "show_cluster_hull", "Show cluster hull", callback=self.graph.update_clusters
        )
        gui.checkBox(
            self._plot_box,
            self,
            "color_by_cluster",
            "Color points by cluster",
            callback=self.__color_by_cluster_changed,
        )
        gui.checkBox(
            self._plot_box, self.graph, "show_ref_data", "Show reference data", callback=self.__show_ref_data_changed
        )

    def __add_annotation_controls(self):
        common_options = {
            'labelWidth': 100,
            'orientation': Qt.Horizontal,
            'sendSelectedValue': True,
            'contentsLength': 14,
        }
        box = gui.vBox(self.controlArea, True)
        ord = (DomainModel.METAS, DomainModel.ATTRIBUTES, DomainModel.CLASSES)
        mod = DomainModel(ord, valid_types=ContinuousVariable)
        gui.comboBox(
            box, self, "attr_x", label="Axis x:", model=mod, callback=self.__axis_attr_changed, **common_options
        )
        gui.comboBox(
            box, self, "attr_y", label="Axis y:", model=mod, callback=self.__axis_attr_changed, **common_options
        )
        gui.comboBox(
            box,
            self,
            "scoring_method",
            label="Scoring method:",
            items=ScoringMethod.items(),
            orientation=Qt.Horizontal,
            contentsLength=13,
            labelWidth=100,
            callback=self.__scoring_combo_changed,
        )
        gui.comboBox(
            box,
            self,
            "statistical_test",
            label="Statistical test:",
            items=StatisticalTest.items(),
            orientation=Qt.Horizontal,
            labelWidth=100,
            callback=self.__scoring_combo_changed,
        )
        gui.doubleSpin(
            box, self, "p_threshold", 0, 1, 0.01, label="FDR threshold:", callback=self.__p_threshold_changed
        )
        hbox = gui.hBox(box)
        gui.checkBox(hbox, self, "use_user_epsilon", "ε for DBSCAN:", callback=self.__epsilon_check_changed)
        self.epsilon_spin = gui.doubleSpin(hbox, self, "epsilon", 0, 10, 0.1, callback=self.__epsilon_changed)
        self.run_button = gui.button(box, self, "Start", self._toggle_run)

    @property
    def effective_variables(self):
        return self.reference_data.domain.attributes

    @property
    def effective_data(self):
        return self.reference_data.transform(
            Domain(self.effective_variables, self.reference_data.domain.class_vars, self.reference_data.domain.metas)
        )

    @property
    def can_annotate(self):
        return (
            self.reference_data
            and self.genes
            and TAX_ID in self.reference_data.attributes
            and ENTREZ_ID in self.genes.domain
            and CELL_TYPE in self.genes.domain
            and TAX_ID in self.genes.attributes
            and self.reference_data.attributes[TAX_ID] == self.genes.attributes[TAX_ID]
            and self.attr_x is not None
            and self.attr_y is not None
            and self.attr_x is not self.attr_y
        )

    def __color_by_cluster_changed(self):
        self.controls.attr_color.setEnabled(not self.color_by_cluster)
        self.graph.update_colors()

    def __show_ref_data_changed(self):
        self.graph.update_reference_item()

    def __axis_attr_changed(self):
        self.__parameter_changed()
        self.setup_plot()

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
        self.scores = self.Scores(z_vals=None, annotations=None, p_vals=None, scores=None)

    def __invalidate_scores_table(self):
        self.scores.table = None

    def __invalidate_scores_annotations(self):
        self.scores.annotations = None
        self.scores.p_vals = None
        self.__invalidate_scores_table()

    def __invalidate_clusters(self):
        epsilon = self.epsilon if self.use_user_epsilon else None
        self.clusters = self.Clusters(table=None, groups=None, epsilon=epsilon, secondary_table=None)

    def __invalidate_clusters_secondary_table(self):
        self.clusters.secondary_table = None

    def _toggle_run(self):
        if self.task is not None:
            self.cancel()
            self.run_button.setText("Resume")
            self.commit()
        else:
            self._run()

    def _run(self):
        self.Information.modified.clear()
        self.graph.update_clusters()  # Remove cluster hulls and labels
        if not self.can_annotate:
            return

        self.run_button.setText("Stop")
        result = Result(scores=self.scores, clusters=self.clusters)
        self.start(
            Runner.run,
            self.reference_data,
            self.secondary_data,
            self.attr_x,
            self.attr_y,
            self.genes,
            self.p_threshold,
            StatisticalTest.values()[self.statistical_test],
            ScoringMethod.values()[self.scoring_method],
            result,
        )

    def on_partial_result(self, value: Tuple[str, Result]):
        which, result = value
        if which == "scores":
            self.scores = result.scores
        elif which == "clusters":
            self.clusters = result.clusters
            if result.clusters.epsilon is not None:
                self.epsilon = result.clusters.epsilon
            self.graph.update_clusters()
            self.graph.update_colors()
            self.graph.reset_view()
        elif which == "secondary_clusters":
            self.clusters.secondary_table = result.clusters.secondary_table
            self.graph.update_colors()

    def on_done(self, result: Result):
        self.scores = result.scores
        self.clusters = result.clusters
        if result.clusters.epsilon is not None:
            self.epsilon = result.clusters.epsilon
        self.run_button.setText("Start")
        self.commit()

    def on_exception(self, ex: Exception):
        self.Error.proj_error(ex)
        self.run_button.setText("Start")

    @Inputs.data
    @check_sql_input
    def set_data(self, data):
        attr_x, attr_y = self.attr_x, self.attr_y
        data_existed = self.reference_data is not None
        effective_data = self.effective_data if data_existed else None
        super().set_data(data)
        self.reference_data = self.data
        if not (
            data_existed and self.reference_data is not None and array_equal(effective_data.X, self.effective_data.X)
        ):
            self.clear()
            self.__invalidate_scores()
            self.__invalidate_clusters()
        elif attr_x is not self.attr_x or attr_y is not self.attr_y:
            self.clear()
            self.__invalidate_clusters()

    def check_data(self):
        self.Warning.missing_tax_id.clear()
        self.Error.no_continuous_vars.clear()
        self.Error.not_enough_inst.clear()
        if self.data:
            if len(self.data) < 2:
                self.Error.not_enough_inst()
                self.data = None
            elif not self.data.domain.has_continuous_attributes(True, True):
                self.Error.no_continuous_vars()
                self.data = None
            elif TAX_ID not in self.data.attributes:
                self.Warning.missing_tax_id()

    def init_attr_values(self):
        super().init_attr_values()
        domain = self.data.domain if self.data and len(self.data) else None
        model = self.controls.attr_x.model()
        model.set_domain(domain)
        self.attr_x = model[0] if model else None
        self.attr_y = model[1] if len(model) >= 2 else self.attr_x

    def enable_controls(self):
        super().enable_controls()
        self.enable_epsilon_spin()
        self.controls.attr_color.setEnabled(not self.color_by_cluster)

    def enable_epsilon_spin(self):
        self.epsilon_spin.setEnabled(self.use_user_epsilon)

    @Inputs.secondary_data
    @check_sql_input
    def set_secondary_data(self, data):
        self.secondary_data = data
        self.clear()
        self.__invalidate_clusters_secondary_table()

    @Inputs.genes
    def set_genes(self, genes: Optional[Table]):
        self.genes = genes
        self.__invalidate_scores_annotations()
        self.__invalidate_clusters()
        self.check_genes()
        self.cancel()

    def check_genes(self):
        self.Warning.missing_tax_id_genes.clear()
        self.Warning.missing_entrez_id.clear()
        self.Warning.missing_cell_type.clear()
        if self.genes:
            if ENTREZ_ID not in self.genes.domain:
                self.Warning.missing_entrez_id()
            elif CELL_TYPE not in self.genes.domain:
                self.Warning.missing_cell_type()
            elif TAX_ID not in self.genes.attributes:
                self.Warning.missing_tax_id_genes()

    def handleNewSignals(self):
        self.check_inputs()
        super().handleNewSignals()
        self._run()

    def check_inputs(self):
        self.Warning.no_genes.clear()
        self.Warning.different_tax_id.clear()
        self.Error.no_reference_data.clear()
        self.Error.proj_error.clear()
        if self.reference_data and not self.genes:
            self.Warning.no_genes()
        elif (
            self.reference_data
            and self.genes
            and TAX_ID in self.reference_data.attributes
            and TAX_ID in self.genes.attributes
            and (self.genes.attributes[TAX_ID] != self.reference_data.attributes[TAX_ID])
        ):
            self.Warning.different_tax_id()
        elif self.secondary_data and not self.reference_data:
            self.Error.no_reference_data()

        self.data = self.reference_data and self.secondary_data or self.reference_data
        self.init_attr_values_secondary()

    def init_attr_values_secondary(self):
        saved_attrs = {
            "attr_color": self.attr_color,
            "attr_shape": self.attr_shape,
            "attr_size": self.attr_size,
            "attr_label": self.attr_label,
        }
        domain = self.data.domain if self.data and len(self.data) else None
        for name, attr in saved_attrs.items():
            model = getattr(self.controls, name).model()
            model.set_domain(domain)
            if domain is not None and attr and attr.name in domain:
                setattr(self, name, domain[attr.name])
            else:
                setattr(self, name, None)

    def get_coordinates_reference_data(self):
        if not self.secondary_data or not self.reference_data:
            return None
        x_data = self.reference_data.get_column_view(self.attr_x)[0]
        y_data = self.reference_data.get_column_view(self.attr_y)[0]
        return np.vstack((x_data, y_data))

    def get_embedding(self):
        def no_embedding_attrs():
            domain = self.secondary_data.domain
            return (
                self.attr_x.name not in domain
                and self.attr_x.compute_value is None
                or self.attr_y.name not in domain
                and self.attr_y.compute_value is None
            )

        self.Warning.missing_compute_value.clear()
        self.Warning.same_axis_features.clear()
        self.Error.sec_proj_error.clear()
        if not self.data:
            self.valid_data = None
            return None

        if self.attr_x is self.attr_y and self.attr_x is not None:
            self.Warning.same_axis_features()

        if self.secondary_data and self.reference_data:
            if no_embedding_attrs():
                self.Warning.missing_compute_value(self.attr_x, self.attr_y)
            try:
                domain = self.reference_data.domain
                self.data = self.secondary_data.transform(domain)
            except Exception as ex:
                self.Error.sec_proj_error(ex)
                return None

        x_data = self.get_column(self.attr_x, filter_valid=False)
        y_data = self.get_column(self.attr_y, filter_valid=False)
        if x_data is None or y_data is None:
            return None

        self.valid_data = np.isfinite(x_data) & np.isfinite(y_data)
        return np.vstack((x_data, y_data)).T

    def get_color_data(self):
        table = self.clusters.secondary_table if self.secondary_data else self.clusters.table
        if not self.color_by_cluster or not table:
            return super().get_color_data()

        attr = table.domain["Clusters"]
        all_data = table.get_column_view(attr)[0]
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
        # the second option is to keep widget backward compatible (Orange < 3.25)
        return (
            self.clusters.table.domain["Clusters"].palette
            if hasattr(self.clusters.table.domain["Clusters"], "palette")
            else ColorPaletteGenerator(number_of_colors=len(colors), rgb_colors=colors)
        )

    def get_cluster_hulls(self):
        if not self.clusters.groups:
            return None
        var = self.clusters.table.domain["Clusters"]
        return [(hull, var.colors[var.values.index(key)]) for key, (_, _, hull) in self.clusters.groups.items()]

    def get_cluster_labels(self):
        if not self.clusters.groups:
            return None
        var = self.clusters.table.domain["Clusters"]
        return [(ann, lab, var.colors[var.values.index(key)]) for key, (ann, lab, _) in self.clusters.groups.items()]

    def send_data(self):
        if self.clusters.secondary_table:
            group_sel, graph = None, self.graph
            ref_proj_data, sec_proj_data = self._get_projection_data()
            data = Table.concatenate([ref_proj_data, sec_proj_data], axis=0)
            data = self.__add_data_source_column(data)

            if graph.selection is not None:
                ref_group_sel = np.zeros(len(ref_proj_data), dtype=int)
                sec_group_sel = np.zeros(len(sec_proj_data), dtype=int)
                sec_group_sel[self.valid_data] = graph.selection
                group_sel = np.hstack([ref_group_sel, sec_group_sel])

            selection = graph.get_selection() + len(ref_proj_data)
            self.Outputs.selected_data.send(self._get_selected_data(data, selection, group_sel))

            # keeping backward compatible - versions of Orange < 3.25 fails on first call
            try:
                self.Outputs.annotated_data.send(self._get_annotated_data(data, group_sel, graph.selection))
            except TypeError:
                self.Outputs.annotated_data.send(self._get_annotated_data(data, selection, group_sel, graph.selection))
        else:
            super().send_data()

    def _get_projection_data(self):
        if not self.scores.table or not self.clusters.table:
            return self.data

        def get_augmented_data(table, scores_x, clusters_x):
            new_table = table.transform(domain)
            if scores_x is not None:
                new_table.metas[:, :n_sco] = scores_x
            new_table.metas[:, n_sco:n_clu] = clusters_x
            return new_table

        n_sco = self.scores.table.X.shape[1]
        n_clu = n_sco + self.clusters.table.X.shape[1]
        domain = self.reference_data.domain
        metas = chain(self.scores.table.domain.attributes, self.clusters.table.domain.attributes, domain.metas)
        domain = Domain(domain.attributes, domain.class_vars, metas)
        if self.clusters.secondary_table:
            domain = self.__concatenate_domains(domain, self.secondary_data.domain)

        augmented_data = get_augmented_data(self.reference_data, self.scores.table.X, self.clusters.table.X)
        if self.clusters.secondary_table:
            sec_aug_data = get_augmented_data(self.secondary_data, None, self.clusters.secondary_table.X)
            return augmented_data, sec_aug_data
        return augmented_data

    @staticmethod
    def __concatenate_domains(ref_domain, sec_domain):
        attributes = list(ref_domain.attributes)
        class_vars = list(ref_domain.class_vars)
        metas = list(ref_domain.metas)
        for attr in sec_domain.attributes:
            if attr.name not in ref_domain:
                attributes.append(attr)
        for attr in sec_domain.class_vars:
            if attr.name not in ref_domain:
                class_vars.append(attr)
        for attr in sec_domain.metas:
            if attr.name not in ref_domain:
                metas.append(attr)
        return Domain(attributes, class_vars, metas)

    def __add_data_source_column(self, data):
        source_attr = DiscreteVariable("Data Source", values=[self.reference_data.name, self.secondary_data.name])
        domain = data.domain
        metas = chain(domain.metas, [source_attr])
        domain = Domain(domain.attributes, domain.class_vars, metas)
        data = data.transform(domain)
        data.metas[: len(self.reference_data), -1] = 0
        data.metas[len(self.reference_data) :, -1] = 1
        return data

    def _get_send_report_caption(self):
        color_vr_name = "Clusters" if self.color_by_cluster else self._get_caption_var_name(self.attr_color)
        return report.render_items_vert(
            (
                ("Color", color_vr_name),
                ("Label", self._get_caption_var_name(self.attr_label)),
                ("Shape", self._get_caption_var_name(self.attr_shape)),
                ("Size", self._get_caption_var_name(self.attr_size)),
                ("Jittering", self.graph.jitter_size != 0 and "{} %".format(self.graph.jitter_size)),
            )
        )

    def clear(self):
        super().clear()
        self.cancel()
        self._invalidated = True

    def onDeleteWidget(self):
        self.shutdown()
        super().onDeleteWidget()


if __name__ == "__main__":
    from Orange.projection import PCA
    from orangecontrib.bioinformatics.utils import serverfiles

    data_path = "https://datasets.biolab.si/sc/aml-1k.tab.gz"
    table_data = Table(data_path)
    table_data.attributes[TAX_ID] = "9606"

    ref_data = table_data[::2]
    pca = PCA(n_components=2)
    pca_model = pca(ref_data)
    proj = pca_model(ref_data)
    new_dom = Domain(
        ref_data.domain.attributes, ref_data.domain.class_vars, chain(ref_data.domain.metas, proj.domain.attributes)
    )
    ref_data = ref_data.transform(new_dom)

    genes_path = serverfiles.localpath_download("marker_genes", "panglao_gene_markers.tab")
    filter_ = FilterString("Organism", FilterString.Equal, "Human")
    table_genes = Values([filter_])(Table(genes_path))
    table_genes.attributes[TAX_ID] = "9606"

    WidgetPreview(OWAnnotateProjection).run(
        set_data=ref_data, set_secondary_data=table_data[1:200:2], set_genes=table_genes
    )
