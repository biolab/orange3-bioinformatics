from typing import Optional

import numpy as np

from Orange.data import Table
from Orange.widgets import gui, settings
from Orange.widgets.widget import Msg
from Orange.widgets.settings import ContextSetting, SettingProvider, DomainContextHandler, PerfectDomainContextHandler
from Orange.widgets.visualize.owscatterplot import OWScatterPlotBase, OWDataProjectionWidget

from orangecontrib.bioinformatics.utils.statistics import score_fold_change
from orangecontrib.bioinformatics.widgets.components import GeneScoringComponent
from orangecontrib.bioinformatics.widgets.utils.data import TableAnnotation
from orangecontrib.bioinformatics.widgets.utils.settings import GeneScoringComponentSettings


class VolcanoGraph(OWScatterPlotBase):
    label_only_selected = settings.Setting(True)

    def set_axis_title(self, axis, title):
        self.plot_widget.setLabel(axis=axis, text=title)


class OWVolcanoPlot(OWDataProjectionWidget):
    name = "Volcano Plot"
    description = "Plots fold change vs. p-value."
    icon = "icons/OWVolcanoPlot.svg"
    priority = 100

    class Warning(OWDataProjectionWidget.Warning):
        insufficient_data = Msg(
            'Insufficient data to compute statistics.' 'More than one measurement per class should be provided '
        )

    class Error(OWDataProjectionWidget.Error):
        exclude_error = Msg('Target labels most exclude/include at least one value.')
        negative_values = Msg('Negative values in the input. The inputs cannot be in ratio scale.')
        data_not_annotated = Msg('The input date is not annotated as expected. Please refer to documentation.')
        gene_column_id_missing = Msg('Can not identify genes column. Please refer to documentation.')

    settingsHandler = GeneScoringComponentSettings()

    # graph settings
    graph = SettingProvider(VolcanoGraph)
    GRAPH_CLASS = VolcanoGraph
    embedding_variables_names = ('log2 (ratio)', '-log10 (P_value)')

    # component settings
    scoring_component = SettingProvider(GeneScoringComponent)

    def __init__(self):
        super().__init__()
        self.original_data: Optional[Table] = None
        self.genes_in_columns: Optional[str] = None
        self.gene_id_column: Optional[str] = None
        self.gene_id_attribute: Optional[str] = None

        self.fold: Optional[np.array] = None
        self.log_p_values: Optional[np.array] = None
        self.valid_data: Optional[np.array] = None

    def _add_controls(self):
        box = gui.vBox(self.controlArea, True, margin=0)
        self.scoring_component = GeneScoringComponent(self, box)
        self.scoring_component.controls_changed.connect(self.setup_plot)

        super()._add_controls()
        self.gui.add_widgets([self.gui.ShowGridLines], self._plot_box)

    def _compute(self):
        self.Error.exclude_error.clear()
        self.fold = None
        self.log_p_values = None

        if self.data:
            x = self.data.X
            score_method = self.scoring_component.get_score_method()

            i1 = self.scoring_component.get_selection_mask()
            i2 = ~i1

            n1, n2 = np.count_nonzero(i1), np.count_nonzero(i2)
            if not n1 or not n2:
                self.Error.exclude_error()
                return

            if n1 < 2 and n2 < 2:
                self.Warning.insufficient_data()

            x1, x2 = x[:, i1], x[:, i2]
            if np.any(x1 < 0.0) or np.any(x2 < 0):
                self.Error.negative_values()
                x1 = np.full_like(x1, np.nan)
                x2 = np.full_like(x2, np.nan)

            with np.errstate(divide='ignore', invalid='ignore'):
                self.fold = score_fold_change(x1, x2, axis=1, log=True)
                _, p_values = score_method(x1, x2, axis=1, threshold=self.scoring_component.get_expression_threshold())
                self.log_p_values = np.log10(p_values)

    def get_embedding(self):
        if self.data is None:
            return None

        if self.fold is None or self.log_p_values is None:
            return

        self.valid_data = np.isfinite(self.fold) & np.isfinite(self.log_p_values)
        return np.array([self.fold, -self.log_p_values]).T

    def send_data(self):
        group_sel, data, graph = None, self._get_projection_data(), self.graph
        if graph.selection is not None:
            group_sel = np.zeros(len(data), dtype=int)
            group_sel[self.valid_data] = graph.selection

        selected_data = self._get_selected_data(data, graph.get_selection(), group_sel)

        if self.genes_in_columns and selected_data:
            selected_data = Table.transpose(selected_data, feature_names_column=self.scoring_component.feature_name)

        self.Outputs.selected_data.send(selected_data)

    def setup_plot(self):
        self._compute()
        super().setup_plot()
        for axis, var in (("bottom", 'log<sub>2</sub> (ratio)'), ("left", '-log<sub>10</sub> (P_value)')):
            self.graph.set_axis_title(axis, var)

    def set_data(self, data):
        self.Warning.clear()
        self.Error.clear()

        self.genes_in_columns = data.attributes.get(TableAnnotation.gene_as_attr_name, None)
        self.gene_id_column = data.attributes.get(TableAnnotation.gene_id_column, None)
        self.gene_id_attribute = data.attributes.get(TableAnnotation.gene_id_attribute, None)

        if self.genes_in_columns:
            self.original_data = data
            # override default meta_attr_name value to avoid unexpected changes.
            data = Table.transpose(data, meta_attr_name=self.scoring_component.feature_name)

        super().set_data(data)
        self.scoring_component.initialize(data)

    def check_data(self):
        self.clear_messages()
        if self.data is not None and (len(self.data) == 0 or len(self.data.domain) == 0):
            self.data = None

        if self.genes_in_columns is None:
            # Note: input data is not annotated properly.
            self.Error.data_not_annotated()
            self.data = None


if __name__ == "__main__":
    from Orange.widgets.utils.widgetpreview import WidgetPreview

    WidgetPreview(OWVolcanoPlot).run()
