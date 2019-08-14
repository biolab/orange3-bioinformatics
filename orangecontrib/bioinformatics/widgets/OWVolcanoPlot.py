import numpy as np

from Orange.widgets import gui, settings
from Orange.widgets.widget import Msg
from Orange.widgets.settings import SettingProvider
from Orange.widgets.visualize.owscatterplot import OWScatterPlotBase, OWDataProjectionWidget

from orangecontrib.bioinformatics.utils.statistics import score_t_test, score_fold_change
from orangecontrib.bioinformatics.widgets.utils.gui import label_selection
from orangecontrib.bioinformatics.widgets.utils.data import GENE_ID_COLUMN, GENE_AS_ATTRIBUTE_NAME


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

        gene_enrichment = Msg('{}, {}.')
        no_selected_gene_sets = Msg('No gene set selected, select them from Gene Sets box.')

    class Error(OWDataProjectionWidget.Error):
        exclude_error = Msg('Target labels most exclude/include at least one value.')
        negative_values = Msg('Negative values in the input. The inputs cannot be in ratio scale.')
        data_not_annotated = Msg('The input date is not annotated as expexted. Please refer to documentation.')
        gene_column_id_missing = Msg('Can not identify genes column. Please refer to documentation.')

    GRAPH_CLASS = VolcanoGraph
    graph = SettingProvider(VolcanoGraph)
    embedding_variables_names = ('log2 (ratio)', '-log10 (P_value)')

    stored_selections = settings.ContextSetting([])
    current_group_index = settings.ContextSetting(0)

    def __init__(self):
        super().__init__()

    def _add_controls(self):
        box = gui.vBox(self.controlArea, "Target Labels")
        self.group_selection_widget = label_selection.LabelSelectionWidget()
        self.group_selection_widget.groupChanged.connect(self.on_target_values_changed)
        self.group_selection_widget.groupSelectionChanged.connect(self.on_target_values_changed)
        box.layout().addWidget(self.group_selection_widget)

        super()._add_controls()
        self.gui.add_widgets([self.gui.ShowGridLines], self._plot_box)

    def get_embedding(self):
        self.Error.exclude_error.clear()

        group, target_indices = self.group_selection_widget.selected_split()

        if self.data and group is not None and target_indices:
            X = self.data.X
            I1 = label_selection.group_selection_mask(self.data, group, target_indices)
            I2 = ~I1

            # print(group)
            if isinstance(group, label_selection.RowGroup):
                X = X.T

            N1, N2 = np.count_nonzero(I1), np.count_nonzero(I2)

            if not N1 or not N2:
                self.Error.exclude_error()
                return

            if N1 < 2 and N2 < 2:
                self.Warning.insufficient_data()

            X1, X2 = X[:, I1], X[:, I2]

            if np.any(X1 < 0.0) or np.any(X2 < 0):
                self.Error.negative_values()
                X1 = np.full_like(X1, np.nan)
                X2 = np.full_like(X2, np.nan)

            with np.errstate(divide="ignore", invalid="ignore"):
                fold = score_fold_change(X1, X2, axis=1, log=True)
                _, p_values = score_t_test(X1, X2, axis=1)
                log_p_values = np.log10(p_values)

            self.valid_data = np.isfinite(fold) & np.isfinite(p_values)
            return np.array([fold, -log_p_values]).T

    def setup_plot(self):
        super().setup_plot()
        for axis, var in (("bottom", 'log<sub>2</sub> (ratio)'), ("left", '-log<sub>10</sub> (P_value)')):
            self.graph.set_axis_title(axis, var)

    def on_target_values_changed(self, index):
        # Save the current selection to persistent settings
        self.current_group_index = index
        selected_indices = [ind.row() for ind in self.group_selection_widget.currentGroupSelection().indexes()]

        if self.current_group_index != -1 and selected_indices:
            self.stored_selections[self.current_group_index] = selected_indices

        self.setup_plot()

    def set_data(self, data):
        self.Warning.clear()
        self.Error.clear()
        super().set_data(data)
        self.group_selection_widget.set_data(self, self.data)

        if self.data:
            if not self.stored_selections:
                self.stored_selections = [[0] for _ in self.group_selection_widget.targets]
            self.group_selection_widget.set_selection()

    def check_data(self):
        self.clear_messages()
        use_attr_names = self.data.attributes.get(GENE_AS_ATTRIBUTE_NAME, None)
        gene_id_column = self.data.attributes.get(GENE_ID_COLUMN, None)

        if self.data is not None and (len(self.data) == 0 or len(self.data.domain) == 0):
            self.data = None

        if use_attr_names is None:
            # Note: input data is not annotated properly.
            self.Error.data_not_annotated()
            self.data = None

        if gene_id_column is None:
            # Note: Can not identify genes column.
            self.Error.gene_column_id_missing()
            self.data = None


if __name__ == "__main__":
    pass
