# Test methods with long descriptive names can omit docstrings
# pylint: disable=missing-docstring
import unittest
from unittest.mock import Mock

import numpy as np

from Orange.data import Table, Variable
from Orange.data.filter import FilterString, Values
from Orange.projection import PCA, MDS
from Orange.widgets.tests.base import WidgetTest, WidgetOutputsTestMixin, \
    ProjectionWidgetTestMixin, simulate

from orangecontrib.bioinformatics.utils import serverfiles
from orangecontrib.bioinformatics.widgets.OWAnnotateProjection import \
    OWAnnotateProjection
from orangecontrib.bioinformatics.widgets.utils.data import TAX_ID


class TestOWAnnotateProjection(WidgetTest, ProjectionWidgetTestMixin,
                               WidgetOutputsTestMixin):
    @classmethod
    def setUpClass(cls):
        super().setUpClass()
        path = "https://datasets.orange.biolab.si/sc/aml-1k.tab.gz"
        Variable._clear_all_caches()
        cls.data = Table(path)[::50]
        cls.data.attributes[TAX_ID] = "9606"

        cls.signal_name = "Data"
        cls.signal_data = cls.data
        cls.same_input_output_domain = False

        genes_path = serverfiles.localpath_download(
            "marker_genes", "panglao_gene_markers.tab")
        f = FilterString("Organism", FilterString.Equal, "Human")
        cls.genes = Values([f])(Table(genes_path))

    def setUp(self):
        self.widget = self.create_widget(OWAnnotateProjection)
        self.send_signal(self.widget.Inputs.projector, PCA())

    def tearDown(self):
        self.widget.cancel()
        self.wait_until_stop_blocking()

    def test_input_projector(self):
        self.send_signal(self.widget.Inputs.projector, None)
        self.send_signal(self.widget.Inputs.genes, self.genes)
        self.send_signal(self.widget.Inputs.data, self.data)
        self.wait_until_stop_blocking()
        self.assertFalse(self.widget.Error.proj_error.is_shown())
        self.assertIsNotNone(self.widget.embedding)
        self.send_signal(self.widget.Inputs.projector, PCA())
        self.assertIsInstance(self.widget.projector, PCA)
        self.wait_until_stop_blocking()
        self.send_signal(self.widget.Inputs.projector, MDS())
        self.wait_until_stop_blocking()
        self.assertTrue(self.widget.Error.proj_error.is_shown())

    def test_input_genes(self):
        self.send_signal(self.widget.Inputs.data, self.data)
        self.assertTrue(self.widget.Warning.no_genes.is_shown())
        self.wait_until_stop_blocking()
        self.send_signal(self.widget.Inputs.genes, self.genes)
        self.assertFalse(self.widget.Warning.no_genes.is_shown())
        self.wait_until_stop_blocking()
        self.send_signal(self.widget.Inputs.genes, None)
        self.assertTrue(self.widget.Warning.no_genes.is_shown())
        self.wait_until_stop_blocking()
        self.send_signal(self.widget.Inputs.data, None)
        self.assertFalse(self.widget.Warning.no_genes.is_shown())

    @unittest.skip("Skip due to timeout")
    def test_scoring_method_control(self):
        def start():
            self.widget.run_button.click()
            self.wait_until_stop_blocking()

        self.send_signal(self.widget.Inputs.genes, self.genes)
        self.send_signal(self.widget.Inputs.data, self.data)
        self.wait_until_stop_blocking()
        output1 = self.get_output(self.widget.Outputs.annotated_data)
        cbox = self.widget.controls.scoring_method
        simulate.combobox_activate_index(cbox, 2)
        start()
        output2 = self.get_output(self.widget.Outputs.annotated_data)
        np.testing.assert_array_equal(output1.X, output2.X)
        np.testing.assert_array_equal(output1.Y, output2.Y)

    @unittest.skip("Skip due to timeout")
    def test_statistical_test_control(self):
        def start():
            self.widget.run_button.click()
            self.wait_until_stop_blocking()

        self.send_signal(self.widget.Inputs.genes, self.genes)
        self.send_signal(self.widget.Inputs.data, self.data)
        self.wait_until_stop_blocking()
        output1 = self.get_output(self.widget.Outputs.annotated_data)
        cbox = self.widget.controls.statistical_test
        simulate.combobox_activate_index(cbox, 1)
        start()
        output2 = self.get_output(self.widget.Outputs.annotated_data)
        np.testing.assert_array_equal(output1.X, output2.X)
        np.testing.assert_array_equal(output1.Y, output2.Y)
        self.assertFalse((output1.metas == output2.metas).all())

    def test_p_threshold_control(self):
        self.send_signal(self.widget.Inputs.genes, self.genes)
        self.send_signal(self.widget.Inputs.data, self.data)
        self.wait_until_stop_blocking()
        output1 = self.get_output(self.widget.Outputs.annotated_data)
        self.widget.controls.p_threshold.valueChanged.emit(0.1)
        self.widget.run_button.click()
        self.wait_until_stop_blocking()
        output2 = self.get_output(self.widget.Outputs.annotated_data)
        np.testing.assert_array_equal(output1.X, output2.X)
        np.testing.assert_array_equal(output1.Y, output2.Y)

    def test_epsilon_control(self):
        self.send_signal(self.widget.Inputs.data, self.data)
        self.assertFalse(self.widget.epsilon_spin.isEnabled())
        self.widget.controls.use_user_epsilon.click()
        self.assertTrue(self.widget.epsilon_spin.isEnabled())
        self.assertFalse(self.widget.Information.modified.is_shown())
        self.widget.controls.use_user_epsilon.click()
        self.assertTrue(self.widget.Information.modified.is_shown())

    def test_output_data(self):
        self.send_signal(self.widget.Inputs.data, self.data)
        self.wait_until_stop_blocking()
        self.send_signal(self.widget.Inputs.genes, self.genes)
        self.wait_until_stop_blocking()
        output = self.get_output(self.widget.Outputs.annotated_data)
        n = len(self.data.domain.metas)
        self.assertGreater(len(output.domain.metas), n + 4)
        self.assertListEqual(
            ["PC1", "PC2", "Clusters", "Annotation"] +
            [m.name for m in self.data.domain.metas] + ["Selected"],
            [m.name for m in output.domain.metas[-n-5:]])
        np.testing.assert_array_equal(output.X, self.data.X)
        np.testing.assert_array_equal(output.Y, self.data.Y)
        np.testing.assert_array_equal(output.metas[:, -n - 1:-1],
                                      self.data.metas)

    def test_button_no_data(self):
        self.widget.run_button.click()
        self.assertEqual(self.widget.run_button.text(), "Start")

    def test_button_with_data(self):
        self.send_signal(self.widget.Inputs.data, self.data)
        self.assertEqual(self.widget.run_button.text(), "Stop")
        self.wait_until_stop_blocking()
        self.assertEqual(self.widget.run_button.text(), "Start")

    def test_button_toggle(self):
        self.send_signal(self.widget.Inputs.data, self.data)
        self.widget.run_button.click()
        self.assertEqual(self.widget.run_button.text(), "Resume")

    def test_plot_once(self):
        self.widget.setup_plot = Mock()
        self.widget.commit = Mock()
        self.send_signal(self.widget.Inputs.genes, self.genes)
        self.widget.commit.assert_called_once()
        self.widget.commit.reset_mock()
        self.send_signal(self.widget.Inputs.data, self.data)
        self.widget.setup_plot.assert_called_once()
        self.widget.commit.assert_called_once()
        self.wait_until_stop_blocking()
        self.widget.setup_plot.reset_mock()
        self.widget.commit.reset_mock()
        self.send_signal(self.widget.Inputs.data_subset, self.data[::10])
        self.widget.setup_plot.assert_not_called()
        self.widget.commit.assert_called_once()

    def test_saved_selection(self):
        self.send_signal(self.widget.Inputs.genes, self.genes)
        self.send_signal(self.widget.Inputs.data, self.data)
        self.wait_until_stop_blocking()

        self.widget.graph.select_by_indices(list(range(0, len(self.data), 10)))
        settings = self.widget.settingsHandler.pack_data(self.widget)
        w = self.create_widget(self.widget.__class__, stored_settings=settings)

        self.send_signal(w.Inputs.genes, self.genes)
        self.send_signal(w.Inputs.data, self.data, widget=w)
        self.wait_until_stop_blocking(widget=w)

        self.assertEqual(np.sum(w.graph.selection), 2)
        np.testing.assert_equal(self.widget.graph.selection, w.graph.selection)

    def test_outputs(self):
        self.send_signal(self.widget.Inputs.genes, self.genes)
        super().test_outputs()

    def test_color_by_cluster(self):
        self.send_signal(self.widget.Inputs.genes, self.genes)
        self.send_signal(self.widget.Inputs.data, self.data)
        self.wait_until_stop_blocking()
        self.assertTrue(self.widget.controls.attr_color.isEnabled())

        simulate.combobox_activate_index(self.widget.controls.attr_color, 0)
        is_grey = [brush.color().name() == "#808080" for brush in
                   self.widget.graph.scatterplot_item.data['brush']]
        self.assertTrue(all(is_grey))
        self.widget.controls.color_by_cluster.click()
        self.assertFalse(self.widget.controls.attr_color.isEnabled())
        is_not_grey = [brush.color().name() != "#808080" for brush in
                       self.widget.graph.scatterplot_item.data['brush']]
        self.assertTrue(any(is_not_grey))

    def test_attr_label_metas(self):
        self.send_signal(self.widget.Inputs.data, self.data)
        self.wait_until_stop_blocking()
        simulate.combobox_activate_item(self.widget.controls.attr_label,
                                        self.data.domain[-1].name)

    def test_attr_models(self):
        self.send_signal(self.widget.Inputs.data, self.data)
        controls = self.widget.controls
        self.assertEqual(len(controls.attr_color.model()), 1006)
        self.assertEqual(len(controls.attr_shape.model()), 5)
        self.assertEqual(len(controls.attr_size.model()), 1002)
        self.assertEqual(len(controls.attr_label.model()), 1008)
        for var in self.data.domain.variables + self.data.domain.metas:
            self.assertIn(var, controls.attr_label.model())
            if var.is_continuous:
                self.assertIn(var, controls.attr_color.model())
                self.assertIn(var, controls.attr_size.model())
                self.assertNotIn(var, controls.attr_shape.model())
            if var.is_discrete:
                self.assertIn(var, controls.attr_color.model())
                self.assertNotIn(var, controls.attr_size.model())
                self.assertIn(var, controls.attr_shape.model())

    def test_subset_data_color(self):
        self.send_signal(self.widget.Inputs.data, self.data)
        self.wait_until_stop_blocking()
        self.send_signal(self.widget.Inputs.data_subset, self.data[:10])
        subset = [brush.color().name() != "#000000" for brush in
                  self.widget.graph.scatterplot_item.data['brush'][:10]]
        other = [brush.color().name() == "#000000" for brush in
                 self.widget.graph.scatterplot_item.data['brush'][10:]]
        self.assertTrue(all(subset))
        self.assertTrue(all(other))


if __name__ == "__main__":
    unittest.main()
