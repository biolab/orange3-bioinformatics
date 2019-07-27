# Test methods with long descriptive names can omit docstrings
# pylint: disable=missing-docstring,arguments-differ
from itertools import chain
import unittest
from unittest.mock import Mock

import numpy as np

from Orange.data import Table, Variable, Domain
from Orange.data.filter import FilterString, Values
from Orange.projection import PCA
from Orange.widgets.tests.base import WidgetTest, WidgetOutputsTestMixin, \
    ProjectionWidgetTestMixin, simulate
from Orange.widgets.unsupervised.owtsne import OWtSNE

from orangecontrib.bioinformatics.utils import serverfiles
from orangecontrib.bioinformatics.widgets.OWAnnotateProjection import \
    OWAnnotateProjection
from orangecontrib.bioinformatics.widgets.utils.data import TAX_ID


class TestOWAnnotateProjection(WidgetTest, ProjectionWidgetTestMixin,
                               WidgetOutputsTestMixin):
    @classmethod
    def setUpClass(cls):
        super().setUpClass()
        Variable._clear_all_caches()
        cls._init_data()
        cls.signal_name = "Reference Data"
        cls.signal_data = cls.data
        cls.same_input_output_domain = False

        genes_path = serverfiles.localpath_download(
            "marker_genes", "panglao_gene_markers.tab")
        filter_ = FilterString("Organism", FilterString.Equal, "Human")
        cls.genes = Values([filter_])(Table(genes_path))
        cls.genes.attributes[TAX_ID] = "9606"

    @classmethod
    def _init_data(cls):
        data_path = "https://datasets.orange.biolab.si/sc/aml-1k.tab.gz"
        table_data = Table(data_path)
        table_data.attributes[TAX_ID] = "9606"

        ref_data = table_data[::2]
        pca = PCA(n_components=2)
        model = pca(ref_data)
        proj = model(ref_data)
        domain = Domain(ref_data.domain.attributes, ref_data.domain.class_vars,
                        chain(ref_data.domain.metas, proj.domain.attributes))
        cls.data = ref_data.transform(domain)
        cls.reference_data = ref_data
        cls.secondary_data = table_data[1:200:2]

    def setUp(self):
        self.widget = self.create_widget(OWAnnotateProjection)

    def tearDown(self):
        self.widget.cancel()
        self.wait_until_stop_blocking()

    def test_input_secondary_data(self):
        self.send_signal(self.widget.Inputs.secondary_data,
                         self.secondary_data)
        self.wait_until_stop_blocking()
        self.assertTrue(self.widget.Error.no_reference_data.is_shown())

        self.send_signal(self.widget.Inputs.data, self.data)
        self.wait_until_stop_blocking()
        self.assertFalse(self.widget.Error.no_reference_data.is_shown())
        opts = self.widget.graph.ref_scatterplot_item.opts
        self.assertEqual(opts["pen"].color().name(), "#c8c8c8")
        self.assertEqual(opts["brush"].color().name(), "#c8c8c8")
        self.assertTrue(self.widget.graph.ref_scatterplot_item.isVisible())
        self.assertTrue(self.widget.graph.scatterplot_item.isVisible())

        self.send_signal(self.widget.Inputs.data, None)
        self.wait_until_stop_blocking()
        self.assertTrue(self.widget.Error.no_reference_data.is_shown())

        self.send_signal(self.widget.Inputs.secondary_data, None)
        self.wait_until_stop_blocking()
        self.assertFalse(self.widget.Error.no_reference_data.is_shown())

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
        n_metas = len(self.data.domain.metas)
        self.assertGreater(len(output.domain.metas), n_metas + 4)
        self.assertListEqual(
            ["Clusters", "Annotation"] +
            [m.name for m in self.data.domain.metas] + ["Selected"],
            [m.name for m in output.domain.metas[-n_metas-3:]])
        np.testing.assert_array_equal(output.X, self.data.X)
        np.testing.assert_array_equal(output.Y, self.data.Y)
        np.testing.assert_array_equal(output.metas[:, -n_metas - 1:-1],
                                      self.data.metas)

    def test_button_no_data(self):
        self.widget.run_button.click()
        self.assertEqual(self.widget.run_button.text(), "Start")

    def test_button_with_data(self):
        self.send_signal(self.widget.Inputs.genes, self.genes)
        self.send_signal(self.widget.Inputs.data, self.data)
        self.assertEqual(self.widget.run_button.text(), "Stop")
        self.wait_until_stop_blocking()
        self.assertEqual(self.widget.run_button.text(), "Start")

    def test_button_toggle(self):
        self.send_signal(self.widget.Inputs.genes, self.genes)
        self.send_signal(self.widget.Inputs.data, self.data)
        self.widget.run_button.click()
        self.assertEqual(self.widget.run_button.text(), "Resume")

    def test_plot_once(self):
        self.widget.setup_plot = Mock()
        self.widget.commit = Mock()
        self.send_signal(self.widget.Inputs.genes, self.genes)
        self.widget.setup_plot.assert_called_once()
        self.widget.setup_plot.reset_mock()
        self.widget.commit.assert_called_once()
        self.widget.commit.reset_mock()
        self.send_signal(self.widget.Inputs.data, self.data)
        self.widget.setup_plot.assert_called_once()
        self.widget.commit.assert_called_once()
        self.wait_until_stop_blocking()
        self.widget.setup_plot.reset_mock()
        self.widget.commit.reset_mock()
        self.send_signal(self.widget.Inputs.secondary_data, self.secondary_data)
        self.widget.setup_plot.assert_called_once()
        self.widget.commit.assert_called_once()

    def test_saved_selection(self):
        self.send_signal(self.widget.Inputs.genes, self.genes)
        self.send_signal(self.widget.Inputs.data, self.data)
        self.wait_until_stop_blocking()

        self.widget.graph.select_by_indices(list(range(0, len(self.data), 10)))
        settings = self.widget.settingsHandler.pack_data(self.widget)
        widget = self.create_widget(self.widget.__class__,
                                    stored_settings=settings)

        self.send_signal(widget.Inputs.genes, self.genes)
        self.send_signal(widget.Inputs.data, self.data, widget=widget)
        self.wait_until_stop_blocking(widget=widget)

        self.assertEqual(np.sum(widget.graph.selection), 50)
        np.testing.assert_equal(self.widget.graph.selection,
                                widget.graph.selection)

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
        self.assertEqual(len(controls.attr_color.model()), 1008)
        self.assertEqual(len(controls.attr_shape.model()), 5)
        self.assertEqual(len(controls.attr_size.model()), 1005)
        self.assertEqual(len(controls.attr_label.model()), 1010)
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
        self.assertRaises(AttributeError,
                          lambda: self.widget.Inputs.data_subset)

    def test_sparse_data(self):
        table = Table("iris").to_sparse()
        self.send_signal(self.widget.Inputs.data, table)

    def test_missing_embedding_data(self):
        self.send_signal(self.widget.Inputs.genes, self.genes)
        self.send_signal(self.widget.Inputs.data, self.secondary_data[::2])
        self.wait_until_stop_blocking()
        secondary_data = self.secondary_data[1::2]
        self.send_signal(self.widget.Inputs.secondary_data, secondary_data)
        self.wait_until_stop_blocking()
        self.assertFalse(self.widget.Warning.missing_compute_value.is_shown())
        secondary_data = secondary_data[:, 2:]
        self.send_signal(self.widget.Inputs.secondary_data, secondary_data)
        self.wait_until_stop_blocking()
        self.assertTrue(self.widget.Warning.missing_compute_value.is_shown())
        self.send_signal(self.widget.Inputs.secondary_data, None)
        self.wait_until_stop_blocking()
        self.assertFalse(self.widget.Warning.missing_compute_value.is_shown())

    def test_tsne_output(self):
        owtsne = self.create_widget(OWtSNE)
        self.send_signal(
            owtsne.Inputs.data, self.reference_data, widget=owtsne)
        self.wait_until_stop_blocking(widget=owtsne, wait=10000)
        tsne_output = self.get_output(owtsne.Outputs.annotated_data, owtsne)

        self.send_signal(self.widget.Inputs.genes, self.genes)
        self.send_signal(self.widget.Inputs.data, tsne_output)
        self.wait_until_stop_blocking()
        self.send_signal(self.widget.Inputs.secondary_data,
                         self.secondary_data)


if __name__ == "__main__":
    unittest.main()
