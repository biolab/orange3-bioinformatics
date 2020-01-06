import os
from unittest.mock import Mock

from Orange.data import Table
from Orange.widgets.tests.base import WidgetTest

from orangecontrib.bioinformatics.widgets.OWGeneSets import OWGeneSets


class TestOWGeneSets(WidgetTest):
    def setUp(self) -> None:
        self._path = os.path.join(os.path.dirname(__file__), 'data')
        self.widget = self.create_widget(OWGeneSets)
        self.markers_file = 'markers-panglao.pkl'

        file_name = os.path.join(self._path, self.markers_file)
        self.data = Table(file_name)

    def test_input_info(self):
        input_sum = self.widget.info.set_input_summary = Mock()

        self.send_signal(self.widget.Inputs.genes, self.data)
        self.wait_until_finished()
        input_sum.assert_called_with("100", "100 unique gene names on input.\n")

        self.send_signal(self.widget.Inputs.custom_sets, self.data)
        self.wait_until_finished()
        input_sum.assert_called_with("100|100", "100 unique gene names on input.\n100 marker genes in 98 sets\n")

        self.send_signal(self.widget.Inputs.genes, None)
        self.send_signal(self.widget.Inputs.custom_sets, None)
        self.wait_until_finished()
        input_sum.assert_called_with(self.widget.info.NoInput)

        self.data.attributes = {}
        self.send_signal(self.widget.Inputs.genes, self.data)
        input_sum.assert_called_with("0", "Input data with incorrect meta data.\nUse Gene Name Matcher widget.")

    def test_output_info(self):
        output_sum = self.widget.info.set_output_summary = Mock()

        self.send_signal(self.widget.Inputs.genes, self.data)
        self.wait_until_finished()
        output_sum.assert_called_with("0", "0 genes on output.\n")

        self.send_signal(self.widget.Inputs.genes, None)
        output_sum.assert_called_with(self.widget.info.NoOutput)
