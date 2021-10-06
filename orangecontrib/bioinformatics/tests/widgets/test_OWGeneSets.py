import os
from unittest.mock import Mock

from orangewidget.tests.base import WidgetTest

from Orange.data import Table, Domain, StringVariable

from orangecontrib.bioinformatics.widgets.OWGeneSets import OWGeneSets
from orangecontrib.bioinformatics.widgets.utils.data import TableAnnotation


class TestOWGeneSets(WidgetTest):
    def setUp(self) -> None:
        self._path = os.path.join(os.path.dirname(__file__), 'data')
        self.widget = self.create_widget(OWGeneSets)

        domain = Domain([], metas=[StringVariable('name'), StringVariable('Entrez ID')])
        self.data = Table.from_list(domain, [['CA1', '759'], ['CA2', '760'], ['CA3', '761']])
        self.data.attributes[TableAnnotation.tax_id] = '9606'
        self.data.attributes[TableAnnotation.gene_as_attr_name] = False
        self.data.attributes[TableAnnotation.gene_id_column] = 'Entrez ID'

    def test_input_info(self):
        input_sum = self.widget.info.set_input_summary = Mock()

        self.send_signal(self.widget.Inputs.data, self.data)
        self.wait_until_finished()
        input_sum.assert_called_with("3", "3 unique gene names on input.\n")

        self.send_signal(self.widget.Inputs.custom_gene_sets, self.data)
        self.wait_until_finished()
        input_sum.assert_called_with("3|3", "3 unique gene names on input.\n3 marker genes in 3 sets\n")

        self.send_signal(self.widget.Inputs.data, None)
        self.send_signal(self.widget.Inputs.custom_gene_sets, None)
        self.wait_until_finished()
        input_sum.assert_called_with(self.widget.info.NoInput)

    def test_output_info(self):
        output_sum = self.widget.info.set_output_summary = Mock()

        self.send_signal(self.widget.Inputs.data, self.data)
        self.wait_until_finished()
        output_sum.assert_called_with("0", "0 genes on output.\n")

        self.send_signal(self.widget.Inputs.data, None)
        output_sum.assert_called_with(self.widget.info.NoOutput)
