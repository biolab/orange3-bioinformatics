import unittest
import os

from Orange.data import Table, Domain, StringVariable
from Orange.widgets.tests.base import WidgetTest
from orangecontrib.bioinformatics.widgets.OWHomologs import OWHomologs
from orangecontrib.bioinformatics.widgets.utils.data import TableAnnotation


class TestOWMHomologs(WidgetTest):
    def setUp(self):
        self.widget = self.create_widget(OWHomologs)

    def test_rows(self):
        name = StringVariable("Entrez ID")
        id = StringVariable("Name")
        attributes = {TableAnnotation.gene_as_attr_name: False, TableAnnotation.tax_id: "10090",
                      TableAnnotation.gene_id_column: "Entrez ID"}

        domain = Domain([], (), [id, name])

        data = [["CD4", "920"],
                ["CD44", "960"],
                ["CD48", "962"],
                ["CD47", "961"],
                 ["CD46", "4179"]]
        genes = Table(domain, data)
        genes.attributes = attributes
        self.widget.target_tax = "9606"
        self.send_signal(self.widget.Inputs.data, genes)

        out_data = self.get_output("Genes", self.widget)

        print(out_data)



if __name__ == '__main__':
    unittest.main()
