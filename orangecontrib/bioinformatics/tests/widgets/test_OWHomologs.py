import unittest

from orangewidget.tests.base import WidgetTest

from Orange.data import Table, Domain, StringVariable, ContinuousVariable

from orangecontrib.bioinformatics.widgets.OWHomologs import HOMOLOG_ID, OWHomologs
from orangecontrib.bioinformatics.widgets.utils.data import TableAnnotation


class TestOWMHomologs(WidgetTest):
    def setUp(self):
        self.widget = self.create_widget(OWHomologs)
        self.expected_results = ['12504', '12505', '12506', '16423', '17221']
        id = StringVariable("Entrez ID")
        name = StringVariable("Name")
        attributes_rows = {
            TableAnnotation.gene_as_attr_name: False,
            TableAnnotation.tax_id: "9606",
            TableAnnotation.gene_id_column: "Entrez ID",
        }
        attributes_columns = {
            TableAnnotation.gene_as_attr_name: True,
            TableAnnotation.tax_id: "9606",
            TableAnnotation.gene_id_attribute: "Entrez ID",
        }

        domain = Domain([], (), [name, id])

        data_row = [["CD4", "920"], ["CD44", "960"], ["CD48", "962"], ["CD47", "961"], ["CD46", "4179"]]
        self.genes_rows = Table(domain, data_row)
        self.genes_rows.attributes = attributes_rows

        domain_columns = Domain([ContinuousVariable(name) for name, id in data_row])
        for col, id in zip(domain_columns.attributes, data_row):
            col.attributes["Entrez ID"] = id[1]

        self.genes_columns = Table(domain_columns, [[1, 2, 1, 2, 1]])
        self.genes_columns.attributes = attributes_columns

    def test_homologs_by_rows(self):

        self.widget.auto_commit = False
        self.widget.target_organism_change(self.widget.taxonomy_ids.index("10090"))

        self.send_signal(self.widget.Inputs.data, self.genes_rows)
        self.widget.auto_commit = True
        self.widget.commit()

        out_data = self.get_output("Genes", self.widget)
        self.assertEqual(len(self.genes_rows), len(out_data))
        mouse_ids = list(out_data.get_column_view(HOMOLOG_ID)[0])
        self.assertListEqual(mouse_ids, self.expected_results)

    def test_auto_commit(self):
        self.send_signal(self.widget.Inputs.data, self.genes_rows)
        out_data = self.get_output("Genes", self.widget)
        self.assertEqual(out_data, None)

    def test_homologs_by_column(self):
        self.widget.auto_commit = False
        self.widget.target_organism_change(self.widget.taxonomy_ids.index("10090"))

        self.send_signal(self.widget.Inputs.data, self.genes_columns)
        self.widget.auto_commit = True
        self.widget.commit()

        out_data = self.get_output("Genes", self.widget)
        mouse_ids = [(list(att.attributes.values())[2]) for att in out_data.domain.attributes]

        self.assertEqual(len(self.genes_rows), len(mouse_ids))

        self.assertListEqual(mouse_ids, self.expected_results)


if __name__ == '__main__':
    unittest.main()
