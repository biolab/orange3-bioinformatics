import os

from AnyQt.QtCore import (  # , QItemSelectionRange, QItemSelection
    Qt,
    QItemSelectionModel,
)

from orangewidget.tests.base import WidgetTest

from Orange.data import Table, Domain, StringVariable

# from orangecontrib.bioinformatics.geneset import GeneSets
from orangecontrib.bioinformatics.widgets.OWGeneSets import OWGeneSets, GeneSetsModel
from orangecontrib.bioinformatics.widgets.utils.data import TableAnnotation


class TestOWGeneSets(WidgetTest):
    def setUp(self) -> None:
        self._path = os.path.join(os.path.dirname(__file__), 'data')
        self.widget = self.create_widget(OWGeneSets)

        domain = Domain([], metas=[StringVariable('name'), StringVariable('Entrez ID')])
        self.data = Table.from_list(
            domain, [['CA1', '759'], ['CA2', '760'], ['CA3', '761']]
        )
        self.data.attributes[TableAnnotation.tax_id] = '9606'
        self.data.attributes[TableAnnotation.gene_as_attr_name] = False
        self.data.attributes[TableAnnotation.gene_id_column] = 'Entrez ID'

    def test_gene_sets_loaded(self):
        self.wait_until_finished()
        self.assertTrue(self.widget.filter_proxy_model.sourceModel().rowCount())
        self.assertTrue(self.widget.filter_proxy_model.sourceModel().columnCount())

    # TODO: return to this test later.
    # on CI this test fails consistently on linux and sometimes for macOS platform,
    # while on windows runs without issues.
    # When running tests locally there is no problem (macOS)
    # def test_gene_set_selection(self):
    #     selection = QItemSelection()
    #
    #     selection.append(
    #         QItemSelectionRange(
    #             self.widget.view.model().index(0, 0),
    #             self.widget.view.model().index(1, 0),
    #         )
    #     )
    #
    #     selection_model = self.widget.view.selectionModel()
    #     selection_flags = QItemSelectionModel.Select | QItemSelectionModel.Rows
    #     selection_model.select(selection, selection_flags)
    #
    #     output = self.get_output(self.widget.Outputs.gene_sets)
    #     self.assertIsInstance(output, GeneSets)
    #     self.assertTrue(len(output) == 2)

    def test_input_data(self):
        self.send_signal(self.widget.Inputs.data, self.data)
        self.wait_until_finished()

        first_row = self.widget.view.model().index(0, 0)
        tenth_row = self.widget.view.model().index(10, 0)

        selection_model = self.widget.view.selectionModel()
        selection_flags = QItemSelectionModel.Select | QItemSelectionModel.Rows
        selection_model.select(first_row, selection_flags)
        selection_model.select(tenth_row, selection_flags)

        selection = selection_model.selectedRows(column=GeneSetsModel.mapped_genes)

        # all genes from the data are found in first gene set
        self.assertTrue(selection[0].data(Qt.DisplayRole))
        # non genes map to 10th gene set
        self.assertFalse(selection[1].data(Qt.DisplayRole))
