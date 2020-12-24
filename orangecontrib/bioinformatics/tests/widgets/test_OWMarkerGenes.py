import unittest
from unittest.mock import Mock

import numpy as np
from orangewidget.tests.base import WidgetTest

from AnyQt.QtCore import QModelIndex, QItemSelectionModel

from Orange.data import Table
from Orange.widgets.tests.utils import simulate

from orangecontrib.bioinformatics.utils import serverfiles
from orangecontrib.bioinformatics.widgets.OWMarkerGenes import SERVER_FILES_DOMAIN, TreeItem, TreeView, OWMarkerGenes


class TestTreeItem(unittest.TestCase):
    def setUp(self) -> None:
        self.item1 = TreeItem("item 1", False, None, None)
        self.item2 = TreeItem("item 2", False, None, None)
        self.item3 = TreeItem("item 3", False, None, None)

    def test_change_parent(self):
        self.item1.change_parent(self.item2)
        self.assertEqual(self.item1.parentItem, self.item2)
        self.assertIn(self.item1, self.item2.childItems)

    def test_change_parent_none(self):
        self.item1.change_parent(None)
        self.assertEqual(self.item1.parentItem, None)
        self.assertListEqual([], self.item2.childItems)

    def change_parents(self):
        TreeItem.change_parents([self.item1, self.item2], self.item3)
        self.assertEqual(self.item1.parentItem, self.item3)
        self.assertIn(self.item1, self.item3.childItems)
        self.assertEqual(self.item2.parentItem, self.item3)
        self.assertIn(self.item2, self.item3.childItems)

    def test_append_child(self):
        self.item1.append_child(self.item3)
        self.assertIn(self.item3, self.item1.childItems)

    def test_remove_children(self):
        TreeItem.change_parents([self.item1, self.item2], self.item3)
        self.assertListEqual([self.item1, self.item2], self.item3.childItems)

        s = self.item3.remove_children(-1, 1)  # should return False (it is not possible)
        self.assertFalse(s)
        self.assertListEqual([self.item1, self.item2], self.item3.childItems)

        s = self.item3.remove_children(2, 1)  # should return False (it is not possible)
        self.assertFalse(s)
        self.assertListEqual([self.item1, self.item2], self.item3.childItems)

        s = self.item3.remove_children(1, 2)  # should return False (it is not possible)
        self.assertFalse(s)
        self.assertListEqual([self.item1, self.item2], self.item3.childItems)

        s = self.item3.remove_children(0, 1)
        self.assertTrue(s)
        self.assertListEqual([self.item2], self.item3.childItems)

        s = self.item3.remove_children(0, 1)
        self.assertTrue(s)
        self.assertListEqual([], self.item3.childItems)

        s = self.item3.remove_children(0, 1)
        self.assertFalse(s)
        self.assertListEqual([], self.item3.childItems)

        # test remove multiple
        TreeItem.change_parents([self.item1, self.item2], self.item3)
        s = self.item3.remove_children(0, 2)
        self.assertTrue(s)
        self.assertListEqual([], self.item3.childItems)

    def test_children_count(self):
        self.assertEqual(0, self.item3.children_count())

        TreeItem.change_parents([self.item1, self.item2], self.item3)
        self.assertEqual(2, self.item3.children_count())

        self.item3.remove_children(0, 1)
        self.assertEqual(1, self.item3.children_count())

        self.item3.remove_children(0, 1)
        self.assertEqual(0, self.item3.children_count())

    def test_row(self):
        # when no parent it is assumed that the item is a root item - its row is 0
        self.assertEqual(0, self.item1.row())
        self.assertEqual(0, self.item2.row())
        self.assertEqual(0, self.item3.row())

        # add parent
        TreeItem.change_parents([self.item1, self.item2], self.item3)
        self.assertEqual(0, self.item1.row())
        self.assertEqual(1, self.item2.row())
        self.assertEqual(0, self.item3.row())

        # different order
        self.item3.remove_children(0, 2)
        TreeItem.change_parents([self.item2, self.item1], self.item3)
        self.assertEqual(1, self.item1.row())
        self.assertEqual(0, self.item2.row())
        self.assertEqual(0, self.item3.row())

    def test_child_from_name_index(self):
        self.assertTupleEqual((None, None), self.item3.child_from_name_index("item 1"))
        TreeItem.change_parents([self.item1, self.item2], self.item3)

        self.assertTupleEqual((0, self.item1), self.item3.child_from_name_index("item 1"))
        self.assertTupleEqual((1, self.item2), self.item3.child_from_name_index("item 2"))

        self.item3.remove_children(0, 2)
        TreeItem.change_parents([self.item2, self.item1], self.item3)

        self.assertTupleEqual((1, self.item1), self.item3.child_from_name_index("item 1"))
        self.assertTupleEqual((0, self.item2), self.item3.child_from_name_index("item 2"))

    def test_contains_text(self):
        # simulate RowItem with dictionary
        self.item1.data_row = {"Name": "item1", "Cell Type": "type1", "Function": "function1"}
        self.item2.data_row = {"Name": "item2", "Cell Type": "type2", "Function": "function2"}
        filter_columns = ["Name", "Cell Type", "Function"]
        TreeItem.change_parents([self.item2, self.item1], self.item3)

        self.assertTrue(self.item1.contains_text("item1", filter_columns))
        self.assertFalse(self.item2.contains_text("item1", filter_columns))
        self.assertFalse(self.item3.contains_text("item1", filter_columns))

        self.assertTrue(self.item1.contains_text("type1", filter_columns))
        self.assertFalse(self.item2.contains_text("type1", filter_columns))
        self.assertFalse(self.item3.contains_text("type1", filter_columns))

        self.assertTrue(self.item1.contains_text("function1", filter_columns))
        self.assertFalse(self.item2.contains_text("function1", filter_columns))
        self.assertFalse(self.item3.contains_text("function1", filter_columns))

        self.assertFalse(self.item1.contains_text("item2", filter_columns))
        self.assertTrue(self.item2.contains_text("item2", filter_columns))
        self.assertFalse(self.item3.contains_text("item2", filter_columns))

        self.assertFalse(self.item1.contains_text("type2", filter_columns))
        self.assertTrue(self.item2.contains_text("type2", filter_columns))
        self.assertFalse(self.item3.contains_text("type2", filter_columns))

        self.assertFalse(self.item1.contains_text("function2", filter_columns))
        self.assertTrue(self.item2.contains_text("function2", filter_columns))
        self.assertFalse(self.item3.contains_text("function2", filter_columns))

        # test case sensitivity
        self.assertTrue(self.item1.contains_text("Item1", filter_columns))
        self.assertFalse(self.item2.contains_text("Item1", filter_columns))

        self.assertTrue(self.item1.contains_text("Type1", filter_columns))
        self.assertFalse(self.item2.contains_text("Type1", filter_columns))

        self.assertFalse(self.item1.contains_text("Item2", filter_columns))
        self.assertTrue(self.item2.contains_text("Item2", filter_columns))

        self.assertFalse(self.item1.contains_text("Type2", filter_columns))
        self.assertTrue(self.item2.contains_text("Type2", filter_columns))

    def test_get_data_rows(self):
        TreeItem.change_parents([self.item2, self.item1], self.item3)

        self.item1.data_row = "dr1"
        self.item2.data_row = "dr2"
        self.item3.data_row = "dr3"

        self.assertListEqual(sorted(["dr1", "dr2", "dr3"]), sorted(self.item3.get_data_rows()))
        self.assertListEqual(["dr1"], self.item1.get_data_rows())
        self.assertListEqual(["dr2"], self.item2.get_data_rows())

        self.item1.data_row = "dr1"
        self.item2.data_row = "dr2"
        self.item3.data_row = None

        self.assertListEqual(sorted(["dr1", "dr2"]), sorted(self.item3.get_data_rows()))
        self.assertListEqual(["dr1"], self.item1.get_data_rows())
        self.assertListEqual(["dr2"], self.item2.get_data_rows())

        self.item1.data_row = "dr1"
        self.item2.data_row = None
        self.item3.data_row = None

        self.assertListEqual(["dr1"], self.item3.get_data_rows())
        self.assertListEqual(["dr1"], self.item1.get_data_rows())
        self.assertListEqual([], self.item2.get_data_rows())

    def test_count_marker_genes(self):
        TreeItem.change_parents([self.item2, self.item1], self.item3)

        self.item1.data_row = "dr1"
        self.item2.data_row = "dr2"
        self.item3.data_row = "dr3"

        self.assertEqual(3, self.item3.count_marker_genes())
        self.assertEqual(1, self.item1.count_marker_genes())
        self.assertEqual(1, self.item2.count_marker_genes())

        self.item1.data_row = "dr1"
        self.item2.data_row = "dr2"
        self.item3.data_row = None

        self.assertEqual(2, self.item3.count_marker_genes())
        self.assertEqual(1, self.item1.count_marker_genes())
        self.assertEqual(1, self.item2.count_marker_genes())

        self.item1.data_row = "dr1"
        self.item2.data_row = None
        self.item3.data_row = None

        self.assertEqual(1, self.item3.count_marker_genes())
        self.assertEqual(1, self.item1.count_marker_genes())
        self.assertEqual(0, self.item2.count_marker_genes())


class TestOWMarkerGenes(WidgetTest):
    @classmethod
    def setUpClass(cls):
        """ Code executed only once for all tests """
        super().setUpClass()
        file_name = "panglao_gene_markers.tab"
        serverfiles.update(SERVER_FILES_DOMAIN, file_name)
        file_path = serverfiles.localpath_download(SERVER_FILES_DOMAIN, file_name)
        cls.panglao = Table.from_file(file_path)

        file_name = "cellMarker_gene_markers.tab"
        serverfiles.update(SERVER_FILES_DOMAIN, file_name)
        file_path = serverfiles.localpath_download(SERVER_FILES_DOMAIN, file_name)
        cls.cell_markers = Table.from_file(file_path)

    def test_minimum_size(self):
        pass

    def setUp(self):
        self.widget = self.create_widget(OWMarkerGenes)

    def move_top_element(self, source_view: TreeView) -> None:
        """
        This function moves top element form the source_view to another.

        Parameters
        ----------
        source_view
            The view from which we move top element.
        """
        parent_index = source_view.model().index(0, 0)
        source_view.selectionModel().select(parent_index, QItemSelectionModel.Select | QItemSelectionModel.Rows)
        self.widget.move_button.click()

    def move_first_child(self, source_view: TreeView) -> None:
        """
        This function moves the first child element from source_view to another.

        Parameters
        ----------
        source_view
            The view from which we move top element.
        """
        parent_index = source_view.model().index(0, 0)
        gene_index = source_view.model().index(0, 0, parent_index)
        source_view.selectionModel().select(gene_index, QItemSelectionModel.Select | QItemSelectionModel.Rows)
        self.widget.move_button.click()

    def test_data_not_empty(self):
        """
        Test that data are present
        """
        w = self.widget
        self.assertIsNotNone(w.data)

    def test_available_sources(self):
        """
        Test if all sources retrieved. Currently there are two sources available. When more available unittest will be
        changed.
        """
        self.assertListEqual(["CellMarker", "Panglao"], list(self.widget.available_sources.keys()))

    def test_source_changed(self):
        """
        Test changing the source of the data.
        """
        # cell marker data
        self.assertEqual("Panglao", self.widget.db_source_cb.currentText())
        self.assertTrue(isinstance(self.widget.data, Table))
        self.assertEqual(len(self.panglao), len(self.widget.data))

        # Panglao data
        simulate.combobox_activate_index(self.widget.controls.source_index, 1)
        self.assertEqual("CellMarker", self.widget.db_source_cb.currentText())
        self.assertTrue(isinstance(self.widget.data, Table))
        self.assertEqual(len(self.cell_markers), len(self.widget.data))

        # cell marker data
        simulate.combobox_activate_index(self.widget.controls.source_index, 0)
        self.assertEqual("Panglao", self.widget.db_source_cb.currentText())
        self.assertTrue(isinstance(self.widget.data, Table))
        self.assertEqual(len(self.panglao), len(self.widget.data))

    def test_organism_changed(self):
        """
        Test changing organisms in combination with changing the source.
        """
        simulate.combobox_activate_index(self.widget.controls.source_index, 0)
        simulate.combobox_activate_index(self.widget.controls.organism_index, 0)
        self.assertEqual("Panglao", self.widget.db_source_cb.currentText())
        self.assertEqual("Human", self.widget.group_cb.currentText())

        human_rows = self.widget.data.get_column_view("Organism")[0] == "Human"
        cell_types = self.widget.data.get_column_view("Cell Type")[0]

        model = self.widget.available_markers_view.model().sourceModel()
        self.assertEqual(sum(self.panglao.get_column_view("Organism")[0] == "Human"), len(model))
        self.assertEqual(len(np.unique(cell_types[human_rows])), len(model.rootItem.childItems))

        simulate.combobox_activate_index(self.widget.controls.source_index, 0)
        simulate.combobox_activate_index(self.widget.controls.organism_index, 1)
        self.assertEqual("Panglao", self.widget.db_source_cb.currentText())
        self.assertEqual("Mouse", self.widget.group_cb.currentText())

        model = self.widget.available_markers_view.model().sourceModel()
        self.assertEqual(sum(self.panglao.get_column_view("Organism")[0] == "Mouse"), len(model))
        self.assertEqual(len(np.unique(cell_types[~human_rows])), len(model.rootItem.childItems))

        simulate.combobox_activate_index(self.widget.controls.source_index, 1)
        simulate.combobox_activate_index(self.widget.controls.organism_index, 0)
        self.assertEqual("CellMarker", self.widget.db_source_cb.currentText())
        self.assertEqual("Human", self.widget.group_cb.currentText())

        human_rows = self.widget.data.get_column_view("Organism")[0] == "Human"
        cell_types = self.widget.data.get_column_view("Cell Type")[0]

        model = self.widget.available_markers_view.model().sourceModel()
        self.assertEqual(sum(self.cell_markers.get_column_view("Organism")[0] == "Human"), len(model))
        self.assertEqual(len(np.unique(cell_types[human_rows])), len(model.rootItem.childItems))

        simulate.combobox_activate_index(self.widget.controls.source_index, 1)
        simulate.combobox_activate_index(self.widget.controls.organism_index, 1)
        self.assertEqual("CellMarker", self.widget.db_source_cb.currentText())
        self.assertEqual("Mouse", self.widget.group_cb.currentText())

        model = self.widget.available_markers_view.model().sourceModel()
        self.assertEqual(sum(self.cell_markers.get_column_view("Organism")[0] == "Mouse"), len(model))
        self.assertEqual(len(np.unique(cell_types[~human_rows])), len(model.rootItem.childItems))

    def test_group_by_changed(self):
        """
        Test changing the group by parameter (Cell Type, Function)
        """
        # group by cell type
        simulate.combobox_activate_index(self.widget.controls.selected_root_attribute, 0)

        human_rows = self.widget.data.get_column_view("Organism")[0] == "Human"
        cell_types = np.unique(self.widget.data.get_column_view("Cell Type")[0][human_rows])
        model = self.widget.available_markers_view.model().sourceModel()
        self.assertEqual(len(cell_types), len(model.rootItem.childItems))

        # group by function
        simulate.combobox_activate_index(self.widget.controls.selected_root_attribute, 1)

        functions = np.unique(self.widget.data.get_column_view("Function")[0][human_rows])
        model = self.widget.available_markers_view.model().sourceModel()
        self.assertEqual(len(functions), len(model.rootItem.childItems))

    def test_move_elements(self):
        """
        Move elements between booth views and observe if numbers matches.
        This test just move the single element.
        """
        # select first three elements
        src_view = self.widget.available_markers_view
        dest_view = self.widget.selected_markers_view
        num_items_source_view = len(src_view.model().sourceModel())
        num_root_items = src_view.model().rowCount(QModelIndex())

        # select first element
        index = src_view.model().index(0, 0)
        src_view.selectionModel().select(index, QItemSelectionModel.Select | QItemSelectionModel.Rows)
        num_items = src_view.model().rowCount(index)  # number of items under selected item

        # move element to the TreeView with selected genes
        self.widget.move_button.click()
        self.assertEqual(len(dest_view.model().sourceModel()), num_items)
        output = self.get_output(self.widget.Outputs.genes)
        self.assertEqual(num_items, len(output))
        self.assertEqual(dest_view.model().rowCount(QModelIndex()), 1)  # one root element added
        index = dest_view.model().index(0, 0)
        self.assertEqual(dest_view.model().rowCount(index), num_items)

        self.assertEqual(src_view.model().rowCount(QModelIndex()), num_root_items - 1)  # one item is remove from view

        # move everything back
        dest_view.selectionModel().select(index, QItemSelectionModel.Select | QItemSelectionModel.Rows)
        self.widget.move_button.click()

        self.assertEqual(num_items_source_view, len(src_view.model().sourceModel()))
        self.assertEqual(src_view.model().rowCount(QModelIndex()), num_root_items)

    def test_move_more_elements(self):
        """
        Move elements between booth views and observe if numbers matches.
        This test keep moving more elements.
        """
        # select first three elements
        src_view = self.widget.available_markers_view
        dest_view = self.widget.selected_markers_view
        num_items_source_view = len(src_view.model().sourceModel())
        num_root_items = src_view.model().rowCount(QModelIndex())

        # select first element
        index = src_view.model().index(0, 0)
        index1 = src_view.model().index(1, 0)
        src_view.selectionModel().select(index, QItemSelectionModel.Select | QItemSelectionModel.Rows)
        src_view.selectionModel().select(index1, QItemSelectionModel.Select | QItemSelectionModel.Rows)
        num_items = src_view.model().rowCount(index)
        num_items1 = src_view.model().rowCount(index1)

        # move elements to the TreeView with selected genes
        self.widget.move_button.click()
        self.assertEqual(len(dest_view.model().sourceModel()), num_items + num_items1)
        output = self.get_output(self.widget.Outputs.genes)
        self.assertEqual(num_items + num_items1, len(output))
        self.assertEqual(dest_view.model().rowCount(QModelIndex()), 2)  # two root element added
        index = dest_view.model().index(0, 0)
        index1 = dest_view.model().index(1, 0)
        self.assertEqual(
            dest_view.model().rowCount(index) + dest_view.model().rowCount(index1), num_items + num_items1
        )

        self.assertEqual(src_view.model().rowCount(QModelIndex()), num_root_items - 2)

        # move everything back one by one
        dest_view.selectionModel().select(index, QItemSelectionModel.Select | QItemSelectionModel.Rows)
        self.widget.move_button.click()
        index1 = dest_view.model().index(0, 0)
        self.assertEqual(dest_view.model().rowCount(index1), num_items1)
        output = self.get_output(self.widget.Outputs.genes)
        self.assertEqual(num_items1, len(output))
        self.assertEqual(len(dest_view.model().sourceModel()), num_items1)
        self.assertEqual(src_view.model().rowCount(QModelIndex()), num_root_items - 1)

        # move the other one back
        dest_view.selectionModel().select(index1, QItemSelectionModel.Select | QItemSelectionModel.Rows)
        self.widget.move_button.click()
        self.assertEqual(num_items_source_view, len(src_view.model().sourceModel()))
        self.assertEqual(0, len(dest_view.model().sourceModel()))
        self.assertEqual(src_view.model().rowCount(QModelIndex()), num_root_items)
        self.assertEqual(dest_view.model().rowCount(QModelIndex()), 0)

    def test_move_genes(self):
        """
        Move elements between booth views and observe if numbers matches.
        In previous test we were moving just root elements. In this test we will move children elements.
        """
        src_view = self.widget.available_markers_view
        dest_view = self.widget.selected_markers_view
        num_items_source_view = len(src_view.model().sourceModel())
        num_root_items = src_view.model().rowCount(QModelIndex())

        # select first element
        parent_index = src_view.model().index(0, 0)
        index = src_view.model().index(0, 0, parent_index)
        index1 = src_view.model().index(1, 0, parent_index)

        src_view.selectionModel().select(index, QItemSelectionModel.Select | QItemSelectionModel.Rows)
        src_view.selectionModel().select(index1, QItemSelectionModel.Select | QItemSelectionModel.Rows)

        # move elements to the TreeView with selected genes
        self.widget.move_button.click()
        self.assertEqual(len(dest_view.model().sourceModel()), 2)
        output = self.get_output(self.widget.Outputs.genes)
        self.assertEqual(2, len(output))
        self.assertEqual(dest_view.model().rowCount(QModelIndex()), 1)  # one parent element added
        parent_index = dest_view.model().index(0, 0)
        self.assertEqual(dest_view.model().rowCount(parent_index), 2)
        self.assertEqual(src_view.model().rowCount(QModelIndex()), num_root_items)

        # move everything back one by one
        index = dest_view.model().index(0, 0, parent_index)
        dest_view.selectionModel().select(index, QItemSelectionModel.Select | QItemSelectionModel.Rows)
        self.widget.move_button.click()

        self.assertEqual(dest_view.model().rowCount(parent_index), 1)
        output = self.get_output(self.widget.Outputs.genes)
        self.assertEqual(1, len(output))
        self.assertEqual(len(dest_view.model().sourceModel()), 1)
        self.assertEqual(src_view.model().rowCount(QModelIndex()), num_root_items)

        # move the other one back
        index1 = dest_view.model().index(0, 0, parent_index)
        dest_view.selectionModel().select(index1, QItemSelectionModel.Select | QItemSelectionModel.Rows)
        self.widget.move_button.click()
        self.assertEqual(num_items_source_view, len(src_view.model().sourceModel()))
        self.assertEqual(src_view.model().rowCount(QModelIndex()), num_root_items)

    def test_description(self):
        """
        Test correctness of the gene description showing in the bottom of widget.
        """
        src_view = self.widget.available_markers_view
        dest_view = self.widget.selected_markers_view

        not_gene = "Select a gene to see information."
        first_gene_info = "Gene name: AKR1C3"

        # first test in the source view
        parent_index = src_view.model().index(0, 0)
        gene_index = src_view.model().index(0, 0, parent_index)
        gene_index1 = src_view.model().index(1, 0, parent_index)

        # no description for parent
        src_view.selectionModel().select(parent_index, QItemSelectionModel.Select | QItemSelectionModel.Rows)
        self.assertEqual(not_gene, self.widget.descriptionlabel.toPlainText())

        # no description for multiple genes
        src_view.selectionModel().select(gene_index, QItemSelectionModel.Select | QItemSelectionModel.Rows)
        src_view.selectionModel().select(gene_index1, QItemSelectionModel.Select | QItemSelectionModel.Rows)
        self.assertEqual(not_gene, self.widget.descriptionlabel.toPlainText())

        # description for first gene
        src_view.selectionModel().select(gene_index, QItemSelectionModel.ClearAndSelect | QItemSelectionModel.Rows)
        self.assertIn(first_gene_info, self.widget.descriptionlabel.toPlainText())

        # move elements
        src_view.selectionModel().select(gene_index, QItemSelectionModel.Select | QItemSelectionModel.Rows)
        src_view.selectionModel().select(gene_index1, QItemSelectionModel.Select | QItemSelectionModel.Rows)
        self.widget.move_button.click()

        # test in the dest view
        parent_index = dest_view.model().index(0, 0)
        gene_index = dest_view.model().index(0, 0, parent_index)
        gene_index1 = dest_view.model().index(1, 0, parent_index)

        # no description for parent
        dest_view.selectionModel().select(parent_index, QItemSelectionModel.Select | QItemSelectionModel.Rows)
        self.assertEqual(not_gene, self.widget.descriptionlabel.toPlainText())

        # no description for multiple genes
        dest_view.selectionModel().select(gene_index, QItemSelectionModel.Select | QItemSelectionModel.Rows)
        dest_view.selectionModel().select(gene_index1, QItemSelectionModel.Select | QItemSelectionModel.Rows)
        self.assertEqual(not_gene, self.widget.descriptionlabel.toPlainText())

        # description for first gene
        dest_view.selectionModel().select(gene_index, QItemSelectionModel.ClearAndSelect | QItemSelectionModel.Rows)
        self.assertIn(first_gene_info, self.widget.descriptionlabel.toPlainText())

    def test_output_info(self):
        """
        Test whether the output info in a status bar is set correctly.
        """
        # mocks
        output_sum = self.widget.info.set_output_summary = Mock()

        # views
        src_view = self.widget.available_markers_view
        dest_view = self.widget.selected_markers_view

        # move elements
        parent_index = src_view.model().index(0, 0)
        gene_index = src_view.model().index(0, 0, parent_index)
        gene_index1 = src_view.model().index(1, 0, parent_index)
        src_view.selectionModel().select(gene_index, QItemSelectionModel.Select | QItemSelectionModel.Rows)
        self.widget.move_button.click()

        output_sum.assert_called_with("Selected: 1")

        # move one more
        src_view.selectionModel().select(gene_index1, QItemSelectionModel.Select | QItemSelectionModel.Rows)
        self.widget.move_button.click()

        output_sum.assert_called_with("Selected: 2")

        # move complete parent
        src_view.selectionModel().select(parent_index, QItemSelectionModel.Select | QItemSelectionModel.Rows)
        self.widget.move_button.click()

        output_sum.assert_called_with(f"Selected: {len(dest_view.model().sourceModel())}")

        # move everthing back
        parent_index = dest_view.model().index(0, 0)
        dest_view.selectionModel().select(parent_index, QItemSelectionModel.Select | QItemSelectionModel.Rows)
        self.widget.move_button.click()

        output_sum.assert_called_with("Selected: 0")

    def test_extend_collapse_both_views(self):
        """
        Test whether the item extends in the other column when extended in one.
        """
        src_view = self.widget.available_markers_view
        dest_view = self.widget.selected_markers_view

        # I do not know how to simulate click so I just call the function that is called on click
        index = src_view.model().index(0, 0)
        src_view.expand_in_other_view(index, is_expanded=True)
        tree_node = src_view.model().sourceModel().node_from_index(src_view.model().mapToSource(index))
        self.assertTrue(tree_node.expanded)

        # move one child left to check simultaneous expand
        self.move_first_child(src_view)
        index_left = src_view.model().index(0, 0)
        index_right = dest_view.model().index(0, 0)
        tree_node_left = src_view.model().sourceModel().node_from_index(src_view.model().mapToSource(index_left))
        tree_node_right = dest_view.model().sourceModel().node_from_index(dest_view.model().mapToSource(index_right))

        src_view.expand_in_other_view(index_left, is_expanded=True)
        self.assertTrue(tree_node_left.expanded)
        self.assertTrue(tree_node_right.expanded)
        self.assertTrue(dest_view.isExpanded(index_right))

        # change back to not expanded
        src_view.expand_in_other_view(index_left, is_expanded=False)
        self.assertFalse(tree_node_left.expanded)
        self.assertFalse(tree_node_right.expanded)
        self.assertFalse(dest_view.isExpanded(index_right))

        dest_view.expand_in_other_view(index_right, is_expanded=True)
        self.assertTrue(tree_node_left.expanded)
        self.assertTrue(tree_node_right.expanded)
        self.assertTrue(src_view.isExpanded(index_left))

        # change back to not expanded
        dest_view.expand_in_other_view(index_right, is_expanded=False)
        self.assertFalse(tree_node_left.expanded)
        self.assertFalse(tree_node_right.expanded)
        self.assertFalse(dest_view.isExpanded(index_left))


if __name__ == '__main__':
    unittest.main()
