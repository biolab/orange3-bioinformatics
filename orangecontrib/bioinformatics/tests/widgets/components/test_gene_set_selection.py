import unittest

from orangewidget.tests.base import WidgetTest

from AnyQt.QtCore import Qt, QPoint
from AnyQt.QtTest import QTest, QSignalSpy
from AnyQt.QtWidgets import QTreeWidget, QTreeWidgetItem, QTreeWidgetItemIterator

from Orange.data import Table, Domain, StringVariable, DiscreteVariable
from Orange.widgets.widget import Msg, OWWidget
from Orange.widgets.settings import SettingProvider
from Orange.widgets.tests.utils import simulate

from orangecontrib.bioinformatics.widgets.components import GeneSetSelection
from orangecontrib.bioinformatics.widgets.utils.data import TableAnnotation


class MockWidget(OWWidget):
    name = 'Mock'
    want_main_area = False
    selection_component = SettingProvider(GeneSetSelection)

    class Error(OWWidget.Error):
        custom_gene_sets_table_format = Msg('FooBar')

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

        self.selection_component = GeneSetSelection(self, self.controlArea)


class TestGeneSetSelection(WidgetTest):
    def setUp(self):
        self.widget = MockWidget()
        self.component = self.widget.selection_component

        self.mouse_tax_id: str = '10090'
        self.human_tax_id: str = '9606'
        self.dicty_tax_id: str = '44689'

        self.default_selection = [('GO', 'biological_process')]
        self.marker_genes = [('Marker Genes', 'Panglao'), ('Marker Genes', 'CellMarker')]
        self.gene_ontologies = [
            ('GO', 'cellular_component'),
            ('GO', 'biological_process'),
            ('GO', 'molecular_function'),
        ]

    def __item_iterator(self):
        iterator = QTreeWidgetItemIterator(self.component.hierarchy_tree_widget, QTreeWidgetItemIterator.All)

        while iterator.value():
            item: QTreeWidgetItem = iterator.value()
            if not item.childCount():
                yield item
            iterator += 1

    def __initialization(self, tax_id: str) -> None:
        self.component.initialize(tax_id)

        # test if gene sets are shown
        gene_set_items = {item.hierarchy for item in self.__item_iterator()}
        self.assertTrue(len(gene_set_items))
        self.assertEqual(self.component.gene_sets.hierarchies(), gene_set_items)

        # we expect gene ontologies
        for go in self.gene_ontologies:
            self.assertIn(go, gene_set_items)

        # test default selection
        selection = self.component._get_selection()
        self.assertTrue(len(selection))
        self.assertEqual(selection, self.default_selection)

    def test_dictyostelium(self):
        self.__initialization(self.dicty_tax_id)
        self.assertEqual(self.component.gene_sets.common_org(), self.dicty_tax_id)

    def test_mus_musculus(self):
        self.__initialization(self.mouse_tax_id)
        self.assertEqual(self.component.gene_sets.common_org(), self.mouse_tax_id)

    def test_homo_sapiens(self):
        self.__initialization(self.human_tax_id)
        self.assertEqual(self.component.gene_sets.common_org(), self.human_tax_id)

    def test_selection(self):
        self.__initialization(self.human_tax_id)

        panglao, _ = self.marker_genes

        new_selection = None
        for item in self.__item_iterator():
            if item.hierarchy == panglao:
                new_selection = item

        tree_widget: QTreeWidget = self.component.hierarchy_tree_widget
        # obtain the rectangular coordinates
        rect = tree_widget.visualItemRect(new_selection)
        emitted_signals: QSignalSpy = QSignalSpy(self.component.selection_changed)

        # simulate click on item but not in checkbox area.
        # checkbox status is not changed and signal is not emitted
        click_item_text = QPoint(rect.x() + 100, rect.y() + 5)
        QTest.mouseClick(tree_widget.viewport(), Qt.LeftButton, Qt.NoModifier, click_item_text)
        self.assertFalse(len(emitted_signals))
        self.assertEqual(self.component._get_selection(), self.default_selection)

        # simulate click on items checkbox.
        # checkbox status is changed and signal is emitted
        click_checkbox = QPoint(rect.x() + 10, rect.y() + 5)
        QTest.mouseClick(tree_widget.viewport(), Qt.LeftButton, Qt.NoModifier, click_checkbox)
        self.assertTrue(len(emitted_signals))
        self.assertNotEqual(self.component._get_selection(), self.default_selection)
        self.assertEqual(set(self.component._get_selection()), set(self.default_selection + [new_selection.hierarchy]))

    def test_custom_gene_sets(self):
        domain = Domain(
            [], metas=[StringVariable('genes')], class_vars=DiscreteVariable('cell_type', values=['foo', 'bar'])
        )
        data = Table.from_list(domain, [[0, 'gene1'], [1, 'gene2'], [1, 'gene3']])
        data.name = 'custom_gene_sets'
        data.attributes[TableAnnotation.tax_id] = self.human_tax_id
        data.attributes[TableAnnotation.gene_as_attr_name] = False
        data.attributes[TableAnnotation.gene_id_column] = 'genes'

        self.__initialization(self.human_tax_id)
        self.component.initialize_custom_gene_sets(data)

        self.assertTrue(self.component.gene_sets_combobox.isEnabled())
        self.assertEqual(self.component.custom_gene_set_hierarchy, (data.name,))
        self.assertEqual(self.component._get_selection(), [(data.name,)])
        self.assertIn((data.name,), self.component.gene_sets.map_hierarchy_to_sets().keys())

        emitted_signals: QSignalSpy = QSignalSpy(self.component.selection_changed)
        simulate.combobox_run_through_all(self.component.gene_sets_combobox)
        self.assertEqual(len(emitted_signals), len(data.domain.metas) + len(data.domain.class_vars))

        self.component.initialize_custom_gene_sets(None)
        self.assertNotIn((data.name,), self.component.gene_sets.map_hierarchy_to_sets().keys())
        self.assertIsNone(self.component.custom_gene_set_hierarchy)
        self.assertTrue(self.component.custom_gs_box.isHidden())


if __name__ == "__main__":
    unittest.main()
