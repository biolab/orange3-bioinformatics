""" Qt component for Gene sets """
from typing import Union
from collections import defaultdict

import numpy as np

from AnyQt.QtCore import Qt
from AnyQt.QtWidgets import QWidget, QGroupBox, QTreeView, QTreeWidget, QTreeWidgetItem, QTreeWidgetItemIterator

from orangecontrib.bioinformatics.geneset import GeneSet, GeneSets, list_all, load_gene_sets

# TODO: better handle stored selection
# TODO: Don't use hardcoded 'Custom sets', use table name if available


class GeneSetsSelection(QWidget):
    def __init__(self, box, parent, settings_var, **kwargs):
        # type: (Union[QGroupBox, QWidget], QWidget, str) -> None
        super().__init__(**kwargs)

        self.parent = parent
        self.stored_selection = settings_var
        # gene sets object
        self.gs_object = GeneSets()  # type: GeneSets

        self.hierarchy_tree_widget = QTreeWidget(self)
        self.hierarchy_tree_widget.setHeaderHidden(True)
        self.hierarchy_tree_widget.setEditTriggers(QTreeView.NoEditTriggers)
        box.layout().addWidget(self.hierarchy_tree_widget)

        self.custom_set_hier = None
        self.default_selection = [
            ('GO', 'molecular_function'),
            ('GO', 'biological_process'),
            ('GO', 'cellular_component'),
        ]

    def clear_custom_sets(self):
        # delete any custom sets if they exists
        self.gs_object.delete_sets_by_hierarchy(self.custom_set_hier)

    def add_custom_sets(self, gene_sets_names, gene_names, hierarchy_title=None, select_customs_flag=False):
        # type: (np.ndarray, np.ndarray) -> None

        self.custom_set_hier = hierarchy_title
        self.clear_custom_sets()

        temp_dict = defaultdict(list)
        for set_name, gene_name in zip(gene_sets_names, gene_names):
            temp_dict[set_name].append(gene_name)

        g_sets = []
        for key, value in temp_dict.items():
            g_sets.append(
                GeneSet(
                    gs_id=key,
                    hierarchy=self.custom_set_hier,
                    organism=self.gs_object.common_org(),
                    name=key,
                    genes=set(value),
                )
            )

        self.gs_object.update(g_sets)
        self.update_gs_hierarchy(select_customs_flag=select_customs_flag)

    def load_gene_sets(self, tax_id):
        # type: (str) -> None
        self.gs_object = GeneSets()
        self.clear()

        gene_sets = list_all(organism=tax_id)
        self.set_hierarchy_model(self.hierarchy_tree_widget, self.hierarchy_tree(gene_sets))

        for gene_set in gene_sets:
            g_sets = load_gene_sets(gene_set, tax_id)
            self.gs_object.update([g_set for g_set in g_sets])

        self.set_selected_hierarchies()

    def clear_gene_sets(self):
        self.gs_object = GeneSets()

    def clear(self):
        # reset hierarchy widget state
        self.hierarchy_tree_widget.clear()

    def update_gs_hierarchy(self, select_customs_flag=False):
        self.clear()
        self.set_hierarchy_model(self.hierarchy_tree_widget, self.hierarchy_tree(self.gs_object.hierarchies()))
        if select_customs_flag:
            self.set_custom_sets()
        else:
            self.set_selected_hierarchies()

    def set_hierarchy_model(self, tree_widget, sets):
        def beautify_displayed_text(text):
            if '_' in text:
                return text.replace('_', ' ').title()
            else:
                return text

        # TODO: maybe optimize this code?
        for key, value in sets.items():
            item = QTreeWidgetItem(tree_widget, [beautify_displayed_text(key)])
            item.setFlags(item.flags() & (Qt.ItemIsUserCheckable | Qt.ItemIsSelectable | Qt.ItemIsEnabled))
            item.setExpanded(True)
            item.hierarchy = key

            if value:
                item.setFlags(item.flags() | Qt.ItemIsTristate)
                self.set_hierarchy_model(item, value)
            else:
                if item.parent():
                    item.hierarchy = (item.parent().hierarchy, key)

            if not item.childCount() and not item.parent():
                item.hierarchy = (key,)

    def get_hierarchies(self, **kwargs):
        """ return selected hierarchy
        """
        only_selected = kwargs.get('only_selected', None)

        sets_to_display = []

        if only_selected:
            iterator = QTreeWidgetItemIterator(self.hierarchy_tree_widget, QTreeWidgetItemIterator.Checked)
        else:
            iterator = QTreeWidgetItemIterator(self.hierarchy_tree_widget)

        while iterator.value():
            # note: if hierarchy value is not a tuple, then this is just top level qTreeWidgetItem that
            #       holds subcategories. We don't want to display all sets from category
            if type(iterator.value().hierarchy) is not str:

                if not only_selected:
                    sets_to_display.append(iterator.value().hierarchy)
                else:
                    if not iterator.value().isDisabled():
                        sets_to_display.append(iterator.value().hierarchy)

            iterator += 1

        return sets_to_display

    def set_selected_hierarchies(self):
        iterator = QTreeWidgetItemIterator(self.hierarchy_tree_widget, QTreeWidgetItemIterator.All)
        defaults = []

        while iterator.value():

            # note: if hierarchy value is not a tuple, then this is just top level qTreeWidgetItem that
            #       holds subcategories. We don't want to display all sets from category
            if type(iterator.value().hierarchy) is not str:
                if iterator.value().hierarchy in self.parent.__getattribute__(self.stored_selection):
                    iterator.value().setCheckState(0, Qt.Checked)
                else:
                    iterator.value().setCheckState(0, Qt.Unchecked)

            # if no items are checked, set defaults
            if iterator.value().hierarchy in self.default_selection:
                defaults.append(iterator.value())

            iterator += 1

        if len(self.get_hierarchies(only_selected=True)) == 0:
            [item.setCheckState(0, Qt.Checked) for item in defaults]

    def set_custom_sets(self):
        iterator = QTreeWidgetItemIterator(self.hierarchy_tree_widget, QTreeWidgetItemIterator.All)

        while iterator.value():

            # note: if hierarchy value is not a tuple, then this is just top level qTreeWidgetItem that
            #       holds subcategories. We don't want to display all sets from category
            if type(iterator.value().hierarchy) is not str:
                if iterator.value().hierarchy == self.custom_set_hier:
                    iterator.value().setCheckState(0, Qt.Checked)
                else:
                    iterator.value().setCheckState(0, Qt.Unchecked)

            iterator += 1

    @staticmethod
    def hierarchy_tree(gene_sets):
        def tree():
            return defaultdict(tree)

        collection = tree()

        def collect(col, set_hierarchy):
            if set_hierarchy:
                collect(col[set_hierarchy[0]], set_hierarchy[1:])

        for hierarchy in gene_sets:
            collect(collection, hierarchy)

        return collection
