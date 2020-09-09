from typing import Set, List, Tuple, Optional, DefaultDict
from collections import defaultdict

from AnyQt.QtCore import Qt, QObject, pyqtSignal
from AnyQt.QtWidgets import QComboBox, QTreeView, QTreeWidget, QTreeWidgetItem, QTreeWidgetItemIterator

from Orange.data import Table, StringVariable, DiscreteVariable
from Orange.widgets.gui import comboBox, widgetBox
from Orange.widgets.widget import OWComponent
from Orange.widgets.settings import Setting
from Orange.widgets.utils.itemmodels import DomainModel

from orangecontrib.bioinformatics import geneset
from orangecontrib.bioinformatics.geneset import GeneSet, GeneSets, load_gene_sets
from orangecontrib.bioinformatics.widgets.utils.data import TableAnnotation


class GeneSetSelection(OWComponent, QObject):
    selection_changed = pyqtSignal()

    # component settings
    custom_gene_set_indicator: int
    custom_gene_set_indicator = Setting(None, schema_only=True)

    selection: list
    selection = Setting([('GO', 'biological_process')], schema_only=True)

    def __init__(self, parent_widget, parent_component):
        QObject.__init__(self)
        OWComponent.__init__(self, widget=parent_widget)
        self.parent_widget = parent_widget
        self.parent_component = parent_component

        # gene sets
        self._gene_sets: GeneSets = GeneSets()

        self.hierarchy_tree_widget: QTreeWidget = QTreeWidget()
        self.hierarchy_tree_widget.setHeaderHidden(True)
        self.hierarchy_tree_widget.setEditTriggers(QTreeView.NoEditTriggers)
        self.hierarchy_tree_widget.itemClicked.connect(self._on_item_clicked)
        # self.hierarchy_tree_widget.setStyleSheet("QTreeWidget::item { border-bottom: 1px solid black;}")

        box = widgetBox(parent_component, 'Gene Sets', margin=0)
        box.setFlat(True)
        box.layout().addWidget(self.hierarchy_tree_widget)

        # custom gene sets
        self.domain_model: DomainModel = DomainModel(valid_types=(DiscreteVariable, StringVariable))
        self.data: Optional[Table] = None
        self.custom_gene_set_hierarchy: Optional[str] = None
        self.gene_sets_combobox: Optional[QComboBox] = None
        self.custom_gs_box = widgetBox(self.parent_component, 'Custom Gene Sets', margin=0)
        self.custom_gs_box.setHidden(True)

    @property
    def gene_sets(self):
        return self._gene_sets

    @property
    def tax_id(self) -> Optional[str]:
        if self.data:
            return self.data.attributes[TableAnnotation.tax_id]

    @property
    def gene_as_attr_name(self) -> Optional[str]:
        if self.data:
            return self.data.attributes[TableAnnotation.gene_as_attr_name]

    @property
    def gene_ids_location(self) -> Optional[str]:
        """
        Return a column where gene ids are located.
        """
        if not self.data:
            return

        if self.gene_as_attr_name:
            return

        return self.data.attributes[TableAnnotation.gene_id_column]

    @property
    def num_of_custom_sets(self) -> int:
        custom_gene_sets = self._gene_sets.map_hierarchy_to_sets().get(self.custom_gene_set_hierarchy, None)
        return len(custom_gene_sets) if custom_gene_sets else 0

    @property
    def num_of_genes(self) -> int:
        return self.data.X.shape[0]

    def _on_item_clicked(self) -> None:
        """
        The problem is that signal 'itemClicked' activates even if we don't click on checkbox area.
        That is why we must check if status of checkboxes changed and only then notify widget of the change.
        """
        _selection = self._get_selection()
        if set(_selection) == set(self.selection):
            return

        self.selection = _selection
        self.selection_changed.emit()

    def _on_custom_gene_set_indicator_changed(self):
        self._gene_sets.delete_sets_by_hierarchy(self.custom_gene_set_hierarchy)
        self._load_custom_gene_sets()
        self._update_tree_widget()

    def _update_tree_widget(self) -> None:
        self._set_tree_widget_model(self.hierarchy_tree_widget, self._gene_sets.hierarchies())
        self._set_selection()

    def _load_custom_gene_sets(self) -> None:
        genes = self._extract_gene_sets_from_data()
        if genes:
            self._gene_sets.update(genes)

    def _load_gene_sets(self, tax_id: str) -> None:
        self._gene_sets = GeneSets()

        for gene_set in geneset.list_all(organism=tax_id):
            self._gene_sets.update(list(load_gene_sets(gene_set, tax_id)))

    def _set_tree_widget_model(self, widget: QTreeWidget, hierarchies: Set[Tuple[str]]) -> None:
        widget.clear()
        tree_widget_model: DefaultDict[str, list] = defaultdict(list)

        for hierarchy in hierarchies:
            domain: str = hierarchy[0]
            sub_domain: str = hierarchy[1] if len(hierarchy) == 2 else None
            tree_widget_model[domain].append(sub_domain)

        for domain, sub_domains in sorted(tree_widget_model.items(), key=lambda x: x[0].lower()):
            parent = QTreeWidgetItem(widget, [domain])
            parent.setFlags(parent.flags() & (Qt.ItemIsUserCheckable | Qt.ItemIsEnabled))
            parent.setCheckState(0, Qt.Unchecked)
            parent.setExpanded(True)
            parent.hierarchy = (domain,)
            parent.setFlags(parent.flags() | Qt.ItemIsTristate)

            for sub_domain in sub_domains:
                if sub_domain is None:
                    continue

                child = QTreeWidgetItem(parent, [sub_domain.replace('_', ' ').title()])
                child.setFlags(child.flags() & (Qt.ItemIsUserCheckable | Qt.ItemIsEnabled))
                child.setCheckState(0, Qt.Unchecked)
                child.hierarchy = (domain, sub_domain)

    def _set_selection(self) -> None:
        iterator = QTreeWidgetItemIterator(self.hierarchy_tree_widget, flags=QTreeWidgetItemIterator.All)

        while iterator.value():
            if iterator.value().hierarchy in self.selection:
                item: QTreeWidgetItem = iterator.value()
                item.setCheckState(0, Qt.Checked)

            iterator += 1

        self.selection_changed.emit()

    def _get_selection(self) -> List[Tuple[str, ...]]:
        """ Return a list of selected hierarchies """
        iterator = QTreeWidgetItemIterator(self.hierarchy_tree_widget, flags=QTreeWidgetItemIterator.Checked)

        selected: List[Tuple[str, ...]] = []
        while iterator.value():
            item: QTreeWidgetItem = iterator.value()

            if not item.childCount():
                selected.append(iterator.value().hierarchy)

            iterator += 1
        return selected

    def _extract_gene_sets_from_data(self) -> Optional[List[GeneSet]]:
        if not self.custom_gene_set_indicator:
            return

        if self.data is None or self.gene_ids_location is None:
            return

        if isinstance(self.custom_gene_set_indicator, DiscreteVariable):
            labels = self.custom_gene_set_indicator.values
            gene_sets_names = [
                labels[int(idx)] for idx in self.data.get_column_view(self.custom_gene_set_indicator)[0]
            ]
        else:
            gene_sets_names, _ = self.data.get_column_view(self.custom_gene_set_indicator)

        gene_names, _ = self.data.get_column_view(self.gene_ids_location)
        domain = (self.data.name if self.data.name else 'Custom sets',)
        self.custom_gene_set_hierarchy = domain

        self.selection = [self.custom_gene_set_hierarchy]

        temp_dict = defaultdict(list)
        for set_name, gene_name in zip(gene_sets_names, gene_names):
            temp_dict[set_name].append(gene_name)

        return [GeneSet(gs_id=key, hierarchy=domain, name=key, genes=set(value)) for key, value in temp_dict.items()]

    def initialize_custom_gene_sets(self, data: Optional[Table]) -> None:
        self.data = None
        self.domain_model.set_domain(None)
        self.parent_widget.Error.custom_gene_sets_table_format.clear()

        if data:
            self.data = data
            if self.gene_as_attr_name:
                self.parent_widget.Error.custom_gene_sets_table_format()
                self.data = None
                return

            self.domain_model.set_domain(self.data.domain)
            if not self.gene_sets_combobox:
                self.gene_sets_combobox = comboBox(
                    self.custom_gs_box,
                    self,
                    'custom_gene_set_indicator',
                    sendSelectedValue=True,
                    model=self.domain_model,
                    callback=self._on_custom_gene_set_indicator_changed,
                )
            self.custom_gs_box.setFlat(True)
            self.custom_gs_box.setHidden(False)
            self.custom_gs_box.layout().addWidget(self.gene_sets_combobox)

            if self.custom_gene_set_indicator in self.domain_model:
                index = self.domain_model.indexOf(self.custom_gene_set_indicator)
                self.custom_gene_set_indicator = self.domain_model[index]
            else:
                self.custom_gene_set_indicator = self.domain_model[0]

            self._load_custom_gene_sets()
            self._update_tree_widget()
        else:
            self.gene_sets.delete_sets_by_hierarchy(self.custom_gene_set_hierarchy)
            self._update_tree_widget()
            self.custom_gs_box.setHidden(True)
            self.custom_gene_set_hierarchy = None

    def initialize(self, tax_id: str):
        """Initialize widget component with provided data."""
        self._load_gene_sets(tax_id)
        if self.data:
            self._load_custom_gene_sets()
        self._update_tree_widget()
