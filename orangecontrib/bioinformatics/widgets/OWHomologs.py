""" OWMarkerGenes """
import sys
from typing import List, Union, Optional

from AnyQt.QtCore import QSize

from Orange.data import Table, Domain, StringVariable
from Orange.widgets import gui, widget
from Orange.widgets.widget import Msg
from Orange.widgets.settings import Setting

from orangecontrib.bioinformatics.ncbi.gene import Gene, GeneMatcher, load_gene_summary
from orangecontrib.bioinformatics.ncbi.taxonomy import (
    COMMON_NAMES_MAPPING,
    common_taxid_to_name,
    species_name_to_taxid,
)
from orangecontrib.bioinformatics.widgets.utils.data import (
    ERROR_ON_MISSING_TAX_ID,
    ERROR_ON_MISSING_GENE_ID,
    ERROR_ON_MISSING_ANNOTATION,
    TableAnnotation,
)

HOMOLOG_SYMBOL = "Homolog"
HOMOLOG_ID = "Homolog ID"


class OWHomologs(widget.OWWidget):
    name = "Homologs"
    icon = 'icons/OWHomologs.svg'
    priority = 120

    class Outputs:
        genes = widget.Output("Genes", Table)

    class Inputs:
        data = widget.Input("Data", Table)

    class Warning(widget.OWWidget.Warning):
        no_genes = Msg("Missing data on input.")
        missing_tax_id = Msg(ERROR_ON_MISSING_TAX_ID)
        mising_gene_as_attribute_name = Msg(ERROR_ON_MISSING_ANNOTATION)
        missing_gene_id = Msg(ERROR_ON_MISSING_GENE_ID)
        mising_gene_id_attribute = Msg(ERROR_ON_MISSING_ANNOTATION)

    want_main_area = False

    auto_commit = Setting(True)
    selected_organism: str = Setting('')

    def __init__(self):
        super().__init__()
        self.taxonomy_names: List[str] = list(COMMON_NAMES_MAPPING.values())
        self.taxonomy_ids: List[str] = list(COMMON_NAMES_MAPPING.keys())
        self.source_tax: Optional[str] = None
        self.target_tax: Optional[str] = None
        self.data: Optional[Table] = None

        info_box = gui.vBox(self.controlArea, "Info")
        self.info_gene_type = gui.widgetLabel(info_box, 'No data on input.')
        self.info_gene_type.setWordWrap(True)
        self.info_gene = gui.widgetLabel(info_box, ' ')
        self.info_gene.setWordWrap(True)
        info_box.setMinimumWidth(200)
        gui.separator(self.controlArea)

        self.combo_box_id = -1
        self.target_organism = gui.comboBox(self.controlArea, self, 'combo_box_id')
        self.target_organism.addItems(self.taxonomy_names)
        self.target_organism.activated[int].connect(self.target_organism_change)

        gui.auto_commit(self.controlArea, self, "auto_commit", "Commit", "Commit Automatically")

        self.info.set_input_summary("0")
        self.info.set_output_summary("0")

    @Inputs.data
    def set_data(self, data: Table) -> None:
        self.Warning.clear()
        self.data = data

        if self.data:
            if TableAnnotation.gene_as_attr_name not in self.data.attributes:
                self.Warning.mising_gene_as_attribute_name()
                self.data = None
                return
            if self.data.attributes[TableAnnotation.gene_as_attr_name]:
                if TableAnnotation.gene_id_attribute not in self.data.attributes:
                    self.Warning.mising_gene_id_attribute()
                    self.data = None
                    return

            else:
                if TableAnnotation.tax_id not in self.data.attributes:
                    self.Warning.missing_tax_id()
                    self.data = None
                    return
                if TableAnnotation.gene_id_column not in self.data.attributes:
                    self.Warning.mising_gene_as_attribute_name()
                    self.data = None
                    return
                if self.data.attributes[TableAnnotation.gene_id_column] not in self.data.domain:
                    self.Warning.missing_gene_id()
                    self.data = None
                    return
        else:
            self.info.set_input_summary("0")
            self.info.set_output_summary("0")
            self.info_gene.clear()
            self.info_gene_type.setText("No data on input.")
            self.Outputs.genes.send(None)

            return

        self.source_tax = data.attributes[TableAnnotation.tax_id]
        taxonomy = common_taxid_to_name(self.source_tax)
        self.target_organism.clear()
        self.target_organism.addItems([tax_name for tax_name in self.taxonomy_names if tax_name != taxonomy])

        if taxonomy == self.selected_organism:
            self.combo_box_id = -1
            self.selected_organism = self.taxonomy_names[0]
            self.target_tax = species_name_to_taxid(self.selected_organism)
        else:
            try:
                self.combo_box_id = self.taxonomy_names.index(self.selected_organism)
            except ValueError:
                self.combo_box_id = -1

            if self.combo_box_id != -1:
                self.target_organism.setCurrentIndex(self.combo_box_id)
                self.selected_organism = self.taxonomy_names[self.combo_box_id]
                self.target_tax = species_name_to_taxid(self.selected_organism)
            else:
                self.target_organism.setCurrentIndex(0)
                self.selected_organism = self.taxonomy_names[0]
                self.target_tax = species_name_to_taxid(self.selected_organism)

        self.info_gene_type.setText(f"Organism: {taxonomy}")
        data_len = (
            len(data.domain.attributes) if self.data.attributes[TableAnnotation.gene_as_attr_name] else len(data)
        )
        self.info_gene.setText(f"Number of genes: {data_len}")
        self.info.set_input_summary(f"{data_len}")

        self.commit()

    def find_homologs(self, genes: List[Union[str, Gene]]) -> List[Optional[Gene]]:
        gm = GeneMatcher(self.source_tax)
        gm.genes = genes

        homologs = [g.homolog_gene(taxonomy_id=self.target_tax) for g in gm.genes]
        homologs = load_gene_summary(self.target_tax, homologs)

        return homologs

    def target_organism_change(self, combo_box_id: int) -> None:
        self.combo_box_id = combo_box_id
        self.selected_organism = self.target_organism.itemText(combo_box_id)
        self.target_tax = species_name_to_taxid(self.selected_organism)

        self.commit()

    def commit(self):
        if self.data:
            if self.data.attributes[TableAnnotation.gene_as_attr_name]:
                domain = self.data.domain.copy()
                table = self.data.transform(domain)

                gene_loc = table.attributes[TableAnnotation.gene_id_attribute]
                genes = [str(attr.attributes.get(gene_loc, None)) for attr in table.domain.attributes]
                homologs = self.find_homologs(genes)

                for homolog, col in zip(homologs, table.domain.attributes):
                    if homolog:
                        col.attributes[HOMOLOG_SYMBOL] = homolog.symbol
                        col.attributes[HOMOLOG_ID] = homolog.gene_id

                table = table.from_table(
                    Domain(
                        [col for col in table.domain.attributes if HOMOLOG_ID in col.attributes],
                        table.domain.class_vars,
                        table.domain.metas,
                    ),
                    table,
                )
                out_table = table if len(table.domain.attributes) > 0 else None
            else:
                genes, _ = self.data.get_column_view(self.data.attributes[TableAnnotation.gene_id_column])

                homologs = self.find_homologs(genes)
                homolog = StringVariable(HOMOLOG_SYMBOL)
                homolog_id = StringVariable(HOMOLOG_ID)
                domain = Domain(
                    self.data.domain.attributes,
                    self.data.domain.class_vars,
                    self.data.domain.metas + (homolog, homolog_id),
                )

                table = self.data.transform(domain)
                col, _ = table.get_column_view(homolog)
                col[:] = [g.symbol if g else "?" for g in homologs]
                col, _ = table.get_column_view(homolog_id)
                col[:] = [g.gene_id if g else "?" for g in homologs]

                # note: filter out rows with unknown homologs
                table = table[table.get_column_view(homolog_id)[0] != "?"]

                out_table = table if len(table) > 0 else None

            self.info.set_output_summary(f"{len(out_table) if out_table else 0}")

            self.Outputs.genes.send(out_table)
        else:
            self.Outputs.genes.send(None)

    def closeEvent(self, event):
        super().closeEvent(event)

    def sizeHint(self):
        return super().sizeHint().expandedTo(QSize(200, 200))


def main(argv=None):
    from AnyQt.QtWidgets import QApplication

    app = QApplication(argv or sys.argv)
    w = OWHomologs()
    w.show()
    w.activateWindow()
    rv = app.exec_()
    w.saveSettings()
    w.onDeleteWidget()
    return rv


if __name__ == "__main__":
    sys.exit(main())
