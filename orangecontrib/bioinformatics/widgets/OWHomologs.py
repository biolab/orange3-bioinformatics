""" OWMarkerGenes """
import sys

from AnyQt.QtCore import QSize

from Orange.widgets import widget, gui
from Orange.widgets.widget import Msg

from Orange.widgets.settings import Setting

from orangecontrib.bioinformatics.widgets.utils.data import TableAnnotation
from orangecontrib.bioinformatics.widgets.utils.data import (
    ERROR_ON_MISSING_ANNOTATION,
    ERROR_ON_MISSING_GENE_ID,
    ERROR_ON_MISSING_TAX_ID,
)

from orangecontrib.bioinformatics.ncbi.taxonomy import COMMON_NAMES
from orangecontrib.bioinformatics.ncbi.gene import GeneMatcher, load_gene_summary

from Orange.data import *


class OWHomologs(widget.OWWidget):
    name = "Homologs"
    icon = 'icons/OWHomologs.svg'
    priority = 170

    class Outputs:
        genes = widget.Output("Genes", Table)

    class Inputs:
        data = widget.Input("Data", Table)

    class Warning(widget.OWWidget.Warning):
        no_genes = Msg(ERROR_ON_MISSING_ANNOTATION)
        missing_tax_id = Msg(ERROR_ON_MISSING_TAX_ID)
        mising_gene_as_attribute_name = Msg(ERROR_ON_MISSING_ANNOTATION)
        missing_gene_id = Msg(ERROR_ON_MISSING_GENE_ID)

    want_main_area = False

    auto_commit = Setting(True)
    selected_organism: str = Setting('')

    def __init__(self):
        super().__init__()
        self.taxonomy_ids = dict(COMMON_NAMES)
        info_box = gui.vBox(self.controlArea, "Info")
        self.info_gene_type = gui.widgetLabel(info_box, 'No data on input.')
        self.info_gene_type.setWordWrap(True)
        self.info_gene = gui.widgetLabel(info_box, ' ')
        self.info_gene.setWordWrap(True)
        info_box.setMinimumWidth(200)
        gui.separator(self.controlArea)
        self.taxonomy = None
        self.data = None
        self.output = None

        self.organism_id = -1
        self.organism = gui.comboBox(self.controlArea, self, 'organism_id')
        self.organism.addItems(self.taxonomy_ids.values())
        self.organism.activated[int].connect(self.organisms)

        gui.auto_commit(self.controlArea, self, "auto_commit", "Commit", "Commit Automatically")

        self.info.set_input_summary(self.info.NoInput)
        self.info.set_output_summary(self.info.NoInput)

    @Inputs.data
    def set_data(self, data):
        self.Warning.no_genes.clear()
        self.taxonomy_ids = dict(COMMON_NAMES)
        self.data = data

        if self.data:
            if TableAnnotation.gene_as_attr_name not in self.data.attributes:
                self.Warning.mising_gene_as_attribute_name()
                return
            if TableAnnotation.tax_id not in self.data.attributes:
                self.Warning.missing_tax_id()
                return
            if TableAnnotation.gene_id_column not in self.data.attributes:
                self.Warning.missing_tax_id_genes()
                return
            if self.data.attributes[TableAnnotation.gene_id_column] not in self.data.domain:
                self.Warning.missing_gene_id()
                return
        else:
            self.Warning.no_genes()
            return

        self.taxonomy = self.taxonomy_ids[data.attributes[TableAnnotation.tax_id]]
        if self.taxonomy == self.selected_organism:
            del self.taxonomy_ids[data.attributes[TableAnnotation.tax_id]]
            self.organism.clear()
            self.organism.addItems(self.taxonomy_ids.values())
        else:
            try:
                self.organism_id = list(self.taxonomy_ids.values()).index(self.selected_organism)
            except ValueError:
                self.organism_id = -1

            del self.taxonomy_ids[data.attributes[TableAnnotation.tax_id]]
            self.organism.clear()
            self.organism.addItems(self.taxonomy_ids.values())

            if self.organism_id != -1:
                self.organism.setCurrentIndex(self.organism_id)
                self.selected_organism = list(self.taxonomy_ids.values())[self.organism_id]

        self.info_gene_type.setText("Organism: " + self.taxonomy)
        self.info_gene.setText("Number of genes: " + str(len(data)))
        self.info.set_input_summary(f"{str(len(data))}")

        self.commit()

    def find_homologs(self, genes, organism):
        gm = GeneMatcher(organism)
        gm.genes = genes

        homologs = [g.homolog_genes(taxonomy_id=self.target_organism) for g in gm.genes]
        homologs = load_gene_summary(self.target_organism, homologs)

        return homologs

    def organisms(self, id):
        self.organism_id = id
        self.selected_organism = self.organism.itemText(id)
        self.commit()

    def commit(self):
        HOMOLOG_SYMBOL = "Homolog"
        HOMOLOG_ID = "Homolog ID"
        out_table = None
        self.target_organism = list(self.taxonomy_ids.keys())[self.organism_id]

        if self.data.attributes[TableAnnotation.gene_as_attr_name]:
            domain = self.data.domain.copy()
            table = self.data.transform(domain)

            gene_loc = table.attributes[TableAnnotation.gene_id_attribute]
            genes = [str(attr.attributes.get(gene_loc, None)) for attr in table.domain.attributes]
            homologs = self.find_homologs(genes, table.attributes[TableAnnotation.tax_id])

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

            homologs = self.find_homologs(genes, self.data.attributes[TableAnnotation.tax_id])
            homolog = StringVariable(HOMOLOG_SYMBOL)
            homolog_id = StringVariable(HOMOLOG_ID)
            domain = Domain(
                self.data.domain.attributes, self.data.domain.class_vars, self.data.domain.metas + (homolog, homolog_id)
            )

            table = self.data.transform(domain)
            col, _ = table.get_column_view(homolog)
            col[:] = [g.symbol if g else "?" for g in homologs]
            col, _ = table.get_column_view(homolog_id)
            col[:] = [g.gene_id if g else "?" for g in homologs]

            # note: filter out rows with unknown homologs
            table = table[table.get_column_view(homolog_id)[0] != "?"]

            out_table = table if len(table) > 0 else None

        if out_table:
            self.info.set_output_summary(f"{len(out_table)}")
        else:
            self.info.set_output_summary(f"{0}")
        self.Outputs.genes.send(out_table)

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
