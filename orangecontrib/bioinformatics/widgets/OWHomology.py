""" OWMarkerGenes """
import sys

from AnyQt.QtCore import QSize

from Orange.widgets import widget, gui
from Orange.widgets.widget import Msg

from Orange.widgets.settings import Setting, DomainContextHandler

from orangecontrib.bioinformatics.widgets.utils.data import TAX_ID, \
            GENE_AS_ATTRIBUTE_NAME, GENE_ID_COLUMN, \
    ERROR_ON_MISSING_ANNOTATION, ERROR_ON_MISSING_GENE_ID, ERROR_ON_MISSING_TAX_ID

from orangecontrib.bioinformatics.ncbi.homologene import match_by_rows
from orangecontrib.bioinformatics.ncbi.taxonomy import COMMON_NAMES

from Orange.data import *



class OWHomology(widget.OWWidget):
    name = "Homology"
    icon = 'icons/OWHomologs.svg'
    priority = 170

    class Outputs:
        genes = widget.Output("Genes", Table)

    class Inputs:
        data = widget.Input("Data", Table)

    class Warning(widget.OWWidget.Warning):
        no_genes = Msg("Missing genes table on input.")
        missing_tax_id = Msg(ERROR_ON_MISSING_TAX_ID)
        mising_gene_as_attribute_name = Msg(ERROR_ON_MISSING_ANNOTATION)
        missing_gene_id = Msg(ERROR_ON_MISSING_GENE_ID)
    want_main_area = False

    auto_commit = Setting(True)
    selected_organism: str = Setting('')
    settingsHandler = DomainContextHandler()

    def __init__(self):
        super().__init__()
        self.ids = dict(COMMON_NAMES)
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
        self.organism.addItems(self.ids.values())
        self.organism.activated[int].connect(self.organisms)

        gui.auto_commit(self.controlArea, self, "auto_commit", "Commit",
                        "Commit Automatically")

        self.info.set_input_summary(self.info.NoInput)
        self.info.set_output_summary(self.info.NoInput)

    @Inputs.data
    def set_data(self, data):
        self.Warning.no_genes.clear()
        self.ids = dict(COMMON_NAMES)
        self.data = data

        if self.data:
            if GENE_AS_ATTRIBUTE_NAME not in self.data.attributes:
                self.Warning.mising_gene_as_attribute_name()
            if TAX_ID not in self.data.attributes:
                self.Warning.missing_tax_id()
            if GENE_ID_COLUMN not in self.data.attributes:
                self.Warning.missing_tax_id_genes()
            if self.data.attributes[GENE_ID_COLUMN] not in self.data.domain:
                self.Warning.missing_gene_id()
        else:
            self.Warning.no_genes()
            return

        self.taxonomy = self.ids[data.attributes[TAX_ID]]
        if self.taxonomy == self.selected_organism:
            del self.ids[data.attributes[TAX_ID]]
            self.organism.clear()
            self.organism.addItems(self.ids.values())
        else:
            try:
                self.organism_id = list(self.ids.values()).index(
                                                    self.selected_organism)
            except ValueError:
                self.organism_id = -1

            del self.ids[data.attributes[TAX_ID]]
            self.organism.clear()
            self.organism.addItems(self.ids.values())

            if self.organism_id != -1:
                self.selected_organism = list(self.ids.values()
                                              )[self.organism_id]

        self.info_gene_type.setText("Organism: "+self.taxonomy)
        self.info_gene.setText("Number of genes: "+str(len(data)))
        self.info.set_input_summary(f"{str(len(data))}")

        self.commit()

    def organisms(self, id):
        self.organism_id = id
        self.selected_organism = self.organism.itemText(id)
        self.commit()

    def commit(self):
        self.output = match_by_rows(self.data,
                                    list(self.ids.keys())[self.organism_id])
        self.Outputs.genes.send(self.output)

    def closeEvent(self, event):
        super().closeEvent(event)

    def sizeHint(self):
        return super().sizeHint().expandedTo(QSize(200, 200))


def main(argv=None):
    from AnyQt.QtWidgets import QApplication

    app = QApplication(argv or sys.argv)
    w = OWHomology()
    w.show()
    w.activateWindow()
    rv = app.exec_()
    w.saveSettings()
    w.onDeleteWidget()
    return rv


if __name__ == "__main__":
    sys.exit(main())
