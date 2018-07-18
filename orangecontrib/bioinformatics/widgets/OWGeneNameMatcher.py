""" GeneNameMatching """
import os
import threading
import numpy as np
import re

from typing import Set

from AnyQt.QtWidgets import (
    QSplitter, QTableView, QWidget, QVBoxLayout, QItemDelegate, QStyledItemDelegate, QHeaderView, QStyleOptionViewItem,
    QStyle, QAbstractItemView, QApplication
)
from AnyQt.QtCore import (
    Qt, QSize, QThreadPool, QSortFilterProxyModel, QAbstractTableModel, QVariant,

)
from AnyQt.QtGui import (
    QIcon, QFont, QAbstractTextDocumentLayout, QFontMetrics, QTextDocument
)

from Orange.widgets.gui import (
    vBox, comboBox, ProgressBar, widgetBox, auto_commit, widgetLabel, checkBox,
    rubber, radioButtons, separator, hBox
)
from Orange.widgets.widget import OWWidget
from Orange.widgets.utils import itemmodels
from Orange.widgets.settings import Setting
from Orange.widgets.utils.signals import Output, Input
from Orange.data import StringVariable, Domain, Table, filter as table_filter

from orangecontrib.bioinformatics.widgets.utils.gui import horizontal_line

from orangecontrib.bioinformatics.widgets.utils.data import (
    TAX_ID, GENE_AS_ATTRIBUTE_NAME, GENE_ID_COLUMN, GENE_ID_ATTRIBUTE
)
from orangecontrib.bioinformatics.widgets.utils.concurrent import Worker
from orangecontrib.bioinformatics.ncbi import taxonomy
from orangecontrib.bioinformatics.ncbi.gene import GeneMatcher, NCBI_ID, GENE_MATCHER_HEADER


class GeneInfoModel(QAbstractTableModel):

    def __init__(self):
        QAbstractTableModel.__init__(self)
        self.__items = np.array([])
        self.__data_matrix = np.array([])

        self.header_labels, self.header_tags = GENE_MATCHER_HEADER

    @property
    def model_items(self):
        return self.__items

    @model_items.setter
    def model_items(self, gene_objects):
        self.__items = np.array(gene_objects)

        for gene in gene_objects:
            # load info from database
            gene.load_ncbi_info()
            # populate data matrix
            if self.__data_matrix.size == 0:
                self.__data_matrix = np.append(self.__data_matrix, self.__gene_object_to_list(gene))
            else:
                self.__data_matrix = np.vstack((self.__data_matrix, self.__gene_object_to_list(gene)))

    def __gene_object_to_list(self, gene_object):
        output_list = []

        for tag in self.header_tags:
            gene_attr = gene_object.__getattribute__(tag)

            if isinstance(gene_attr, dict):
                # note: db_refs are stored as dicts
                gene_attr = ', '.join('{}: {}'.format(key, val)
                                      for (key, val) in gene_attr.items()) if gene_attr else ' '
            elif isinstance(gene_attr, list):
                # note: synonyms are stored as lists
                gene_attr = ', '.join(gene_attr) if gene_attr else ' '

            output_list.append(gene_attr)

        return output_list

    def rowCount(self, *args, **kwargs):
        return self.__data_matrix.shape[0]

    def columnCount(self, *args, **kwargs):
        return self.__data_matrix.shape[1]

    def headerData(self, section, orientation, role=Qt.DisplayRole):
        if role == Qt.DisplayRole and orientation == Qt.Horizontal:
            return self.header_labels[section]
        return QAbstractTableModel.headerData(self, section, orientation, role)

    def data(self, model_index, role=None):
        # check if data is set
        if self.__data_matrix.size == 0 or len(self.__items) == 0:
            return QVariant()

        # return empty QVariant if model index is unknown
        if not model_index.isValid() \
                or not (0 <= model_index.row() < self.rowCount()) \
                or not (0 <= model_index.column() < self.columnCount()):
            return QVariant()

        if role == Qt.DisplayRole:
            # note: Data is not displayed if QVariant is returned, why?
            #       return QVariant(self.__data_matrix[model_index.row()][model_index.column()])
            return '{}'.format(self.__data_matrix[model_index.row()][model_index.column()])


class OWGeneNameMatcher(OWWidget):
    name = "Gene Name Matcher"
    description = "Tool for working with genes"
    icon = "../widgets/icons/OWGeneInfo.svg"
    priority = 5
    want_main_area = True

    use_attr_names = Setting(True)
    selected_organism = Setting(11)

    selected_filter = Setting(0)
    gene_as_attr_name = Setting(0)
    filter_unknown = Setting(True)
    include_entrez_id = Setting(True)
    # include_ensembl_id = Setting(True)
    auto_commit = Setting(True)

    class Inputs:
        data_table = Input("Data", Table)

    class Outputs:
        custom_data_table = Output("Data", Table)

    class Information(OWWidget.Information):
        pass

    def sizeHint(self):
        return QSize(1280, 960)

    def __init__(self):
        super().__init__()
        # ATTRIBUTES #

        # input data
        self.input_data = None
        self.input_genes = None
        self.tax_id = None
        self.column_candidates = []
        self.selected_gene_col = None

        # input options
        self.organisms = []

        # gene matcher
        self.gene_matcher = None

        # threads
        self.threadpool = QThreadPool(self)
        self.workers = None

        # progress bar
        self.progress_bar = None

        # filter
        self.filter_labels = ['Unique', 'Partial', 'Unknown']

        # GUI SECTION #

        # Control area
        self.info_box = widgetLabel(
            widgetBox(self.controlArea, "Info", addSpace=True), "Initializing\n"
        )

        organism_box = vBox(self.controlArea, 'Organism')
        self.organism_select_combobox = comboBox(organism_box, self,
                                                 'selected_organism',
                                                 callback=self.on_input_option_change)

        self.get_available_organisms()
        self.organism_select_combobox.setCurrentIndex(self.selected_organism)

        box = widgetBox(self.controlArea, 'Gene names')
        self.gene_columns_model = itemmodels.DomainModel(valid_types=(StringVariable, ))
        self.gene_column_combobox = comboBox(box, self, 'selected_gene_col',
                                             model=self.gene_columns_model,
                                             sendSelectedValue=True,
                                             callback=self.on_input_option_change)

        self.attr_names_checkbox = checkBox(box, self, 'use_attr_names', 'Use attribute names',
                                            disables=[(-1, self.gene_column_combobox)],
                                            callback=self.on_input_option_change)

        self.gene_column_combobox.setDisabled(bool(self.use_attr_names))

        output_box = vBox(self.controlArea, 'Output settings')

        # TODO: will widget support transposing tables?
        # radioButtonsInBox(output_box, self, "gene_as_attr_name", ["Genes in rows", "Genes in columns"],
        # callback=self.on_output_option_change)

        # separator(output_box)
        checkBox(output_box, self, 'filter_unknown', 'Filter unknown genes', callback=self.on_output_option_change)
        separator(output_box)
        output_box.layout().addWidget(horizontal_line())
        checkBox(output_box, self, 'include_entrez_id', 'Include Entrez ID', callback=self.on_output_option_change)

        # TODO: provide support for ensembl ids as output option
        # checkBox(output_box, self, 'include_ensembl_id', 'Include Ensembl ID', callback=self.on_output_option_change)

        auto_commit(self.controlArea, self, "auto_commit", label="Commit")

        rubber(self.controlArea)

        # Main area
        filter_box = hBox(self.mainArea, 'Filter results')
        self.radio_group = radioButtons(filter_box, self, value='selected_filter',
                                        btnLabels=self.filter_labels,
                                        orientation=Qt.Horizontal, callback=self.on_filter_changed)
        rubber(self.radio_group)
        self.mainArea.layout().addWidget(filter_box)

        self.proxy_model = QSortFilterProxyModel()
        # left side list view
        self.table_view = QTableView()
        self.table_view.setModel(self.proxy_model)
        self.table_view.horizontalHeader().setStretchLastSection(True)
        # self.table_view.horizontalHeader().setSectionResizeMode(QHeaderView.Stretch)
        # self.table_view.selectionModel().selectionChanged.connect(self.__selection_changed)

        self.mainArea.layout().addWidget(self.table_view, 1)

    def __reset_widget_state(self):
        self.Outputs.custom_data_table.send(None)
        # self.proxy_model.setSourceModel(None)
        # self.extended_view.reset_genes_model()
        # self.extended_view.reset_info_model()

    def __selection_changed(self):
        genes = [model_index.data() for model_index in self.extended_view.get_selected_gens()]
        self.extended_view.set_info_model(genes)

    def _update_info_box(self):

        if self.input_genes and self.gene_matcher:
            num_genes = len(self.gene_matcher.genes)
            known_genes = len(self.gene_matcher.get_known_genes())

            info_text = 'Genes on input:  {}\n' \
                        'Known genes :    {} ({:.2f} %)\n'.format(num_genes, known_genes, known_genes * 100 / num_genes)

        else:
            info_text = 'No genes on input'

        self.info_box.setText(info_text)

    def _progress_advance(self):
        # GUI should be updated in main thread. That's why we are calling advance method here
        if self.progress_bar:
            self.progress_bar.advance()

    def _handle_matcher_results(self):
        assert threading.current_thread() == threading.main_thread()

        if self.progress_bar:
            self.progress_bar.finish()
            self.setStatusMessage('')

        # if no known genes, clean up and return
        if not len(self.gene_matcher.get_known_genes()):
            self._update_info_box()
            self.__reset_widget_state()
            return

        self._update_info_box()
        self.table_model = GeneInfoModel()
        self.table_model.model_items = self.gene_matcher.genes
        self.proxy_model.setSourceModel(self.table_model)
        # self.table_view.resizeRowsToContents()
        # self.commit()

    def get_available_organisms(self):
        available_organism = sorted([(tax_id, taxonomy.name(tax_id)) for tax_id in taxonomy.common_taxids()],
                                    key=lambda x: x[1])

        self.organisms = [tax_id[0] for tax_id in available_organism]
        self.organism_select_combobox.addItems([tax_id[1] for tax_id in available_organism])

    def gene_names_from_table(self):
        """ Extract and return gene names from `Orange.data.Table`.
        """
        self.input_genes = []
        if self.input_data:
            if self.use_attr_names:
                self.input_genes = [str(attr.name).strip() for attr in self.input_data.domain.attributes]
            elif self.selected_gene_col:
                if self.selected_gene_col in self.input_data.domain:
                    self.input_genes = [str(e[self.selected_gene_col]) for e in self.input_data
                                        if not np.isnan(e[self.selected_gene_col])]

    def _update_gene_matcher(self):
        self.gene_names_from_table()

        if not self.input_genes:
            self._update_info_box()

        if not self.gene_matcher:
            self.gene_matcher = GeneMatcher(self.get_selected_organism(), case_insensitive=True)

        self.gene_matcher.genes = self.input_genes
        self.gene_matcher.organism = self.get_selected_organism()

    def get_selected_organism(self):
        return self.organisms[self.selected_organism]

    def match_genes(self):
        if self.gene_matcher:
            # init progress bar
            self.progress_bar = ProgressBar(self, iterations=len(self.gene_matcher.genes))
            # status message
            self.setStatusMessage('Gene matcher running')

            worker = Worker(self.gene_matcher.run_matcher, progress_callback=True)
            worker.signals.progress.connect(self._progress_advance)
            worker.signals.finished.connect(self._handle_matcher_results)

            # move download process to worker thread
            self.threadpool.start(worker)

    def on_input_option_change(self):
        self.__reset_widget_state()
        self._update_gene_matcher()
        self.match_genes()

    @Inputs.data_table
    def handle_input(self, data):
        # self.__reset_widget_state()
        self.gene_columns_model.set_domain(None)

        if data:
            self.input_data = data

            self.gene_columns_model.set_domain(self.input_data.domain)

            if self.gene_columns_model:
                self.selected_gene_col = self.gene_columns_model[0]

            self.tax_id = str(self.input_data.attributes.get(TAX_ID, ''))
            self.use_attr_names = self.input_data.attributes.get(GENE_AS_ATTRIBUTE_NAME, self.use_attr_names)

            if self.tax_id in self.organisms:
                self.selected_organism = self.organisms.index(self.tax_id)

            self.on_input_option_change()

    @staticmethod
    def get_gene_id_identifier(gene_id_strings):
        # type: (Set[str]) -> str

        if not len(gene_id_strings):
            return NCBI_ID

        regex = re.compile(r'Entrez ID \(.*?\)')
        filtered = filter(regex.search, gene_id_strings)

        return NCBI_ID + ' ({})'.format(len(set(filtered)) + 1)

    def __handle_ids(self, data_table):
        """
        If 'use_attr_names' is True, genes from the input data are in columns.
        """
        if self.use_attr_names:
            # set_of_attributes = set([key for attr in data_table.domain[:] for key in attr.attributes.keys()
            # if key.startswith(NCBI_ID)])
            # gene_id = self.get_gene_id_identifier(set_of_attributes)
            gene_id = NCBI_ID

            for gene in self.gene_matcher.genes:
                if gene.ncbi_id:
                    data_table.domain[gene.input_name].attributes[gene_id] = str(gene.ncbi_id)
        else:
            set_of_variables = set([var.name for var in data_table.domain.variables + data_table.domain.metas
                                    if var.name.startswith(NCBI_ID)])

            gene_id = self.get_gene_id_identifier(set_of_variables)

            temp_domain = Domain([], metas=[StringVariable(gene_id)])
            temp_data = [[str(gene.ncbi_id) if gene.ncbi_id else '?'] for gene in self.gene_matcher.genes]
            temp_table = Table(temp_domain, temp_data)

            # if columns differ, then concatenate.
            if NCBI_ID in data_table.domain:
                if gene_id != NCBI_ID and not np.array_equal(np.array(temp_data).ravel(),
                                                             data_table.get_column_view(NCBI_ID)[0]):

                    data_table = Table.concatenate([data_table, temp_table])
                else:
                    gene_id = NCBI_ID
            else:
                data_table = Table.concatenate([data_table, temp_table])

        return data_table, gene_id

    def __apply_filters(self, data_table):
        set_of_attributes = set([key for attr in data_table.domain[:] for key in attr.attributes.keys()
                                 if key == NCBI_ID])

        gene_id = NCBI_ID if NCBI_ID in data_table.domain or set_of_attributes else None

        if self.include_entrez_id:
            data_table, gene_id = self.__handle_ids(data_table)

        if self.filter_unknown:
            known_input_genes = [gene.input_name for gene in self.gene_matcher.get_known_genes()]

            if self.use_attr_names:
                temp_domain = Domain(
                    [attr for attr in data_table.domain.attributes if attr.name in known_input_genes],
                    metas=data_table.domain.metas,
                    class_vars=data_table.domain.class_vars
                )
                data_table = data_table.transform(temp_domain)
            else:

                # create filter from selected column for genes
                only_known = table_filter.FilterStringList(self.selected_gene_col, known_input_genes)
                # apply filter to the data
                data_table = table_filter.Values([only_known])(data_table)

        return data_table, gene_id

    def commit(self):
        self.Outputs.custom_data_table.send(None)

        if not self.input_data:
            return

        if not self.use_attr_names and not self.gene_columns_model:
            return

        output_data_table = self.input_data.transform(self.input_data.domain.copy())
        output_data_table, gene_id = self.__apply_filters(output_data_table.copy())

        # handle table attributes
        output_data_table.attributes[TAX_ID] = self.get_selected_organism()
        output_data_table.attributes[GENE_AS_ATTRIBUTE_NAME] = bool(self.use_attr_names)

        if not bool(self.use_attr_names):
            output_data_table.attributes[GENE_ID_COLUMN] = gene_id
        else:
            output_data_table.attributes[GENE_ID_ATTRIBUTE] = gene_id

        self.Outputs.custom_data_table.send(output_data_table)
        # gene_objs = [self.proxy_model.index(row, 0).data() for row in range(self.proxy_model.rowCount())]

    def on_output_option_change(self):
        self.commit()

    def on_filter_changed(self):
        self.proxy_model.invalidateFilter()
        self.extended_view.genes_view.resizeRowsToContents()
