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
    Qt, QSize, QThreadPool, QSortFilterProxyModel, QAbstractListModel, QVariant,

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
from orangecontrib.bioinformatics.ncbi.gene import GeneMatcher, NCBI_ID


class FilterProxyModel(QSortFilterProxyModel):
    def __init__(self, parent=None):
        super(FilterProxyModel, self).__init__(parent)
        self.ow = parent

    def filterAcceptsRow(self, source_row, source_parent):
        model = self.sourceModel()
        if isinstance(model, GeneMatcherModel):
            rows = model.get_rows()
            unique, partial, unknown = range(len(self.ow.filter_labels))

            possible_hits = rows[source_row].possible_hits
            ncbi_id = rows[source_row].ncbi_id

            if self.ow.selected_filter == unique and ncbi_id is not None:
                return True
            elif self.ow.selected_filter == partial and possible_hits:
                return True
            elif self.ow.selected_filter == unknown and ncbi_id is None and not possible_hits:
                return True

        return False


class HTMLDelegate(QStyledItemDelegate):
    """
    https://stackoverflow.com/questions/1956542/how-to-make-item-view-render-rich-html-text-in-qt
    https://stackoverflow.com/questions/2375763/how-to-open-an-url-in-a-qtableview
    """

    def sizeHint(self, option, index):
        options = QStyleOptionViewItem(option)
        gene_obj = index.data(Qt.DisplayRole)
        self.initStyleOption(options, index)

        doc = QTextDocument()
        doc.setHtml(gene_obj.to_html())
        doc.setTextWidth(options.rect.width() - 10)

        return QSize(doc.idealWidth(), doc.size().height())

    def paint(self, painter, option, index):
        options = QStyleOptionViewItem(option)
        gene_obj = index.data(Qt.DisplayRole)
        self.initStyleOption(options, index)

        style = QApplication.style() if options.widget is None else options.widget.style()

        doc = QTextDocument()
        doc.setHtml(gene_obj.to_html())
        doc.setTextWidth(option.rect.width() - 10)

        options.text = ""
        style.drawControl(QStyle.CE_ItemViewItem, options, painter)

        ctx = QAbstractTextDocumentLayout.PaintContext()

        text_rect = style.subElementRect(QStyle.SE_ItemViewItemText, options)
        painter.save()
        painter.translate(text_rect.topLeft())
        painter.setClipRect(text_rect.translated(-text_rect.topLeft()))
        doc.documentLayout().draw(painter, ctx)

        painter.restore()


class GeneItemDelegate(QItemDelegate):
    """
    https://stackoverflow.com/questions/6905147/qt-qlistwidgetitem-multiple-lines
    """

    def sizeHint(self, option, index):
        return QSize(option.rect.width(), 60)

    def paint(self, painter, option, index):

        image_space = 5
        icon = index.data(Qt.DecorationRole)
        gene_obj = index.data(Qt.DisplayRole)

        bold_font = QFont()
        bold_font.setBold(True)
        fm = QFontMetrics(QFont())

        input_name_str = 'Input name:  '
        type_of_match_str = 'Input type:  '
        gene_id_str = 'Gene ID:  '

        if not icon.isNull():
            # paint icon
            icon.paint(painter, option.rect.adjusted(image_space, image_space, -image_space, -image_space),
                       Qt.AlignVCenter | Qt.AlignLeft)
            # if image is set, change space variable
            image_space = 55

        # paint gene object data

        # input string
        r = option.rect.adjusted(image_space, 7, 0, 0)  # left, top, width and height

        painter.setFont(bold_font)
        painter.drawText(r.left(), r.top(), r.width(), r.height(),
                         Qt.AlignLeft,
                         input_name_str)

        painter.setFont(QFont())
        painter.drawText(r.left() + fm.width(input_name_str), r.top(), r.width(), r.height(),
                         Qt.AlignLeft,
                         str(gene_obj.input_name))

        # gene id string
        r = option.rect.adjusted(image_space, 22, 0, 0)  # left, top, width and height

        painter.setFont(bold_font)
        painter.drawText(r.left(), r.top(), r.width(), r.height(),
                         Qt.AlignLeft,
                         type_of_match_str)

        painter.setFont(QFont())
        painter.drawText(r.left() + fm.width(input_name_str), r.top(), r.width(), r.height(),
                         Qt.AlignLeft,
                         str(gene_obj.type_of_match) if gene_obj.type_of_match else 'Unknown')

        # type of match string
        r = option.rect.adjusted(image_space, 37, 0, 0)  # left, top, width and height

        painter.setFont(bold_font)
        painter.drawText(r.left(), r.top(), r.width(), r.height(),
                         Qt.AlignLeft,
                         gene_id_str)

        painter.setFont(QFont())
        painter.drawText(r.left() + fm.width(input_name_str), r.top(), r.width(), r.height(),
                         Qt.AlignLeft,
                         str(gene_obj.ncbi_id) if gene_obj.ncbi_id else 'Unknown')


class GeneMatcherModel(QAbstractListModel):

    def __init__(self, show_icon=True):
        QAbstractListModel.__init__(self)
        self.icon_path = ''
        self.show_icon = show_icon
        self.__items = []

    def add_rows(self, rows):
        self.__items = rows

    def get_rows(self):
        return self.__items

    def __handle_icon(self, gene_obj):
        if gene_obj.ncbi_id is None:
            if gene_obj.possible_hits:
                self.icon_path = 'icons/gene_icon_orange.svg'
            else:
                self.icon_path = 'icons/gene_icon_red.svg'
        else:
            self.icon_path = 'icons/gene_icon_green.svg'

        return QIcon(os.path.join(os.path.dirname(__file__), self.icon_path))

    def rowCount(self, *args, **kwargs):
        return len(self.__items)

    def data(self, model_index, role=None):
        # check if data is set
        if not self.__items:
            return QVariant()

        # return empty QVariant if model index is unknown
        if not model_index.isValid() or not (0 <= model_index.row() < len(self.__items)):
            return QVariant()

        gene_obj = self.__items[model_index.row()]

        if role == Qt.DisplayRole:
            return gene_obj
        elif role == Qt.DecorationRole and self.show_icon:
            return self.__handle_icon(gene_obj)


class ExtendedTableView(QWidget):

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.ow = kwargs.get('parent', None)

        # set layout
        layout = QVBoxLayout()
        layout.setContentsMargins(0, 0, 0, 0)
        self.setLayout(layout)

        # set splitter
        self.splitter = QSplitter()
        self.splitter.setOrientation(Qt.Horizontal)

        # data models
        self.genes_model = None
        self.info_model = None

        # left side list view
        self.genes_view = QTableView()
        self.genes_view.horizontalHeader().hide()

        self.genes_view.setItemDelegate(GeneItemDelegate())
        self.genes_view.horizontalHeader().setSectionResizeMode(QHeaderView.Stretch)

        # right side list view
        self.info_view = QTableView()
        self.info_view.setItemDelegate(HTMLDelegate())
        self.info_view.horizontalHeader().hide()

        self.info_view.horizontalHeader().setSectionResizeMode(QHeaderView.Stretch)

        self.splitter.addWidget(self.genes_view)
        self.splitter.addWidget(self.info_view)

        # self.splitter.setStretchFactor(0, 60)
        # self.splitter.setStretchFactor(1, 40)

        self.layout().addWidget(self.splitter)

    def set_genes_model(self, rows):
        self.genes_model = GeneMatcherModel()
        self.genes_model.add_rows(rows)

    def get_selected_gens(self):
        # return a list of QModelIndex
        return self.genes_selection_model().selectedRows()

    def reset_genes_model(self):
        if self.genes_model:
            self.genes_model.deleteLater()
            self.genes_model = None

    def genes_selection_model(self):
        return self.genes_view.selectionModel()

    def reset_info_model(self):
        if self.info_model:
            self.info_model.deleteLater()
            self.info_model = None
            self.info_view.setModel(None)

    def set_info_model(self, rows):
        unique, partial, unknown = range(len(self.ow.filter_labels))

        if self.ow.selected_filter == unique:
            # create model
            self.info_model = GeneMatcherModel(show_icon=False)
            # add rows
            self.info_model.add_rows(rows)
            # add model to the view
            self.info_view.setModel(self.info_model)
            # disable selection of gene info cards
            self.info_view.setSelectionMode(QAbstractItemView.NoSelection)
            # call sizeHint function
            self.info_view.resizeRowsToContents()
        else:
            self.reset_info_model()


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

        self.proxy_model = FilterProxyModel(self)
        self.extended_view = ExtendedTableView(parent=self)
        self.extended_view.genes_view.setModel(self.proxy_model)
        self.extended_view.genes_selection_model().selectionChanged.connect(self.__selection_changed)

        self.mainArea.layout().addWidget(self.extended_view, 1)

    def __reset_widget_state(self):
        self.Outputs.custom_data_table.send(None)
        self.proxy_model.setSourceModel(None)
        self.extended_view.reset_genes_model()
        self.extended_view.reset_info_model()

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
        self.extended_view.set_genes_model(self.gene_matcher.genes)
        self.proxy_model.setSourceModel(self.extended_view.genes_model)
        self.extended_view.genes_view.resizeRowsToContents()
        self.commit()

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
        self.__reset_widget_state()
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
