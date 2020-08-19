""" Gene Expression Omnibus datasets widget """
import sys
from types import SimpleNamespace
from typing import Any, Optional, DefaultDict
from functools import lru_cache
from collections import OrderedDict, defaultdict

import numpy as np
import requests

from AnyQt.QtGui import QFont, QColor
from AnyQt.QtCore import Qt, QSize, QVariant, QModelIndex
from AnyQt.QtWidgets import (
    QStyle,
    QSplitter,
    QTableView,
    QTreeWidget,
    QTreeWidgetItem,
    QAbstractItemView,
    QAbstractScrollArea,
)

from Orange.data import Table
from Orange.widgets.gui import (
    LinkRole,
    IndicatorItemDelegate,
    LinkStyledItemDelegate,
    rubber,
    lineEdit,
    separator,
    widgetBox,
    auto_commit,
    widgetLabel,
    radioButtonsInBox,
)
from Orange.widgets.utils import itemmodels
from Orange.widgets.widget import Msg, OWWidget
from Orange.widgets.settings import Setting
from Orange.widgets.utils.signals import Output
from Orange.widgets.utils.concurrent import TaskState, ConcurrentWidgetMixin

from orangecontrib.bioinformatics.geo import is_cached, pubmed_url, local_files
from orangecontrib.bioinformatics.geo.dataset import GDSInfo, get_samples, dataset_download


class Result(SimpleNamespace):
    gds_dataset: Optional[Table] = None


def run_download_task(gds_id: str, samples: DefaultDict[str, list], transpose: bool, state: TaskState):
    res = Result()
    current_iter = 0
    max_iter = 102

    def callback():
        nonlocal current_iter
        current_iter += 1
        state.set_progress_value(100 * (current_iter / max_iter))

    state.set_status("Downloading...")
    res.gds_dataset = dataset_download(gds_id, samples, transpose=transpose, callback=callback)
    return res


class GEODatasetsModel(itemmodels.PyTableModel):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

        self.setHorizontalHeaderLabels(
            ['', 'ID', 'PubMedID', 'Organism', 'Samples', 'Features', 'Genes', 'Subsets', 'Title']
        )

        (
            self.indicator_col,
            self.gds_id_col,
            self.pubmed_id_col,
            self.organism_col,
            self.samples_col,
            self.features_col,
            self.genes_col,
            self.subsets_col,
            self.title_col,
        ) = range(9)

        self.info = None
        self.table = None
        self._sort_column = self.gds_id_col
        self._sort_order = Qt.AscendingOrder

        self.font = QFont()
        self.font.setUnderline(True)
        self.color = QColor(Qt.blue)

        @lru_cache(maxsize=10000)
        def _row_instance(row, column):
            return self[int(row)][int(column)]

        self._row_instance = _row_instance

    def initialize(self, info: OrderedDict):
        self.info = info

        def _gds_to_row(gds: dict):
            gds_id = gds['name']
            title = gds['title']
            organism = gds['sample_organism']
            samples = len(get_samples(gds))
            features = gds['variables']
            genes = gds['genes']
            subsets = len(gds['subsets'])

            pubmed = gds.get('pubmed_id', '')
            pubmed_id = pubmed
            if isinstance(pubmed, list) and len(pubmed) > 0:
                pubmed_id = pubmed[0]

            return [
                ' ' if is_cached(gds_id) else '',
                gds_id,
                pubmed_id,
                organism,
                samples,
                features,
                genes,
                subsets,
                title,
            ]

        self.table = np.asarray([_gds_to_row(gds) for gds in info.values()])
        self.show_table()

    def _argsortData(self, data: np.ndarray, order) -> Optional[np.ndarray]:
        if not data.size:
            return

        # store order choice.
        self._sort_column = column = self.sortColumn()
        self._sort_order = self.sortOrder()

        if column == self.gds_id_col:
            data = np.char.replace(data, 'GDS', '')
            data = data.astype(int)

        elif column in (self.samples_col, self.features_col, self.genes_col, self.subsets_col, self.pubmed_id_col):
            data[data == ''] = '0'
            data = data.astype(int)

        indices = np.argsort(data, kind='mergesort')
        if order == Qt.DescendingOrder:
            indices = indices[::-1]
        return indices

    def columnCount(self, parent=QModelIndex()):
        return 0 if parent.isValid() else self._table.shape[1]

    def data(
        self,
        index,
        role,
        _str=str,
        _Qt_DisplayRole=Qt.DisplayRole,  # noqa: N803
        _Qt_EditRole=Qt.EditRole,
        _Qt_FontRole=Qt.FontRole,
        _Qt_ForegroundRole=Qt.ForegroundRole,
        _LinkRolee=LinkRole,
        _recognizedRoles=frozenset([Qt.DisplayRole, Qt.EditRole, Qt.FontRole, Qt.ForegroundRole, LinkRole]),
    ):

        if role not in _recognizedRoles:
            return None

        row, col = index.row(), index.column()
        if not 0 <= row <= self.rowCount():
            return None
        row = self.mapToSourceRows(row)

        try:
            # value = self[row][col]
            value = self._row_instance(row, col)
        except IndexError:
            return

        if role == Qt.DisplayRole:
            return QVariant(str(value))
        elif role == Qt.ToolTipRole:
            return QVariant(str(value))

        if col == self.pubmed_id_col:
            if role == _Qt_ForegroundRole:
                return self.color
            elif role == _Qt_FontRole:
                return self.font
            elif role == _LinkRolee:
                return pubmed_url.format(value)

    def get_row_index(self, gds_name):
        # test = self._table[self._table[:, 1] == gds_name, :]
        rows, _ = np.where(np.isin(self._table, gds_name))
        if rows is not None and len(rows) > 0:
            return self.mapFromSourceRows(rows[0])

    def filter_table(self, filter_pattern: str):
        selection = np.full(self.table.shape, True)
        for search_word in filter_pattern.split():
            match_result = np.core.defchararray.find(np.char.lower(self.table), search_word.lower()) >= 0
            selection = selection & match_result
        return selection

    def update_cache_indicator(self):
        self.table[:, 0] = [' ' if is_cached(gid) else '' for gid in self.table[:, self.gds_id_col]]

    def show_table(self, filter_pattern=''):
        # clear cache if model changes
        self._row_instance.cache_clear()
        self.wrap(self.table[self.filter_table(filter_pattern).any(axis=1), :])
        self.sort(self._sort_column, self._sort_order)


class OWGEODatasets(OWWidget, ConcurrentWidgetMixin):
    name = "GEO Data Sets"
    description = "Access to Gene Expression Omnibus data sets."
    icon = "icons/OWGEODatasets.svg"
    priority = 10

    class Warning(OWWidget.Warning):
        using_local_files = Msg("Can't connect to serverfiles. Using cached files.")

    class Error(OWWidget.Error):
        no_connection = Msg("Widget can't connect to serverfiles.")

    class Outputs:
        gds_data = Output("Expression Data", Table)

    search_pattern = Setting('')
    auto_commit = Setting(True)
    genes_as_rows = Setting(False)
    selected_gds = Setting(None)
    gds_selection_states = Setting({})
    splitter_settings = Setting(
        (
            b'\x00\x00\x00\xff\x00\x00\x00\x00\x00\x00\x00\x02\x00\x00\x01'
            b'\xea\x00\x00\x00\xd7\x01\x00\x00\x00\x07\x01\x00\x00\x00\x02',
            b'\x00\x00\x00\xff\x00\x00\x00\x00\x00\x00\x00\x02\x00\x00\x01'
            b'\xb5\x00\x00\x02\x10\x01\x00\x00\x00\x07\x01\x00\x00\x00\x01',
        )
    )

    def __init__(self):
        OWWidget.__init__(self)
        ConcurrentWidgetMixin.__init__(self)

        try:
            self.gds_info: Optional[GDSInfo] = GDSInfo()
        except requests.exceptions.ConnectionError:
            self.gds_info = {}
            self.Error.no_connection()
            return

        self.gds_data: Optional[Table] = None

        # Control area
        box = widgetBox(self.controlArea, 'Info', addSpace=True)
        self.infoBox = widgetLabel(box, 'Initializing\n\n')

        box = widgetBox(self.controlArea, 'Output', addSpace=True)
        radioButtonsInBox(box, self, 'genes_as_rows', ['Samples in rows', 'Genes in rows'], callback=self._run)
        separator(box)

        rubber(self.controlArea)
        auto_commit(self.controlArea, self, 'auto_commit', '&Commit', box=False)

        # Main Area

        # Filter widget
        self.filter = lineEdit(
            self.mainArea, self, 'search_pattern', 'Filter:', callbackOnType=True, callback=self._apply_filter
        )
        self.mainArea.layout().addWidget(self.filter)

        splitter_vertical = QSplitter(Qt.Vertical, self.mainArea)

        self.mainArea.layout().addWidget(splitter_vertical)

        # set table view
        self.table_view = QTableView(splitter_vertical)
        self.table_view.setShowGrid(False)
        self.table_view.setSortingEnabled(True)
        self.table_view.sortByColumn(1, Qt.AscendingOrder)
        self.table_view.setAlternatingRowColors(True)
        self.table_view.verticalHeader().setVisible(False)
        self.table_view.setSelectionBehavior(QAbstractItemView.SelectRows)
        self.table_view.setSelectionMode(QAbstractItemView.SingleSelection)
        self.table_view.viewport().setMouseTracking(True)
        self.table_view.setSizeAdjustPolicy(QAbstractScrollArea.AdjustToContents)

        self.table_model = GEODatasetsModel()
        self.table_model.initialize(self.gds_info)
        self.table_view.setModel(self.table_model)

        self.table_view.horizontalHeader().setStretchLastSection(True)
        self.table_view.resizeColumnsToContents()

        v_header = self.table_view.verticalHeader()
        option = self.table_view.viewOptions()
        size = self.table_view.style().sizeFromContents(QStyle.CT_ItemViewItem, option, QSize(20, 20), self.table_view)

        v_header.setDefaultSectionSize(size.height() + 2)
        v_header.setMinimumSectionSize(5)

        # set item delegates
        self.table_view.setItemDelegateForColumn(
            self.table_model.pubmed_id_col, LinkStyledItemDelegate(self.table_view)
        )
        self.table_view.setItemDelegateForColumn(self.table_model.gds_id_col, LinkStyledItemDelegate(self.table_view))
        self.table_view.setItemDelegateForColumn(
            self.table_model.indicator_col, IndicatorItemDelegate(self.table_view, role=Qt.DisplayRole)
        )

        splitter_horizontal = QSplitter(Qt.Horizontal, splitter_vertical)

        # Description Widget
        box = widgetBox(splitter_horizontal, 'Description')
        self.description_widget = widgetLabel(box, '')
        self.description_widget.setWordWrap(True)
        rubber(box)

        # Sample Annotations Widget
        box = widgetBox(splitter_horizontal, 'Sample Annotations')
        self.annotations_widget = QTreeWidget(box)
        self.annotations_widget.setHeaderLabels(['Type (Sample annotations)', 'Sample count'])
        self.annotations_widget.setRootIsDecorated(True)
        box.layout().addWidget(self.annotations_widget)
        self._annotations_updating = False
        self.annotations_widget.itemChanged.connect(self.on_annotation_selection_changed)
        self.splitters = splitter_vertical, splitter_horizontal

        for sp, setting in zip(self.splitters, self.splitter_settings):
            sp.splitterMoved.connect(self._splitter_moved)
            sp.restoreState(setting)

        self.table_view.selectionModel().selectionChanged.connect(self.on_gds_selection_changed)
        self._apply_filter()

        self.commit()

    def _splitter_moved(self, *args):
        self.splitter_settings = [bytes(sp.saveState()) for sp in self.splitters]

    def _set_description_widget(self):
        self.description_widget.setText(self.selected_gds.get('description', 'Description not available.'))

    def _set_annotations_widget(self, gds):
        self._annotations_updating = True
        self.annotations_widget.clear()

        annotations = defaultdict(set)
        subsets_count = {}

        for desc in gds['subsets']:
            annotations[desc['type']].add(desc['description'])
            subsets_count[desc['description']] = str(len(desc['sample_id']))

        for _type, subsets in annotations.items():
            key = (gds["name"], _type)
            parent = QTreeWidgetItem(self.annotations_widget, [_type])
            parent.key = key
            for subset in subsets:
                key = (gds['name'], _type, subset)
                item = QTreeWidgetItem(parent, [subset, subsets_count.get(subset, '')])
                item.setFlags(item.flags() | Qt.ItemIsUserCheckable)
                item.setCheckState(0, self.gds_selection_states.get(key, Qt.Checked))
                item.key = key

        self._annotations_updating = False
        self.annotations_widget.expandAll()
        for i in range(self.annotations_widget.columnCount()):
            self.annotations_widget.resizeColumnToContents(i)

    def _set_selection(self):
        if self.selected_gds is not None:
            index = self.table_model.get_row_index(self.selected_gds.get('name'))
            if index is not None:
                self.table_view.selectionModel().blockSignals(True)
                self.table_view.selectRow(index)
                self._handle_selection_changed()
                self.table_view.selectionModel().blockSignals(False)

    def _handle_selection_changed(self):
        if self.table_model.table is not None:
            selection = self.table_view.selectionModel().selectedRows(self.table_model.gds_id_col)
            selected_gds_name = selection[0].data() if len(selection) > 0 else None

            if selected_gds_name:
                self.selected_gds = self.table_model.info.get(selected_gds_name)
                self._set_annotations_widget(self.selected_gds)
                self._set_description_widget()
            else:
                self.annotations_widget.clear()
                self.description_widget.clear()

            self.update_info()

    def _apply_filter(self):
        if self.table_model.table is not None:
            self.table_model.show_table(filter_pattern=str(self.search_pattern))
            self._set_selection()
            self.update_info()

    def _run(self):
        self.Warning.using_local_files.clear()
        if self.selected_gds is not None:
            self.gds_data = None
            self.start(
                run_download_task, self.selected_gds.get('name'), self.get_selected_samples(), self.genes_as_rows
            )

    def on_gds_selection_changed(self):
        self._handle_selection_changed()
        self.commit()

    def on_annotation_selection_changed(self):
        if self._annotations_updating:
            return
        for i in range(self.annotations_widget.topLevelItemCount()):
            item = self.annotations_widget.topLevelItem(i)
            if 'key' in item.__dict__:
                self.gds_selection_states[item.key] = item.checkState(0)
            for j in range(item.childCount()):
                child = item.child(j)
                if 'key' in child.__dict__:
                    self.gds_selection_states[child.key] = child.checkState(0)

        self.commit()

    def update_info(self):
        all_gds = len(self.table_model.info)
        text = "{} datasets\n{} datasets cached\n".format(all_gds, len(local_files.listfiles()))
        filtered = self.table_view.model().rowCount()
        if all_gds != filtered:
            text += "{} after filtering".format(filtered)
        self.infoBox.setText(text)

    def get_selected_samples(self):
        """
        Return the currently selected sample annotations.

        The return value is a list of selected (sample type, sample value)
        tuples.

        .. note:: if some Sample annotation type has no selected values.
                  this method will return all values for it.

        TODO: this could probably be simplified.

        """

        def childiter(item):
            """ Iterate over the children of an QTreeWidgetItem instance.
            """
            for i in range(item.childCount()):
                yield item.child(i)

        samples = []
        unused_types = []
        used_types = []

        for stype in childiter(self.annotations_widget.invisibleRootItem()):
            selected_values = []
            all_values = []
            for sval in childiter(stype):
                value = (str(stype.text(0)), str(sval.text(0)))
                if self.gds_selection_states.get(sval.key, True):
                    selected_values.append(value)
                all_values.append(value)
            if selected_values:
                samples.extend(selected_values)
                used_types.append(str(stype.text(0)))
            else:
                # If no sample of sample type is selected we don't filter on it.
                samples.extend(all_values)
                unused_types.append(str(stype.text(0)))

        _samples = defaultdict(list)
        for sample, sample_type in samples:
            _samples[sample].append(sample_type)
        return _samples

    def commit(self):
        self._run()

    def on_exception(self, ex: Exception):
        self.Warning.using_local_files()

    def on_done(self, result: Result):
        assert isinstance(result.gds_dataset, Table)
        self.gds_data = result.gds_dataset

        if self.gds_info:
            self.table_model.update_cache_indicator()
            self._apply_filter()

        self.Outputs.gds_data.send(self.gds_data)

    def on_partial_result(self, result: Any) -> None:
        pass

    def onDeleteWidget(self):
        self.shutdown()
        super().onDeleteWidget()

    def send_report(self):
        self.report_items(
            "GEO Dataset",
            [
                ("ID", self.selected_gds['name']),
                ("Title", self.selected_gds['title']),
                ("Organism", self.selected_gds['sample_organism']),
            ],
        )
        self.report_items(
            "Data",
            [
                ("Samples", self.selected_gds['sample_count']),
                ("Features", self.selected_gds['variables']),
                ("Genes", self.selected_gds['genes']),
            ],
        )
        self.report_name("Sample annotations")
        subsets = defaultdict(list)
        for subset in self.selected_gds['subsets']:
            subsets[subset['type']].append((subset['description'], len(subset['sample_id'])))
        self.report_html += "<ul>"
        for _type in subsets:
            self.report_html += "<b>" + _type + ":</b></br>"
            for desc, count in subsets[_type]:
                self.report_html += 9 * "&nbsp" + "<b>{}:</b> {}</br>".format(desc, count)
        self.report_html += "</ul>"


if __name__ == "__main__":

    def main_test():
        from AnyQt.QtWidgets import QApplication

        app = QApplication([])
        w = OWGEODatasets()
        w.show()
        w.raise_()
        r = app.exec_()
        w.saveSettings()
        return r

    sys.exit(main_test())
