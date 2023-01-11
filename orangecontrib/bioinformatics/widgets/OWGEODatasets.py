""" Gene Expression Omnibus datasets widget """
import sys
from types import SimpleNamespace
from typing import Any, Optional, DefaultDict
from collections import defaultdict

import requests

from AnyQt.QtCore import Qt
from AnyQt.QtWidgets import (
    QSplitter,
    QTreeWidget,
    QTreeWidgetItem,
    QAbstractItemView,
    QAbstractScrollArea,
)

from orangewidget.utils.itemdelegates import DataDelegate

from Orange.data import Table
from Orange.widgets.gui import (
    LinkRole,
    IndicatorItemDelegate,
    rubber,
    lineEdit,
    separator,
    widgetBox,
    auto_commit,
    widgetLabel,
    radioButtonsInBox,
)
from Orange.widgets.widget import Msg, OWWidget
from Orange.widgets.settings import Setting
from Orange.widgets.utils.signals import Output
from Orange.widgets.utils.tableview import TableView
from Orange.widgets.utils.concurrent import TaskState, ConcurrentWidgetMixin

from orangecontrib.bioinformatics.geo import is_cached, pubmed_url, local_files
from orangecontrib.bioinformatics.geo.dataset import (
    GDSInfo,
    get_samples,
    dataset_download,
)
from orangecontrib.bioinformatics.widgets.utils.gui import FilterProxyModel
from orangecontrib.bioinformatics.widgets.utils.itemmodels import TableModel
from orangecontrib.bioinformatics.widgets.utils.itemdelegates import (
    LinkStyledItemDelegate,
)


class Result(SimpleNamespace):
    gds_dataset: Optional[Table] = None


def run_download_task(
    gds_id: str, samples: DefaultDict[str, list], transpose: bool, state: TaskState
):
    res = Result()
    current_iter = 0
    max_iter = 102

    def callback():
        nonlocal current_iter
        current_iter += 1
        state.set_progress_value(100 * (current_iter / max_iter))

    state.set_status("Downloading...")
    res.gds_dataset = dataset_download(
        gds_id, samples, transpose=transpose, callback=callback
    )
    return res


class GEODatasetsModel(TableModel):
    (
        indicator_col,
        gds_id_col,
        pubmed_id_col,
        organism_col,
        samples_col,
        features_col,
        genes_col,
        subsets_col,
        title_col,
    ) = range(9)

    def __init__(self, gds: GDSInfo):
        items = list(gds.values())
        pubmed_ids = [g.get("pubmed_id", None) for g in items]

        def title(name: str):
            return {Qt.DisplayRole: name}

        def pubmed_id(row: int) -> Optional[str]:
            pubmed = pubmed_ids[row]
            if isinstance(pubmed, list):
                return pubmed[0] if len(pubmed) > 0 else None
            else:
                return pubmed

        def render_pubmed_url(row):
            id_ = pubmed_id(row)
            if id_:
                return pubmed_url.format(id_)
            return None

        columns = [
            TableModel.Column(
                title(""),
                {
                    Qt.DisplayRole: lambda row: " "
                    if is_cached(items[row]["name"])
                    else "",
                    Qt.UserRole: lambda row: items[row],
                },
            ),
            TableModel.Column(
                title("ID"),
                {
                    Qt.DisplayRole: lambda row: items[row]["name"],
                    Qt.UserRole: lambda row: items[row],
                },
            ),
            TableModel.Column(
                title("PubMedID"),
                {
                    Qt.DisplayRole: pubmed_id,
                    LinkRole: render_pubmed_url,
                    Qt.ToolTipRole: render_pubmed_url,
                },
            ),
            TableModel.Column(
                title("Organism"),
                {
                    Qt.DisplayRole: lambda row: items[row]["sample_organism"],
                },
            ),
            TableModel.Column(
                title("Samples"),
                {
                    Qt.DisplayRole: lambda row: len(get_samples(items[row])),
                },
            ),
            TableModel.Column(
                title("Features"),
                {
                    Qt.DisplayRole: lambda row: items[row]["variables"],
                },
            ),
            TableModel.Column(
                title("Genes"),
                {
                    Qt.DisplayRole: lambda row: items[row]["genes"],
                },
            ),
            TableModel.Column(
                title("Subsets"),
                {
                    Qt.DisplayRole: lambda row: len(items[row]["subsets"]),
                },
            ),
            TableModel.Column(
                title("Title"),
                {
                    Qt.DisplayRole: lambda row: items[row]["title"],
                    Qt.ToolTipRole: lambda row: items[row]["title"],
                },
            ),
        ]
        super().__init__(len(items), columns)
        self.info = gds

    def sortColumnData(self, column):
        if column == self.gds_id_col:
            items = [g["name"] for g in self.info.values()]
            items = [int(name.strip('GDS')) for name in items]
            return items
        if column == self.pubmed_id_col:
            items = [
                self.data(self.index(i, column), Qt.DisplayRole)
                for i in self.mapFromSourceRows(range(self.rowCount()))
            ]
            items = [int(pid or 0) for pid in items]
            return items
        return super().sortColumnData(column)

    def update_cache_indicator(self):
        self.dataChanged.emit(self.index(0, 0), self.index(self.rowCount() - 1, 0))


class FilterProxyModel(FilterProxyModel):
    def sort(self, column: int, order: Qt.SortOrder = Qt.AscendingOrder) -> None:
        """
        Reimplemented.

        Dispatch the sorting to the source model.
        """
        source = self.sourceModel()
        if source is not None:
            source.sort(column, order)


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
        self.__updating_filter = False

        # Control area
        box = widgetBox(self.controlArea, 'Info', addSpace=True)
        self.infoBox = widgetLabel(box, 'Initializing\n\n')

        box = widgetBox(self.controlArea, 'Output', addSpace=True)
        radioButtonsInBox(
            box,
            self,
            'genes_as_rows',
            ['Samples in rows', 'Genes in rows'],
            callback=self._run,
        )
        separator(box)

        rubber(self.controlArea)
        auto_commit(self.controlArea, self, 'auto_commit', '&Commit', box=False)

        # Main Area

        # Filter widget
        self.filter = lineEdit(
            self.mainArea,
            self,
            'search_pattern',
            'Filter:',
            callbackOnType=True,
            callback=self._apply_filter,
        )
        self.mainArea.layout().addWidget(self.filter)

        splitter_vertical = QSplitter(Qt.Vertical, self.mainArea)

        self.mainArea.layout().addWidget(splitter_vertical)

        # set table view
        self.table_view = TableView(
            splitter_vertical,
            showGrid=False,
            sortingEnabled=True,
            alternatingRowColors=True,
            selectionBehavior=QAbstractItemView.SelectRows,
            selectionMode=QAbstractItemView.SingleSelection,
            sizeAdjustPolicy=QAbstractScrollArea.AdjustToContents,
        )
        self.table_view.sortByColumn(1, Qt.AscendingOrder)
        self.table_view.verticalHeader().setVisible(False)
        self.table_view.viewport().setMouseTracking(True)

        self.table_model = GEODatasetsModel(self.gds_info)
        self.proxy_model = FilterProxyModel()
        self.proxy_model.setSourceModel(self.table_model)
        self.table_view.setModel(self.proxy_model)

        self.table_view.horizontalHeader().setStretchLastSection(True)
        self.table_view.resizeColumnsToContents()

        # set item delegates
        self.table_view.setItemDelegate(DataDelegate(self.table_view))
        self.table_view.setItemDelegateForColumn(
            self.table_model.pubmed_id_col, LinkStyledItemDelegate(self.table_view)
        )
        self.table_view.setItemDelegateForColumn(
            self.table_model.gds_id_col, LinkStyledItemDelegate(self.table_view)
        )
        self.table_view.setItemDelegateForColumn(
            self.table_model.indicator_col,
            IndicatorItemDelegate(self.table_view, role=Qt.DisplayRole),
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
        self.annotations_widget.setHeaderLabels(
            ['Type (Sample annotations)', 'Sample count']
        )
        self.annotations_widget.setRootIsDecorated(True)
        box.layout().addWidget(self.annotations_widget)
        self._annotations_updating = False
        self.annotations_widget.itemChanged.connect(
            self.on_annotation_selection_changed
        )
        self.splitters = splitter_vertical, splitter_horizontal

        for sp, setting in zip(self.splitters, self.splitter_settings):
            sp.splitterMoved.connect(self._splitter_moved)
            sp.restoreState(setting)

        self.table_view.selectionModel().selectionChanged.connect(
            self.__on_selection_changed
        )
        self._apply_filter()

        self.commit()

    def _splitter_moved(self, *args):
        self.splitter_settings = [bytes(sp.saveState()) for sp in self.splitters]

    def _set_description_widget(self):
        self.description_widget.setText(
            self.selected_gds.get('description', 'Description not available.')
        )

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
            index = None
            ids = self.proxy_model.match(
                self.proxy_model.index(0, 1),
                Qt.DisplayRole,
                self.selected_gds["name"],
                1,
                Qt.MatchCaseSensitive | Qt.MatchExactly,
            )
            if ids:
                index = ids[0]

            if index is not None:
                self.table_view.selectionModel().blockSignals(True)
                self.table_view.selectRow(index.row())
                self._handle_selection_changed()
                self.table_view.selectionModel().blockSignals(False)

    def _handle_selection_changed(self):
        if self.table_model is not None:
            selection = self.table_view.selectionModel().selectedRows(
                self.table_model.gds_id_col
            )
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
        keys = self.search_pattern.lower().split()

        def filter(gds):
            source = (gds["name"], gds["title"], gds["sample_organism"])
            source = [t.lower() for t in source]
            return all(any(k in s for s in source) for k in keys)

        changed = False

        def current_changed():
            nonlocal changed
            changed = True

        selmodel = self.table_view.selectionModel()
        selmodel.selectionChanged.connect(current_changed)
        try:
            self.__updating_filter = True
            self.proxy_model.set_filters(
                [
                    FilterProxyModel.Filter(
                        self.table_model.gds_id_col, Qt.UserRole, filter
                    )
                ]
            )
            if changed:
                self.table_view.clearSelection()
                self._handle_selection_changed()
        finally:
            selmodel.selectionChanged.disconnect(current_changed)
            self.__updating_filter = False

        self._set_selection()

    def _run(self):
        self.Warning.using_local_files.clear()
        if self.selected_gds is not None:
            self.gds_data = None
            self.start(
                run_download_task,
                self.selected_gds.get('name'),
                self.get_selected_samples(),
                self.genes_as_rows,
            )

    def __on_selection_changed(self):
        if not self.__updating_filter:
            self.on_gds_selection_changed()

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
        text = "{} datasets\n{} datasets cached\n".format(
            all_gds, len(local_files.listfiles())
        )
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
            """Iterate over the children of an QTreeWidgetItem instance."""
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
            subsets[subset['type']].append(
                (subset['description'], len(subset['sample_id']))
            )
        self.report_html += "<ul>"
        for _type in subsets:
            self.report_html += "<b>" + _type + ":</b></br>"
            for desc, count in subsets[_type]:
                self.report_html += 9 * "&nbsp" + "<b>{}:</b> {}</br>".format(
                    desc, count
                )
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
