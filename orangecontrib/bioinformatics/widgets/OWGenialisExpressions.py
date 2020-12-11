import io
from enum import IntEnum
from typing import Any, Dict, List, Tuple, Optional, NamedTuple
from datetime import datetime

import numpy as np
import pandas as pd
from resdk.resources.data import Data
from dateutil.relativedelta import relativedelta

from AnyQt.QtCore import (
    Qt,
    QSize,
    QObject,
    QAbstractAnimation,
    QPropertyAnimation,
    QParallelAnimationGroup,
    pyqtSignal,
)
from AnyQt.QtWidgets import (
    QFrame,
    QLabel,
    QDialog,
    QComboBox,
    QLineEdit,
    QTableView,
    QHBoxLayout,
    QHeaderView,
    QPushButton,
    QScrollArea,
    QSizePolicy,
    QToolButton,
    QVBoxLayout,
    QAbstractItemView,
)

from Orange.data import Table, Domain
from Orange.widgets import gui, widget, settings
from Orange.widgets.widget import Msg, Output, StateInfo, OWComponent
from Orange.data.pandas_compat import table_from_frame
from Orange.widgets.credentials import CredentialManager
from Orange.widgets.utils.concurrent import TaskState, ConcurrentWidgetMixin
from Orange.widgets.utils.itemmodels import PyTableModel

from orangecontrib.bioinformatics.resolwe import ResolweAPI, ResolweAuthException, ResolweDataObjectsNotFound, connect
from orangecontrib.bioinformatics.ncbi.gene import GeneMatcher
from orangecontrib.bioinformatics.preprocess import ZScore, LogarithmicScale, QuantileTransform, QuantileNormalization
from orangecontrib.bioinformatics.ncbi.taxonomy import species_name_to_taxid
from orangecontrib.bioinformatics.resolwe.resapi import DEFAULT_URL, RESOLWE_URLS, CREDENTIAL_MANAGER_SERVICE
from orangecontrib.bioinformatics.widgets.utils.data import TableAnnotation
from orangecontrib.bioinformatics.resolwe.sample_metadata import merge_clinical_data

PAST_HOUR = relativedelta(hours=-1)
PAST_DAY = relativedelta(days=-1)
PAST_WEEK = relativedelta(weeks=-1)
PAST_MONTH = relativedelta(months=-1)


class FilterByDateModified(IntEnum):
    any_time = 0
    past_hour = 1
    past_day = 2
    past_week = 3
    past_month = 4

    @staticmethod
    def values():
        now = datetime.now()
        return [None, now + PAST_HOUR, now + PAST_DAY, now + PAST_WEEK, now + PAST_MONTH]

    @staticmethod
    def labels():
        return ['Any time', 'Past hour', 'Past 24 hours', 'Past week', 'Past month']


class SortBy(IntEnum):
    relevance = 0
    name_a_to_z = 1
    name_z_to_a = 2
    contributor_a_to_z = 3
    contributor_z_to_a = 4
    newest_first = 5
    oldest_first = 6

    @staticmethod
    def values():
        return ['', 'name', '-name', 'contributor__last_name', '-contributor__last_name', '-modified', 'modified']

    @staticmethod
    def labels():
        return [
            'Relevance',
            'Name A-Z',
            'Name Z-A',
            'Contributor A-Z',
            'Contributor Z-A',
            'Date modified (newest first)',
            'Date modified (oldest first)',
        ]


class ItemsPerPage(IntEnum):
    min = 0
    med = 1
    max = 2

    @staticmethod
    def values():
        return [20, 50, 100]


class QuantileTransformDist(IntEnum):
    normal = 0
    uniform = 1

    @staticmethod
    def values():
        return ['normal', 'uniform']


class TableHeader(IntEnum):
    id = 0
    slug = 1
    name = 2
    samples = 3
    species = 4
    created = 5
    modified = 6
    contributor = 7
    description = 8
    tags = 9

    @staticmethod
    def labels():
        return [
            'Id',
            'Slug',
            'Name',
            'Samples',
            'Species',
            'Created',
            'Modified',
            'Contributor',
            'Description',
            'Tags',
        ]


class CollapsibleFilterComponent(OWComponent, QObject):
    options_changed = pyqtSignal()

    filter_by_full_text: str = settings.Setting('', schema_only=True)
    filter_by_name: str = settings.Setting('', schema_only=True)
    filter_by_contrib: str = settings.Setting('', schema_only=True)
    filter_by_owner: str = settings.Setting('', schema_only=True)
    filter_by_modified: int = settings.Setting(FilterByDateModified.any_time, schema_only=True)
    sort_by: int = settings.Setting(SortBy.newest_first, schema_only=True)

    FILTER_FULL_TEXT_LABEL = 'Search'
    TOGGLE_BTN_LABEL = 'Narrow your search'
    FILTER_NAME_LABEL = 'Filter by name'
    FILTER_CONTRIB_LABEL = 'Filter by contributor'
    FILTER_OWNER_LABEL = 'Filter by owner'
    FILTER_MODIFIED_LABEL = 'Filter by date modified'
    SORTING_LABEL = 'Sorting'

    def __init__(self, parent_widget, parent_component):
        QObject.__init__(self)
        OWComponent.__init__(self, widget=parent_widget)

        box = gui.widgetBox(parent_component, margin=0)

        self.filter_full_text = gui.lineEdit(
            box,
            self,
            'filter_by_full_text',
            label=self.FILTER_FULL_TEXT_LABEL,
            callback=self.on_filter_full_text_changed,
        )

        self.toggle_animation = QParallelAnimationGroup()
        self.toggle_button = QToolButton()
        self.toggle_button.setCheckable(True)
        self.toggle_button.setChecked(False)
        self.toggle_button.setToolButtonStyle(Qt.ToolButtonTextBesideIcon)
        self.toggle_button.setArrowType(Qt.RightArrow)
        self.toggle_button.setText(self.TOGGLE_BTN_LABEL)
        self.toggle_button.setStyleSheet('QToolButton {border: none; padding-top: 5px; }')
        self.toggle_button.setIconSize(QSize(15, 15))
        self.toggle_button.pressed.connect(self.on_toggle)

        self.collapsible_components = QScrollArea()
        self.collapsible_components.setMaximumHeight(0)
        self.collapsible_components.setMinimumHeight(0)
        self.collapsible_components.setSizePolicy(QSizePolicy.Expanding, QSizePolicy.Fixed)
        self.collapsible_components.setFrameShape(QFrame.NoFrame)

        box = gui.widgetBox(parent_component, margin=0)
        box.layout().addWidget(self.toggle_button)
        box.layout().addWidget(self.collapsible_components)

        self.toggle_animation.addAnimation(QPropertyAnimation(box, b"minimumHeight"))
        self.toggle_animation.addAnimation(QPropertyAnimation(box, b"maximumHeight"))
        self.toggle_animation.addAnimation(QPropertyAnimation(self.collapsible_components, b"maximumHeight"))

        layout = QHBoxLayout()
        left_box = gui.widgetBox(None, self, margin=0, flat=True)
        mid_box = gui.widgetBox(None, self, margin=0, flat=True)
        right_box = gui.widgetBox(None, self, margin=0, flat=True)

        self.filter_name = gui.lineEdit(
            left_box, self, 'filter_by_name', label=self.FILTER_NAME_LABEL, callback=self.on_filter_changed, addSpace=5
        )
        self.filter_contrib = gui.lineEdit(
            mid_box,
            self,
            'filter_by_contrib',
            label=self.FILTER_CONTRIB_LABEL,
            callback=self.on_filter_changed,
            addSpace=5,
        )
        self.filter_owner = gui.lineEdit(
            right_box,
            self,
            'filter_by_owner',
            label=self.FILTER_OWNER_LABEL,
            callback=self.on_filter_changed,
            addSpace=5,
        )
        self.filter_modified = gui.comboBox(
            left_box,
            self,
            'filter_by_modified',
            label=self.FILTER_MODIFIED_LABEL,
            callback=self.on_filter_changed,
            items=FilterByDateModified.labels(),
        )
        self.sorting = gui.comboBox(
            mid_box, self, 'sort_by', label=self.SORTING_LABEL, callback=self.on_filter_changed, items=SortBy.labels()
        )

        gui.rubber(left_box)
        gui.rubber(mid_box)
        gui.rubber(right_box)
        layout.addWidget(left_box)
        layout.addWidget(mid_box)
        layout.addWidget(right_box)
        self.collapsible_components.setLayout(layout)

        collapsed_height = box.layout().sizeHint().height() - self.collapsible_components.maximumHeight()
        content_height = layout.sizeHint().height()

        for i in range(self.toggle_animation.animationCount()):
            animation = self.toggle_animation.animationAt(i)
            animation.setDuration(100)
            animation.setStartValue(collapsed_height)
            animation.setEndValue(collapsed_height + content_height)

        content_animation = self.toggle_animation.animationAt(self.toggle_animation.animationCount() - 1)
        content_animation.setDuration(100)
        content_animation.setStartValue(0)
        content_animation.setEndValue(content_height)

    def on_toggle(self):
        """ Start animation """
        checked = self.toggle_button.isChecked()
        self.toggle_button.setArrowType(Qt.DownArrow if not checked else Qt.RightArrow)
        self.toggle_animation.setDirection(QAbstractAnimation.Forward if not checked else QAbstractAnimation.Backward)
        self.toggle_animation.start()

    def on_filter_full_text_changed(self):
        self.sort_by = SortBy.relevance if self.filter_full_text else SortBy.newest_first
        self.on_filter_changed()

    def on_filter_changed(self):
        self.options_changed.emit()


class PaginationComponent(OWComponent, QObject):
    options_changed = pyqtSignal()

    items_per_page: int = settings.Setting(ItemsPerPage.min, schema_only=True)
    current_page: int = settings.Setting(1, schema_only=True)
    offset: int = settings.Setting(0, schema_only=True)

    PAGE_LIMIT_LABEL = 'Items per page  '

    def __init__(self, parent_widget, parent_component):
        QObject.__init__(self)
        OWComponent.__init__(self, widget=parent_widget)

        box = gui.widgetBox(parent_component, orientation=Qt.Horizontal)

        self.page_limit = gui.radioButtons(
            box,
            self,
            'items_per_page',
            [str(val) for val in ItemsPerPage.values()],
            orientation=Qt.Horizontal,
            callback=self.on_limit_changed,
            label=self.PAGE_LIMIT_LABEL,
        )

        self.page_left_btn = QToolButton()
        self.page_left_btn.setStyleSheet('QToolButton {border: none;}')
        self.page_left_btn.setArrowType(Qt.LeftArrow)
        self.page_left_btn.pressed.connect(self.left_btn_pressed)

        self.offset_label = gui.label(None, self, str(self.current_page), labelWidth=15)
        self.offset_label.setAlignment(Qt.AlignCenter)

        self.page_right_btn = QToolButton()
        self.page_right_btn.setStyleSheet('QToolButton {border: none;}')
        self.page_right_btn.setArrowType(Qt.RightArrow)
        self.page_right_btn.pressed.connect(self.right_btn_pressed)

        parent_widget.pagination_availability.connect(self._handle_paginate_buttons)

        box.layout().addStretch(1)
        box.layout().addWidget(self.page_left_btn)
        box.layout().addWidget(self.offset_label)
        box.layout().addWidget(self.page_right_btn)

    def _handle_paginate_buttons(self, next_page: bool, previous_page: bool):
        self.page_left_btn.setEnabled(previous_page)
        self.page_right_btn.setEnabled(next_page)

    def reset_pagination(self):
        self.offset = 0
        self.current_page = 1
        self.offset_label.setText(str(self.current_page))

    def left_btn_pressed(self):
        new_offset = self.offset - ItemsPerPage.values()[self.items_per_page]

        if new_offset >= 0:
            self.current_page -= 1
            self.offset = new_offset
            self.offset_label.setText(str(self.current_page))

        self.options_changed.emit()

    def right_btn_pressed(self):
        self.current_page += 1
        self.offset += ItemsPerPage.values()[self.items_per_page]
        self.offset_label.setText(str(self.current_page))
        self.options_changed.emit()

    def on_limit_changed(self):
        self.reset_pagination()
        self.options_changed.emit()


class NormalizationComponent(OWComponent, QObject):
    """ Gene expression normalization component """

    options_changed = pyqtSignal()

    quantile_norm: bool
    quantile_norm = settings.Setting(False, schema_only=True)

    log_norm: bool
    log_norm = settings.Setting(True, schema_only=True)

    z_score_norm: bool
    z_score_norm = settings.Setting(False, schema_only=True)

    z_score_axis: int
    z_score_axis = settings.Setting(0, schema_only=True)

    quantile_transform: bool
    quantile_transform = settings.Setting(False, schema_only=True)

    quantile_transform_axis: int
    quantile_transform_axis = settings.Setting(0, schema_only=True)

    quantile_transform_dist: int
    quantile_transform_dist = settings.Setting(0, schema_only=True)

    BOX_TITLE = 'Normalization'
    QUANTILE_NORM_LABEL = 'Quantile normalization'
    LOG_NORM_LABEL = 'Log2(x+1)'
    Z_SCORE_LABEL = 'Z-score'
    QUANTILE_TRANSFORM_LABEL = 'Quantile transform'

    def __init__(self, parent_widget, parent_component):
        QObject.__init__(self)
        OWComponent.__init__(self, widget=parent_widget)

        box = gui.widgetBox(parent_component, self.BOX_TITLE, margin=3)
        gui.checkBox(box, self, 'quantile_norm', self.QUANTILE_NORM_LABEL, callback=self.options_changed.emit)
        gui.checkBox(box, self, 'log_norm', self.LOG_NORM_LABEL, callback=self.options_changed.emit)
        gui.checkBox(box, self, 'z_score_norm', self.Z_SCORE_LABEL, callback=self.on_z_score_selected)
        self.z_score_axis_btn = gui.radioButtons(
            gui.indentedBox(box),
            self,
            'z_score_axis',
            btnLabels=['columns', 'rows'],
            callback=self.options_changed.emit,
            orientation=Qt.Horizontal,
        )
        self.z_score_axis_btn.setHidden(not bool(self.z_score_norm))

        gui.checkBox(
            box,
            self,
            'quantile_transform',
            self.QUANTILE_TRANSFORM_LABEL,
            callback=self.on_quantile_transform_selected,
        )
        self.quantile_transform_axis_btn = gui.radioButtons(
            gui.indentedBox(box),
            self,
            'quantile_transform_axis',
            btnLabels=['columns', 'rows'],
            callback=self.options_changed.emit,
            orientation=Qt.Horizontal,
        )
        self.quantile_transform_dist_btn = gui.radioButtons(
            gui.indentedBox(box),
            self,
            'quantile_transform_dist',
            btnLabels=QuantileTransformDist.values(),
            callback=self.options_changed.emit,
            orientation=Qt.Horizontal,
        )

        self.quantile_transform_axis_btn.setHidden(not bool(self.quantile_transform))
        self.quantile_transform_dist_btn.setHidden(not bool(self.quantile_transform))

    def on_z_score_selected(self):
        self.z_score_axis_btn.setHidden(not bool(self.z_score_norm))
        self.options_changed.emit()

    def on_quantile_transform_selected(self):
        self.quantile_transform_axis_btn.setHidden(not bool(self.quantile_transform))
        self.quantile_transform_dist_btn.setHidden(not bool(self.quantile_transform))
        self.options_changed.emit()


class SignInForm(QDialog):
    def __init__(self, flags, *args, **kwargs):
        super().__init__(flags, *args, **kwargs)
        self.cm: CredentialManager = CredentialManager(CREDENTIAL_MANAGER_SERVICE)

        self.setWindowTitle('Sign in')
        self.setFixedSize(400, 250)

        self.server_cb_label = QLabel('Server *')
        self.server_cb = QComboBox(self)
        self.server_cb.addItems(RESOLWE_URLS)
        self.server_cb.setEditable(True)

        self.username_label = QLabel('Username *')
        self.username_line_edit = QLineEdit(self)
        self.username_line_edit.setPlaceholderText('Enter correct username')
        self.username_line_edit.returnPressed.connect(self.sign_in)
        self.username_line_edit.textChanged.connect(self.handle_sign_in_btn)

        self.password_label = QLabel('Password *')
        self.password_line_edit = QLineEdit(self)
        self.password_line_edit.setPlaceholderText('Enter correct password')
        self.password_line_edit.returnPressed.connect(self.sign_in)
        self.password_line_edit.textChanged.connect(self.handle_sign_in_btn)
        self.password_line_edit.setEchoMode(QLineEdit.Password)

        self.sign_in_btn = QPushButton('Sign in', self)
        self.sign_in_btn.setDisabled(True)
        self.sign_in_btn.clicked.connect(self.sign_in)

        self.error_msg = QLabel('Unable to log in with provided credentials.')
        self.error_msg.setStyleSheet('color:red')
        self.error_msg.hide()

        layout = QVBoxLayout(self)
        layout.addWidget(self.server_cb_label)
        layout.addWidget(self.server_cb)
        layout.addWidget(self.username_label)
        layout.addWidget(self.username_line_edit)
        layout.addWidget(self.password_label)
        layout.addWidget(self.password_line_edit)
        layout.addWidget(self.error_msg)
        layout.addStretch()
        layout.addWidget(self.sign_in_btn)

        self.resolwe_instance = None

    def handle_sign_in_btn(self):
        self.sign_in_btn.setEnabled(
            True if self.username_line_edit.text() and self.password_line_edit.text() else False
        )

    def sign_in(self):
        self.server_cb_label.setStyleSheet(None)
        self.username_label.setStyleSheet(None)
        self.password_label.setStyleSheet(None)
        self.error_msg.hide()

        server = self.server_cb.currentText()
        username = self.cm.username if self.cm.username else self.username_line_edit.text()
        password = self.cm.password if self.cm.password else self.password_line_edit.text()

        if not server:
            self.server_cb_label.setStyleSheet('color:red')
            return

        if not username:
            self.username_label.setStyleSheet('color:red')
            return

        if not password:
            self.password_label.setStyleSheet('color:red')
            return

        try:
            self.resolwe_instance = connect(username, password, url=server)
        except ResolweAuthException:
            self.error_msg.show()
            return

        self.cm.username = username
        self.cm.password = password
        self.accept()


class GenialisExpressionsModel(PyTableModel):
    def __init__(self, parent_widget, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.parent = parent_widget

    def flags(self, index):
        """
        Disable the row selection by clicking on the first column.
        """
        return Qt.ItemIsEnabled if index.column() == 0 else Qt.ItemIsEnabled | Qt.ItemIsSelectable

    def data(self, index, role=Qt.DisplayRole):
        if not index.isValid():
            return

        row, column = self.mapToSourceRows(index.row()), index.column()

        role_value = self._roleData.get(row, {}).get(column, {}).get(role)
        if role_value is not None:
            return role_value

        try:
            value = self[row][column]
            tag = self[row][TableHeader.tags]
            slug = self[row][TableHeader.slug]
            url = f'{self.parent.res.url}/{tag}/search/collection/{slug}'
        except IndexError:
            return

        if role == Qt.DisplayRole:
            if isinstance(value, datetime):
                return str(value.date().strftime('%m/%d/%Y'))
            elif column == TableHeader.description:
                if value:
                    return f'{str(value[:100])} ...'
            return str(value)

        elif role == Qt.ToolTipRole:
            if isinstance(value, datetime):
                return str(value.strftime('%b %d, %Y %H:%M'))
            if tag and column == TableHeader.id:
                return url
            return str(value)

        # elif role == Qt.TextAlignmentRole and isinstance(value, Number):
        #     return Qt.AlignRight | Qt.AlignVCenter

        elif role == gui.LinkRole:
            if tag:
                return url

    def set_data(self, collections: List[Dict[str, str]], col_to_species: Dict[str, str]):
        def model_row(collection: Dict[str, str]) -> List[Any]:
            first_name = collection.get('contributor', {}).get('first_name', '')
            last_name = collection.get('contributor', {}).get('last_name', '')
            user_name = collection.get('contributor', {}).get('username', '')
            contributor = f'{first_name} {last_name}'.strip()

            tags = collection.get('tags', [])
            tags = [tag for tag in tags if 'community' in tag]
            if len(tags) == 1:
                _, tag = tags[0].split(':')
            else:
                tag = ''

            return [
                collection['id'],
                collection.get('slug', ''),
                collection.get('name', ''),
                collection.get('entity_count', ''),
                col_to_species.get(collection['id'], ''),
                datetime.strptime(collection.get('created', ''), '%Y-%m-%dT%H:%M:%S.%f%z'),
                datetime.strptime(collection.get('modified', ''), '%Y-%m-%dT%H:%M:%S.%f%z'),
                contributor if contributor else user_name,
                collection.get('description', ''),
                tag,
            ]

        self.wrap([model_row(result) for result in collections])


class Expression(NamedTuple):
    type: str
    name: str


class Process(NamedTuple):
    type: str
    name: str


class InputAnnotation(NamedTuple):
    source: str
    species: str
    build: str


class DataOutputOptions(NamedTuple):
    expression: Tuple[Expression]
    process: Tuple[Process]
    input_annotation: Tuple[InputAnnotation]


def available_data_output_options(data_objects: List[Data]) -> DataOutputOptions:
    """
    Traverse the data objects in the selected collection and store the
    information regarding available expression types, process types and
    input annotations.

    The method returns a named tuple (`DataOutputOptions`) which used for
    creating radio buttons in the control area.
    """

    expression_types = sorted({data.output['exp_type'] for data in data_objects})
    expression_types = (Expression('rc', 'Read Counts'),) + tuple(
        Expression(exp_type, exp_type) for exp_type in expression_types
    )

    process_types = sorted({(data.process.type, data.process.name) for data in data_objects})
    process_types = tuple(Process(proc_type, proc_name) for proc_type, proc_name in process_types)

    input_annotations = sorted(
        {(data.output['source'], data.output['species'], data.output['build']) for data in data_objects}
    )
    input_annotations = tuple(InputAnnotation(source, species, build) for source, species, build in input_annotations)

    return DataOutputOptions(expression=expression_types, process=process_types, input_annotation=input_annotations)


class OWGenialisExpressions(widget.OWWidget, ConcurrentWidgetMixin):
    name = 'Genialis Expressions'
    priority = 30
    want_main_area = True
    want_control_area = True
    icon = '../widgets/icons/OWGenialisExpressions.svg'

    pagination_availability = pyqtSignal(bool, bool)

    norm_component = settings.SettingProvider(NormalizationComponent)
    pagination_component = settings.SettingProvider(PaginationComponent)
    filter_component = settings.SettingProvider(CollapsibleFilterComponent)

    exp_type: int
    exp_type = settings.Setting(0, schema_only=True)

    proc_type: int
    proc_type = settings.Setting(0, schema_only=True)

    input_annotation: int
    input_annotation = settings.Setting(0, schema_only=True)

    selected_orange_data_object: int
    selected_orange_data_object = settings.Setting(0, schema_only=True)

    auto_commit: bool
    auto_commit = settings.Setting(False, schema_only=True)

    class Outputs:
        table = Output('Expressions', Table)

    class Warning(widget.OWWidget.Warning):
        no_expressions = Msg('Expression data objects not found.')
        no_data_objects = Msg('No expression data matches the selected options.')
        unexpected_feature_type = Msg('Can not import expression data, unexpected feature type "{}".')
        multiple_feature_type = Msg('Can not import expression data, multiple feature types found.')

    def __init__(self):
        super().__init__()
        ConcurrentWidgetMixin.__init__(self)

        self._res: Optional[ResolweAPI] = None

        # Store collection ID from currently selected row
        self.selected_collection_id: Optional[str] = None
        # Store data output options
        self.data_output_options: Optional[DataOutputOptions] = None
        # List of all available data objects (expressions)
        self.data_objects: Optional[List[Data]] = None
        # Subset of data objects based on DataOutputOptions
        self.filtered_data_objects: Optional[List[Data]] = None
        # List of orange data objects (clinical data)
        self.orange_data_objects: Dict[int, Data] = {}
        # Cache output data table
        self.data_table: Optional[Table] = None
        # Cache clinical metadata
        self.clinical_metadata: Optional[Table] = None

        # Control area
        self.info_box = gui.widgetLabel(gui.widgetBox(self.controlArea, "Info", margin=3), 'No data on output.')

        self.exp_type_combo = gui.comboBox(
            self.controlArea, self, 'exp_type', label='Expression Type', callback=self.on_output_option_changed
        )
        self.proc_type_combo = gui.comboBox(
            self.controlArea, self, 'proc_type', label='Process Name', callback=self.on_output_option_changed
        )
        self.input_anno_combo = gui.comboBox(
            self.controlArea,
            self,
            'input_annotation',
            label='Expression source',
            callback=self.on_output_option_changed,
        )
        self.orange_data_objects_combo = gui.comboBox(
            self.controlArea,
            self,
            'selected_orange_data_object',
            label='Clinical Data',
            callback=self.on_clinical_data_changed,
        )

        self.norm_component = NormalizationComponent(self, self.controlArea)
        self.norm_component.options_changed.connect(self.on_normalization_changed)

        gui.rubber(self.controlArea)
        box = gui.widgetBox(self.controlArea, 'Sign in')
        self.user_info = gui.label(box, self, '')
        self.server_info = gui.label(box, self, '')

        box = gui.widgetBox(box, orientation=Qt.Horizontal)
        self.sign_in_btn = gui.button(box, self, 'Sign In', callback=self.sign_in, autoDefault=False)
        self.sign_out_btn = gui.button(box, self, 'Sign Out', callback=self.sign_out, autoDefault=False)

        self.commit_button = gui.auto_commit(self.controlArea, self, 'auto_commit', '&Commit', box=False)
        self.commit_button.button.setAutoDefault(False)

        # Main area
        self.table_view = QTableView()
        self.table_view.setAlternatingRowColors(True)
        self.table_view.viewport().setMouseTracking(True)
        self.table_view.setShowGrid(False)
        self.table_view.verticalHeader().hide()
        self.table_view.horizontalHeader().setSectionResizeMode(QHeaderView.ResizeToContents)
        self.table_view.horizontalHeader().setStretchLastSection(True)
        self.table_view.setSelectionBehavior(QAbstractItemView.SelectRows)
        self.table_view.setSelectionMode(QAbstractItemView.SingleSelection)
        # self.table_view.setStyleSheet('QTableView::item:selected{background-color: palette(highlight); color: palette(highlightedText);};')

        self.model = GenialisExpressionsModel(self)
        self.model.setHorizontalHeaderLabels(TableHeader.labels())
        self.table_view.setModel(self.model)
        self.table_view.selectionModel().selectionChanged.connect(self.on_selection_changed)

        self.filter_component = CollapsibleFilterComponent(self, self.mainArea)
        self.filter_component.options_changed.connect(self.on_filter_changed)
        self.mainArea.layout().addWidget(self.table_view)
        self.pagination_component = PaginationComponent(self, self.mainArea)
        self.pagination_component.options_changed.connect(self.update_collections_view)

        self.sign_in(silent=True)

    @property
    def res(self):
        return self._res

    @res.setter
    def res(self, value: ResolweAPI):
        if isinstance(value, ResolweAPI):
            self._res = value
            self.update_user_status()
            self.update_collections_view()
            self.__invalidate()
            self.Outputs.table.send(None)

    def __invalidate(self):
        self.data_objects = None
        self.data_table = None
        self.selected_collection_id = None
        self.clinical_metadata = None
        self.orange_data_objects = {}

        self.data_output_options = None
        self.exp_type_combo.clear()
        self.proc_type_combo.clear()
        self.input_anno_combo.clear()
        self.orange_data_objects_combo.clear()

        self.Outputs.table.send(None)
        self.Warning.no_expressions.clear()
        self.Warning.multiple_feature_type.clear()
        self.Warning.unexpected_feature_type.clear()
        self.Warning.no_data_objects.clear()
        self.info.set_output_summary(StateInfo.NoOutput)
        self.update_info_box()

    def update_user_status(self):
        user = self.res.get_currently_logged_user()

        if user:
            user_info = f"{user[0].get('first_name', '')} {user[0].get('last_name', '')}".strip()
            user_info = f"User: {user_info if user_info else user[0].get('username', '')}"
            self.sign_in_btn.setEnabled(False)
            self.sign_out_btn.setEnabled(True)
        else:
            user_info = 'User: Anonymous'
            self.sign_in_btn.setEnabled(True)
            self.sign_out_btn.setEnabled(False)

        self.user_info.setText(user_info)
        self.server_info.setText(f'Server: {self.res.url[8:]}')

    def update_info_box(self):

        if self.data_table:
            total_genes = len(self.data_table.domain.attributes)
            known_genes = len([col for col in self.data_table.domain.attributes if len(col.attributes)])

            info_text = (
                '{} genes on output\n'
                '{} genes match Entrez database\n'
                '{} genes with match conflicts\n'.format(total_genes, known_genes, total_genes - known_genes)
            )

        else:
            info_text = 'No data on output.'

        self.info_box.setText(info_text)

    def sign_in(self, silent=False):
        dialog = SignInForm(self)

        if silent:
            dialog.sign_in()
            if dialog.resolwe_instance is not None:
                self.res = dialog.resolwe_instance
            else:
                self.res = connect(url=DEFAULT_URL)

        if not silent and dialog.exec_():
            self.res = dialog.resolwe_instance

    def sign_out(self):
        # Use public credentials when user signs out
        self.res = connect(url=DEFAULT_URL)
        # Remove username and
        cm = CredentialManager(CREDENTIAL_MANAGER_SERVICE)
        del cm.username
        del cm.password

    def on_filter_changed(self):
        self.pagination_component.reset_pagination()
        self.update_collections_view()

    def get_query_parameters(self) -> Dict[str, str]:
        params = {
            'limit': ItemsPerPage.values()[self.pagination_component.items_per_page],
            'offset': self.pagination_component.offset,
            'ordering': SortBy.values()[self.filter_component.sort_by],
        }

        if self.filter_component.filter_by_full_text:
            params.update({'text': self.filter_component.filter_by_full_text})

        if self.filter_component.filter_by_name:
            params.update({'name__icontains': self.filter_component.filter_by_name})

        if self.filter_component.filter_by_contrib:
            params.update({'contributor_name': self.filter_component.filter_by_contrib})

        if self.filter_component.filter_by_owner:
            params.update({'owners_name': self.filter_component.filter_by_owner})

        last_modified = FilterByDateModified.values()[self.filter_component.filter_by_modified]
        if last_modified:
            params.update({'modified__gte': last_modified.isoformat()})

        return params

    def get_collections(self) -> Tuple[Dict[str, str], Dict[str, str]]:
        # Get response from the server
        collections = self.res.get_collections(**self.get_query_parameters())
        # Loop trough collections and store ids
        collection_ids = [collection['id'] for collection in collections.get('results', [])]
        # Get species by collection ids
        collection_to_species = self.res.get_species(collection_ids)

        return collections, collection_to_species

    def update_collections_view(self):
        collections, collection_to_species = self.get_collections()

        # Pass the results to data model
        self.model.set_data(collections.get('results', []), collection_to_species)
        self.table_view.setItemDelegateForColumn(TableHeader.id, gui.LinkStyledItemDelegate(self.table_view))
        self.table_view.setColumnHidden(TableHeader.slug, True)
        self.table_view.setColumnHidden(TableHeader.tags, True)

        # Check pagination parameters and emit pagination_availability signal
        next_page = True if collections.get('next') else False
        previous_page = True if collections.get('previous') else False
        self.pagination_availability.emit(next_page, previous_page)

    def normalize(self, table: Table) -> Optional[Table]:
        if not table:
            return

        if self.norm_component.quantile_norm:
            table = QuantileNormalization()(table)

        if self.norm_component.log_norm:
            table = LogarithmicScale()(table)

        if self.norm_component.z_score_norm:
            table = ZScore(axis=self.norm_component.z_score_axis)(table)

        if self.norm_component.quantile_transform:
            axis = self.norm_component.quantile_transform_axis
            quantiles = table.X.shape[int(not axis)]
            distribution = QuantileTransformDist.values()[self.norm_component.quantile_transform_dist]
            table = QuantileTransform(axis=axis, n_quantiles=quantiles, output_distribution=distribution)(table)

        return table

    def commit(self):
        self.Warning.no_data_objects.clear()
        self.cancel()

        if not self.filtered_data_objects:
            self.Warning.no_data_objects()
            self.Outputs.table.send(None)
            self.info.set_output_summary(StateInfo.NoOutput)
            return

        self.start(self.runner)

    def filter_data_objects(self):
        options = self.data_output_options
        exp_type = options.expression[self.exp_type].type
        proc_type = options.process[self.proc_type].type
        source = options.input_annotation[self.input_annotation].source
        species = options.input_annotation[self.input_annotation].species
        build = options.input_annotation[self.input_annotation].build

        def is_selected(obj):
            return all(
                [
                    obj.process.type == proc_type,
                    obj.output['source'] == source,
                    obj.output['species'] == species,
                    obj.output['build'] == build,
                    obj.output['exp_type'] == exp_type if exp_type != 'rc' else True,
                ]
            )

        return [data_object for data_object in self.data_objects if is_selected(data_object)]

    def fetch_clinical_metadata(self):
        data_object = self.orange_data_objects.get(self.selected_orange_data_object, None)
        if data_object is not None:
            return self.res.get_clinical_metadata(data_object.id, data_object.output['table']['file'])

    def on_output_option_changed(self):
        self.data_table = None
        self.filtered_data_objects = self.filter_data_objects()
        self.commit()

    def on_clinical_data_changed(self):
        self.clinical_metadata = self.fetch_clinical_metadata()
        self.commit()

    def on_normalization_changed(self):
        self.commit()

    def on_selection_changed(self):
        self.__invalidate()

        collection_id: str = self.get_selected_row_data(TableHeader.id)
        if not collection_id:
            return

        self.selected_collection_id = collection_id
        self.data_objects = self.res.get_expression_data_objects(collection_id)
        self.orange_data_objects = {
            idx: data for idx, data in enumerate(self.res.get_orange_data_objects(collection_id))
        }

        self.data_output_options = available_data_output_options(self.data_objects)

        self.exp_type_combo.addItems(exp_name for _, exp_name in self.data_output_options.expression)
        if self.exp_type >= len(self.data_output_options.expression):
            self.exp_type = 0
        self.exp_type_combo.setCurrentIndex(self.exp_type)

        self.proc_type_combo.addItems(proc_name for _, proc_name in self.data_output_options.process)
        if self.proc_type >= len(self.data_output_options.process):
            self.proc_type = 0
        self.proc_type_combo.setCurrentIndex(self.proc_type)

        self.input_anno_combo.addItems(
            f'{species}, {build}' for _, species, build in self.data_output_options.input_annotation
        )
        if self.input_annotation >= len(self.data_output_options.input_annotation):
            self.input_annotation = 0
        self.input_anno_combo.setCurrentIndex(self.input_annotation)

        self.orange_data_objects_combo.addItems(data.name for data in self.orange_data_objects.values())
        self.clinical_metadata = self.fetch_clinical_metadata()

        if not self.data_objects:
            self.Warning.no_expressions()
            return

        # Note: This here is to handle an edge case where we get
        #       different 'feature_type' data object in a collection.
        #       For now we raise a warning, but in the future we should
        #       discuss about how to properly handle different types of features.
        feature_types = {data.output['feature_type'] for data in self.data_objects}

        if len(feature_types) == 1 and 'gene' not in feature_types:
            self.Warning.unexpected_feature_type(feature_types.pop())
            self.data_objects = []
            return

        if len(feature_types) > 1:
            self.Warning.multiple_feature_type()
            self.data_objects = []
            return

        self.on_output_option_changed()

    def get_selected_row_data(self, column: int) -> Optional[str]:
        selection_model = self.table_view.selectionModel()
        rows = selection_model.selectedRows(column=column)
        if not rows:
            return

        return rows[0].data()

    def on_done(self, table: Table):
        if table:
            samples, genes = table.X.shape
            self.info.set_output_summary(f'Samples: {samples} Genes: {genes}')
            self.update_info_box()
            self.Outputs.table.send(table)

    def on_exception(self, ex):
        if isinstance(ex, ResolweDataObjectsNotFound):
            self.Warning.no_data_objects()
            self.Outputs.table.send(None)
            self.data_table = None
            self.info.set_output_summary(StateInfo.NoOutput)
            self.update_info_box()
        else:
            raise ex

    def on_partial_result(self, result: Any) -> None:
        pass

    def onDeleteWidget(self):
        self.shutdown()
        super().onDeleteWidget()

    def sizeHint(self):
        return QSize(1280, 620)

    def runner(self, state: TaskState) -> Table:
        data_objects = self.filtered_data_objects
        exp_type = self.data_output_options.expression[self.exp_type].type
        species = self.data_output_options.input_annotation[self.input_annotation].species
        output_field = 'exp' if exp_type != 'rc' else exp_type
        collection_id = self.selected_collection_id

        table = self.data_table
        step, steps = 0, 4

        def set_progress():
            nonlocal step
            step += 1
            state.set_progress_value(100 * (step / steps))

        if not table:
            steps += len(data_objects)
            state.set_status('Downloading ...')
            data_frames = []

            for data_object in data_objects:
                set_progress()

                response = self.res.download(data_object.id, data_object.output[output_field]['file'])
                with io.BytesIO() as f:
                    f.write(response.content)
                    f.seek(0)
                    # expressions to data frame
                    df = pd.read_csv(f, sep='\t', compression='gzip')
                    df = df.set_index('Gene').T.reset_index(drop=True)
                    data_frames.append(df)

            state.set_status('Concatenating samples ...')
            df = pd.concat(data_frames, axis=0)

            state.set_status('To data table ...')
            table = table_from_frame(df)
            set_progress()

            state.set_status('Matching genes ...')
            tax_id = species_name_to_taxid(species)
            gm = GeneMatcher(tax_id)
            table = gm.match_table_attributes(table, rename=True)
            table.attributes[TableAnnotation.tax_id] = tax_id
            table.attributes[TableAnnotation.gene_as_attr_name] = True
            table.attributes[TableAnnotation.gene_id_attribute] = 'Entrez ID'

            self.data_table = table
            set_progress()

        state.set_status('Adding metadata ...')
        metadata = self.res.aggregate_sample_metadata([data.sample for data in data_objects], collection_id)
        domain = Domain(table.domain.attributes, table.domain.class_vars, list(metadata.keys()))
        table = table.transform(domain)
        for key, value in metadata.items():
            table[:, key] = np.reshape(value, (-1, 1))

        if self.clinical_metadata:
            table = merge_clinical_data(table, self.clinical_metadata)
        set_progress()

        state.set_status('Normalizing ...')
        table = self.normalize(table)
        set_progress()

        return table


if __name__ == "__main__":
    from orangewidget.utils.widgetpreview import WidgetPreview

    WidgetPreview(OWGenialisExpressions).run()
