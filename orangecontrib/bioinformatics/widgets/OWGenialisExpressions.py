import io
from enum import IntEnum
from typing import Any, Dict, List, Optional
from numbers import Number
from datetime import datetime
from collections import defaultdict

import numpy as np
import pandas as pd
from scipy.stats import zscore
from resdk.resources.data import Data
from dateutil.relativedelta import relativedelta
from resdk.resources.sample import Sample

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

from Orange.data import Table, Domain, StringVariable
from Orange.widgets import gui, widget, settings
from Orange.widgets.widget import Msg, Output, StateInfo, OWComponent
from Orange.data.pandas_compat import table_from_frame
from Orange.widgets.utils.concurrent import TaskState, ConcurrentWidgetMixin
from Orange.widgets.utils.itemmodels import PyTableModel

from orangecontrib.bioinformatics.resolwe import ResolweAPI, ResolweAuthException, connect
from orangecontrib.bioinformatics.ncbi.gene import GeneMatcher
from orangecontrib.bioinformatics.ncbi.taxonomy import species_name_to_taxid
from orangecontrib.bioinformatics.resolwe.resapi import DEFAULT_URL, RESOLWE_URLS, SAMPLE_DESCRIPTOR_LABELS
from orangecontrib.bioinformatics.widgets.utils.data import TableAnnotation

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

    def _reset_pagination(self):
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
        self._reset_pagination()
        self.options_changed.emit()


class SignInForm(QDialog):
    def __init__(self, flags, *args, **kwargs):
        super().__init__(flags, *args, **kwargs)
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
        username = self.username_line_edit.text()
        password = self.password_line_edit.text()

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

        elif role == Qt.TextAlignmentRole and isinstance(value, Number):
            return Qt.AlignRight | Qt.AlignVCenter

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


def runner(
    res: ResolweAPI,
    data_objects: List[Data],
    exp_type: str,
    species: str,
    state: TaskState,
    log_norm=False,
    z_score_norm=False,
) -> Table:
    data_frames = []
    metadata = defaultdict(list)

    def parse_sample_descriptor(sample: Sample) -> None:
        general = sample.descriptor.get('general', {})

        for label in SAMPLE_DESCRIPTOR_LABELS:
            metadata[label].append([general.get(label, '')])

        metadata['sample_name'].append([sample.name])

    if exp_type != 'rc':
        output_field = 'exp'
        data_objects = [obj for obj in data_objects if obj.output['exp_type'] == exp_type]
    else:
        output_field = 'rc'

    step, steps = 0, len(data_objects) + 3

    def set_progress():
        nonlocal step
        step += 1
        state.set_progress_value(100 * (step / steps))

    state.set_status('Downloading ...')
    for data_object in data_objects:
        set_progress()
        parse_sample_descriptor(data_object.sample)
        metadata['expression_type'].append([exp_type.upper()])

        response = res.get_expressions(data_object.id, data_object.output[output_field]['file'])
        with io.BytesIO() as f:
            f.write(response.content)
            f.seek(0)
            # expressions to data frame
            df = pd.read_csv(f, sep='\t', compression='gzip')
            df = df.set_index('Gene').T.reset_index(drop=True)
            data_frames.append(df)

    state.set_status('Concatenating samples ...')
    df = pd.concat(data_frames, axis=0)

    if log_norm:
        state.set_status('Applying log2 normalization ...')
        df = np.log2(df + 1)

    if z_score_norm:
        state.set_status('Applying z-score normalization ...')
        columns = df.select_dtypes(include=[np.number]).columns
        df[columns] = zscore(df[columns], axis=1)

    state.set_status('To data table ...')
    table = table_from_frame(df)
    set_progress()

    state.set_status('Adding metadata ...')
    metas = [StringVariable(label) for label in metadata.keys()]
    domain = Domain(table.domain.attributes, table.domain.class_vars, metas)
    table = table.transform(domain)

    for key, value in metadata.items():
        table[:, key] = value
    set_progress()

    state.set_status('Matching genes ...')
    tax_id = species_name_to_taxid(species)
    gm = GeneMatcher(tax_id)
    gm.match_table_attributes(table)
    table.attributes[TableAnnotation.tax_id] = tax_id
    table.attributes[TableAnnotation.gene_as_attr_name] = True
    table.attributes[TableAnnotation.gene_id_attribute] = 'Entrez ID'
    set_progress()

    return table


class OWMGenialisExpressions(widget.OWWidget, ConcurrentWidgetMixin):
    name = 'Genialis Expressions'
    priority = 180
    want_main_area = True
    want_control_area = True
    icon = '../widgets/icons/OWGenialisExpressions.svg'

    pagination_availability = pyqtSignal(bool, bool)

    pagination_component = settings.SettingProvider(PaginationComponent)
    filter_component = settings.SettingProvider(CollapsibleFilterComponent)

    exp_type: int
    exp_type = settings.Setting(None, schema_only=True)

    log_norm: bool
    log_norm = settings.Setting(False, schema_only=True)

    z_score_norm: bool
    z_score_norm = settings.Setting(False, schema_only=True)

    auto_commit: bool
    auto_commit = settings.Setting(False, schema_only=True)

    class Outputs:
        table = Output('Expressions', Table)

    class Warning(widget.OWWidget.Warning):
        no_expressions = Msg('Expression data objects not found.')

    def __init__(self):
        super().__init__()
        ConcurrentWidgetMixin.__init__(self)
        self._res = connect(url=DEFAULT_URL)

        # Control area
        box = gui.widgetBox(self.controlArea, 'Sign in')
        self.user_info = gui.label(box, self, '')
        self.server_info = gui.label(box, self, '')

        box = gui.widgetBox(box, orientation=Qt.Horizontal)
        self.sign_in_btn = gui.button(box, self, 'Sign in', callback=self.sign_in, autoDefault=False)
        self.sign_out_btn = gui.button(box, self, 'Sign out', callback=self.sign_out, autoDefault=False)

        self.exp_type_box = gui.widgetBox(self.controlArea, 'Expression type')
        self.exp_type_options = gui.radioButtons(
            self.exp_type_box, self, 'exp_type', callback=self.on_exp_type_changed
        )
        self.set_exp_type_options()

        box = gui.widgetBox(self.controlArea, 'Normalization')
        gui.checkBox(box, self, 'log_norm', 'log2(x+1)', callback=self.on_normalization_changed)
        gui.checkBox(box, self, 'z_score_norm', 'z-score', callback=self.on_normalization_changed)

        gui.rubber(self.controlArea)
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
        self.filter_component.options_changed.connect(self.update_collections_view)
        self.mainArea.layout().addWidget(self.table_view)
        self.pagination_component = PaginationComponent(self, self.mainArea)
        self.pagination_component.options_changed.connect(self.update_collections_view)

        self.data_objects: Optional[List[Data]] = None
        self.update_collections_view()
        self.update_user_status()

    def __invalidate(self):
        self.data_objects = None
        self.Warning.no_expressions.clear()
        self.info.set_output_summary(StateInfo.NoOutput)

    def set_exp_type_options(self, options: List[str] = None) -> None:
        self.clear_exp_type_options()
        gui.appendRadioButton(self.exp_type_options, 'Read Counts', disabled=not bool(options))

        for option in options if options is not None else []:
            gui.appendRadioButton(self.exp_type_options, option)

        if len(self.exp_type_options.buttons) > 1:
            self.exp_type = 1

    def clear_exp_type_options(self) -> None:
        for btn in self.exp_type_options.buttons:
            btn.deleteLater()
        self.exp_type_options.buttons = []

    @property
    def res(self):
        return self._res

    @res.setter
    def res(self, value):
        if isinstance(value, ResolweAPI):
            self._res = value
            self.update_user_status()
            self.update_collections_view()
            self.__invalidate()
            self.Outputs.table.send(None)

    def update_user_status(self):
        user = self.res.get_currently_logged_user()

        if user:
            user_info = f"{user[0].get('first_name', '')} {user[0].get('last_name', '')}".strip()
            user_info = f"User: {user_info if user_info else user[0].get('username', '')}"
            self.sign_out_btn.setEnabled(True)
        else:
            user_info = 'User: Anonymous'
            self.sign_out_btn.setEnabled(False)

        self.user_info.setText(user_info)
        self.server_info.setText(f'Server: {self.res.url[8:]}')

    def sign_in(self):
        dialog = SignInForm(self)
        if dialog.exec_():
            self.res = dialog.resolwe_instance

    def sign_out(self):
        self.res = connect(url=DEFAULT_URL)

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

    def update_collections_view(self):
        # Get response from the server
        collections = self.res.get_collections(**self.get_query_parameters())
        # Loop trough collections and store ids
        collection_ids = [collection['id'] for collection in collections.get('results', [])]
        # Get species by collection ids
        collection_to_species = self.res.get_species(collection_ids)

        # Pass the results to data model
        self.model.set_data(collections.get('results', []), collection_to_species)
        self.table_view.setItemDelegateForColumn(TableHeader.id, gui.LinkStyledItemDelegate(self.table_view))
        self.table_view.setColumnHidden(TableHeader.slug, True)
        self.table_view.setColumnHidden(TableHeader.tags, True)

        # Check pagination parameters and emit pagination_availability signal
        next_page = True if collections.get('next') else False
        previous_page = True if collections.get('previous') else False
        self.pagination_availability.emit(next_page, previous_page)

    def commit(self):
        self.cancel()

        if self.data_objects:
            exp_type: str = 'rc' if self.exp_type == 0 else self.exp_type_options.group.checkedButton().text()
            collection_species: str = self.get_selected_row_data(TableHeader.species)

            self.start(
                runner,
                self.res,
                self.data_objects,
                exp_type,
                collection_species,
                log_norm=self.log_norm,
                z_score_norm=self.z_score_norm,
            )

    def on_exp_type_changed(self):
        if self.data_objects:
            self.commit()

    def on_normalization_changed(self):
        if self.data_objects:
            self.commit()

    def on_selection_changed(self):
        self.__invalidate()

        collection_id: str = self.get_selected_row_data(TableHeader.id)
        if not collection_id:
            return

        self.data_objects = self.res.get_expression_data_objects(collection_id)
        self.set_exp_type_options(list({data.output['exp_type'] for data in self.data_objects}))

        if not self.data_objects:
            self.Warning.no_expressions()
            return

        self.commit()

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
            self.Outputs.table.send(table)

    def on_partial_result(self, result: Any) -> None:
        pass

    def onDeleteWidget(self):
        self.shutdown()
        super().onDeleteWidget()


if __name__ == "__main__":
    from orangewidget.utils.widgetpreview import WidgetPreview

    WidgetPreview(OWMGenialisExpressions).run()
