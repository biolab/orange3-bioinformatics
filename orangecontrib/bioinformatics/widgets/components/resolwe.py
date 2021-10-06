from enum import IntEnum
from datetime import datetime

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
    QHBoxLayout,
    QPushButton,
    QScrollArea,
    QSizePolicy,
    QToolButton,
    QVBoxLayout,
)

from Orange.widgets import gui, settings
from Orange.widgets.widget import OWComponent
from Orange.widgets.credentials import CredentialManager

from orangecontrib.bioinformatics.resolwe import (
    GENESIS_PLATFORM,
    RESOLWE_PLATFORM,
    ResolweAuthError,
    genapi,
    resapi,
    connect,
)


def get_credential_manager(server_type: str) -> CredentialManager:
    if server_type == RESOLWE_PLATFORM:
        service_name = resapi.CREDENTIAL_MANAGER_SERVICE
    elif server_type == GENESIS_PLATFORM:
        service_name = genapi.CREDENTIAL_MANAGER_SERVICE
    else:
        raise ValueError('Unexpected server type. Available options: resolwe or genesis.')

    return CredentialManager(service_name)


class SignIn(QDialog):
    def __init__(self, flags, *args, server_type: str = 'resolwe', **kwargs):
        super().__init__(flags, *args, **kwargs)

        self.server_type: str = server_type
        self.cm: CredentialManager = get_credential_manager(server_type)

        self.setWindowTitle('Sign in')
        self.setFixedSize(400, 250)

        self.server_cb_label = QLabel('Server *')
        self.server_cb = QComboBox(self)
        self.server_cb.addItems(resapi.RESOLWE_URLS if server_type == RESOLWE_PLATFORM else [genapi.DEFAULT_URL])
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
            self.resolwe_instance = connect(username, password, url=server, server_type=self.server_type)
        except ResolweAuthError:
            self.error_msg.show()
            return

        self.cm.username = username
        self.cm.password = password
        self.accept()


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
            left_box, self, 'filter_by_name', label=self.FILTER_NAME_LABEL, callback=self.on_filter_changed
        )
        self.filter_contrib = gui.lineEdit(
            mid_box,
            self,
            'filter_by_contrib',
            label=self.FILTER_CONTRIB_LABEL,
            callback=self.on_filter_changed,
        )
        self.filter_owner = gui.lineEdit(
            right_box,
            self,
            'filter_by_owner',
            label=self.FILTER_OWNER_LABEL,
            callback=self.on_filter_changed,
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
        """Start animation"""
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
    """Gene expression normalization component"""

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
