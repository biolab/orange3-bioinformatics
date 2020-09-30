import random
import unittest
from datetime import datetime, timezone
from unittest.mock import call, patch

from orangewidget.tests.base import GuiTest, WidgetTest

from AnyQt.QtCore import pyqtSignal
from AnyQt.QtTest import QSignalSpy

from Orange.widgets.widget import OWWidget
from Orange.widgets.settings import SettingProvider
from Orange.widgets.credentials import CredentialManager
from Orange.widgets.tests.utils import simulate

from orangecontrib.bioinformatics.widgets.OWGenialisExpressions import (
    DEFAULT_URL,
    SortBy,
    SignInForm,
    ItemsPerPage,
    PaginationComponent,
    FilterByDateModified,
    ResolweAuthException,
    OWGenialisExpressions,
    CollapsibleFilterComponent,
)


class MockWidget(OWWidget):
    name = 'MockWidget'
    want_main_area = False

    filter_component = SettingProvider(CollapsibleFilterComponent)
    pagination_component = SettingProvider(PaginationComponent)

    pagination_availability = pyqtSignal(bool, bool)

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.filter_component = CollapsibleFilterComponent(self, self.controlArea)
        self.pagination_component = PaginationComponent(self, self.controlArea)
        self.sign_in_dialog = SignInForm(self)


class TestCollapsibleFilterComponent(WidgetTest):
    def setUp(self):
        self.widget = MockWidget()
        self.component = self.widget.filter_component
        self.emitted_signals: QSignalSpy = QSignalSpy(self.component.options_changed)

    def test_filter_by_modified(self):
        # check if correct option is set by default
        self.assertEqual(FilterByDateModified.any_time, self.component.filter_by_modified)
        # run trough all the options
        simulate.combobox_run_through_all(self.component.filter_modified)
        self.assertEqual(FilterByDateModified.past_month, self.component.filter_by_modified)

        # check if signal is emitted on every option changed
        self.assertEqual(len(FilterByDateModified.labels()), len(self.emitted_signals))

    def test_sorting(self):
        # check if correct option is set by default
        self.assertEqual(SortBy.newest_first, self.component.sort_by)
        # run trough all the options
        simulate.combobox_run_through_all(self.component.sorting)
        self.assertEqual(SortBy.oldest_first, self.component.sort_by)

        # check if signal is emitted on every option changed
        self.assertEqual(len(SortBy.labels()), len(self.emitted_signals))

    def test_full_text_search(self):
        search_parameter = 'foobar'
        # Default sorting should be set initially
        self.assertEqual(self.component.sort_by, SortBy.newest_first)
        # simulate user input
        self.component.filter_full_text.setText(search_parameter)
        self.component.on_filter_full_text_changed()
        # check if setting is correctly set
        self.assertEqual(self.component.filter_by_full_text, search_parameter)
        # check if sorting is set accordingly
        self.assertEqual(SortBy.relevance, self.component.sort_by)


class TestPaginationComponent(WidgetTest):
    def setUp(self):
        self.widget = MockWidget()
        self.component = self.widget.pagination_component

    def test_pagination_availability(self):
        self.widget.pagination_availability.emit(False, False)
        self.assertFalse(self.component.page_left_btn.isEnabled())
        self.assertFalse(self.component.page_right_btn.isEnabled())

        self.widget.pagination_availability.emit(True, False)
        self.assertFalse(self.component.page_left_btn.isEnabled())
        self.assertTrue(self.component.page_right_btn.isEnabled())

        self.widget.pagination_availability.emit(False, True)
        self.assertTrue(self.component.page_left_btn.isEnabled())
        self.assertFalse(self.component.page_right_btn.isEnabled())

        self.widget.pagination_availability.emit(True, True)
        self.assertTrue(self.component.page_left_btn.isEnabled())
        self.assertTrue(self.component.page_right_btn.isEnabled())

    def test_page_limit(self):
        # test default option
        self.assertTrue(self.component.page_limit.buttons[ItemsPerPage.min].isChecked())
        self.assertFalse(self.component.page_limit.buttons[ItemsPerPage.med].isChecked())
        self.assertFalse(self.component.page_limit.buttons[ItemsPerPage.max].isChecked())
        self.assertEqual(0, self.component.offset)
        self.assertEqual(1, self.component.current_page)

        # change current page
        self.component.right_btn_pressed()
        self.assertEqual(ItemsPerPage.values()[ItemsPerPage.min], self.component.offset)
        self.assertEqual(2, self.component.current_page)

        # change page limit
        self.component.items_per_page = ItemsPerPage.max
        self.component.on_limit_changed()
        self.assertEqual(0, self.component.offset)
        self.assertEqual(1, self.component.current_page)
        self.assertFalse(self.component.page_limit.buttons[ItemsPerPage.min].isChecked())
        self.assertFalse(self.component.page_limit.buttons[ItemsPerPage.med].isChecked())
        self.assertTrue(self.component.page_limit.buttons[ItemsPerPage.max].isChecked())

        # change current page
        self.component.right_btn_pressed()
        self.assertEqual(ItemsPerPage.values()[ItemsPerPage.max], self.component.offset)
        self.assertEqual(2, self.component.current_page)

    def test_pagination(self):
        # change page limit
        self.component.items_per_page = ItemsPerPage.med
        self.component.on_limit_changed()
        self.assertEqual(0, self.component.offset)
        self.assertEqual(1, self.component.current_page)
        self.assertFalse(self.component.page_limit.buttons[ItemsPerPage.min].isChecked())
        self.assertTrue(self.component.page_limit.buttons[ItemsPerPage.med].isChecked())
        self.assertFalse(self.component.page_limit.buttons[ItemsPerPage.max].isChecked())

        current_offset = 0

        for page in [2, 3, 4]:
            self.component.right_btn_pressed()
            current_offset += ItemsPerPage.values()[ItemsPerPage.med]
            self.assertEqual(current_offset, self.component.offset)
            self.assertEqual(page, self.component.current_page)

        for page in [3, 2, 1]:
            self.component.left_btn_pressed()
            current_offset -= ItemsPerPage.values()[ItemsPerPage.med]
            self.assertEqual(current_offset, self.component.offset)
            self.assertEqual(page, self.component.current_page)

        # test left pagination limit
        for _ in range(3):
            self.component.left_btn_pressed()
            self.assertEqual(0, self.component.offset)
            self.assertEqual(1, self.component.current_page)


class TestSignInForm(GuiTest):
    @patch(
        'orangecontrib.bioinformatics.widgets.OWGenialisExpressions.CREDENTIAL_MANAGER_SERVICE',
        'resolwe_credentials_test',
    )
    @patch('orangecontrib.bioinformatics.widgets.OWGenialisExpressions.connect')
    def test_dialog_success(self, mocked_connect):
        widget = MockWidget()
        dialog = SignInForm(widget)
        self.assertFalse(dialog.sign_in_btn.isEnabled())
        dialog.sign_in()
        mocked_connect.assert_not_called()
        self.assertIsNone(dialog.resolwe_instance)

        dialog.username_line_edit.setText('foo')
        dialog.sign_in()
        mocked_connect.assert_not_called()

        dialog.password_line_edit.setText('bar')
        self.assertTrue(dialog.sign_in_btn.isEnabled())
        dialog.sign_in()

        mocked_connect.assert_called_once()
        self.assertTrue(dialog.error_msg.isHidden())
        self.assertIsNotNone(dialog.resolwe_instance)

        # cleanup
        cm = CredentialManager('resolwe_credentials_test')
        del cm.username
        del cm.password

    @patch(
        'orangecontrib.bioinformatics.widgets.OWGenialisExpressions.CREDENTIAL_MANAGER_SERVICE',
        'resolwe_credentials_test',
    )
    @patch('orangecontrib.bioinformatics.widgets.OWGenialisExpressions.connect', side_effect=ResolweAuthException())
    def test_dialog_fail(self, mocked_connect):
        widget = MockWidget()
        dialog = SignInForm(widget)
        dialog.username_line_edit.setText('foo')
        dialog.password_line_edit.setText('bar')

        self.assertTrue(dialog.error_msg.isHidden())
        dialog.sign_in()
        self.assertFalse(dialog.error_msg.isHidden())

        mocked_connect.assert_called_once()

        # cleanup
        cm = CredentialManager('resolwe_credentials_test')
        del cm.username
        del cm.password


def collection_generator(num_of_collections):
    return [
        {
            'id': int(i),
            'name': f'Collection name {i}',
            'contributor': {'first_name': f'Foo_{i}', 'last_name': f'Bar_{i}', 'username': f'FooBar_{i}'},
            'description': f'Test description {1}',
            'data_count': random.randint(0, 10),
            'entity_count': random.randint(0, 10),
            'tags': ['community:expressions'],
            'created': datetime.strftime(datetime.now().replace(tzinfo=timezone.utc), '%Y-%m-%dT%H:%M:%S.%f%z'),
            'modified': datetime.strftime(datetime.now().replace(tzinfo=timezone.utc), '%Y-%m-%dT%H:%M:%S.%f%z'),
        }
        for i in range(num_of_collections)
    ]


class TestOWGenialisExpressions(WidgetTest):
    @patch('orangecontrib.bioinformatics.resolwe.ResolweAPI', spec=True)
    @patch('orangecontrib.bioinformatics.widgets.OWGenialisExpressions.OWGenialisExpressions.sign_in')
    def test_signin_popup(self, mock_sign_in, _):
        self.widget = self.create_widget(OWGenialisExpressions)
        self.widget.sign_in_btn.click()
        mock_sign_in.assert_has_calls([call(silent=True), call(False)])

    @patch('orangecontrib.bioinformatics.resolwe.ResolweAPI', spec=True)
    def test_widget_initialization(self, mock_resolwe_api):
        count = 2
        mock_resolwe_api.return_value.get_collections.return_value = {
            'count': count,
            'next': None,
            'previous': None,
            'results': collection_generator(count),
        }
        self.widget = self.create_widget(OWGenialisExpressions)

        # Check if calls are made with defaults
        offset = 0
        limit = ItemsPerPage.values()[ItemsPerPage.min]
        ordering = SortBy.values()[SortBy.newest_first]
        mock_resolwe_api.assert_called_once_with(None, None, DEFAULT_URL)
        mock_resolwe_api.return_value.get_collections.assert_called_once_with(
            limit=limit, offset=offset, ordering=ordering
        )
        mock_resolwe_api.return_value.get_species.assert_called_once_with(list(range(count)))

        self.assertEqual(count, self.widget.model.rowCount())

    @patch('orangecontrib.bioinformatics.resolwe.ResolweAPI', spec=True)
    def test_sign_out(self, mock_resolwe_api):
        # collections before/after sign out
        count_before = 10
        count_after = 2

        mock_resolwe_api.return_value.get_collections.return_value = {
            'count': count_before,
            'next': None,
            'previous': None,
            'results': collection_generator(count_before),
        }
        self.widget = self.create_widget(OWGenialisExpressions)
        self.assertEqual(count_before, self.widget.model.rowCount())

        # check if state resets after sign out
        mock_resolwe_api.return_value.get_collections.return_value = {
            'count': count_after,
            'next': None,
            'previous': None,
            'results': collection_generator(count_after),
        }
        self.widget.sign_out_btn.click()
        self.assertEqual(count_after, self.widget.model.rowCount())

    @patch('orangecontrib.bioinformatics.resolwe.ResolweAPI', spec=True)
    def test_query_params(self, mock_resolwe_api):
        self.widget = self.create_widget(OWGenialisExpressions)

        full_text = 'test_fts'
        self.widget.filter_component.filter_full_text.setText(full_text)

        col_name = 'test_name'
        self.widget.filter_component.filter_name.setText(col_name)

        col_owner = 'test_owner'
        self.widget.filter_component.filter_owner.setText(col_owner)

        col_contrib = 'test_contrib'
        self.widget.filter_component.filter_contrib.setText(col_contrib)

        self.widget.filter_component.filter_by_modified = FilterByDateModified.past_hour

        self.widget.on_filter_changed()
        _, query_params = mock_resolwe_api.return_value.get_collections.call_args
        self.assertIn('limit', query_params)
        self.assertIn('offset', query_params)
        self.assertIn('ordering', query_params)
        self.assertIn('ordering', query_params)
        self.assertIn('text', query_params)
        self.assertIn('name__icontains', query_params)
        self.assertIn('contributor_name', query_params)
        self.assertIn('owners_name', query_params)
        self.assertIn('modified__gte', query_params)
        self.assertIn(full_text, query_params.values())
        self.assertIn(col_name, query_params.values())
        self.assertIn(col_owner, query_params.values())
        self.assertIn(col_contrib, query_params.values())


if __name__ == "__main__":
    unittest.main()
