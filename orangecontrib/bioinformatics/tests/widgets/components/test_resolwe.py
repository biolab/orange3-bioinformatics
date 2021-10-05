import unittest
from unittest.mock import patch

from AnyQt.QtTest import QSignalSpy

from Orange.widgets.tests.base import GuiTest, WidgetTest
from Orange.widgets.credentials import CredentialManager
from Orange.widgets.tests.utils import simulate

from orangecontrib.bioinformatics.resolwe import ResolweAuthError
from orangecontrib.bioinformatics.widgets.components.resolwe import SignIn
from orangecontrib.bioinformatics.widgets.OWGenialisExpressions import SortBy, ItemsPerPage, FilterByDateModified
from orangecontrib.bioinformatics.tests.widgets.test_OWGenialisExpressions import MockWidget


class TestSignInForm(GuiTest):
    @patch(
        'orangecontrib.bioinformatics.widgets.components.resolwe.get_credential_manager',
        return_value=CredentialManager('resolwe_credentials_test'),
    )
    @patch('orangecontrib.bioinformatics.widgets.components.resolwe.connect')
    def test_dialog_success(self, mocked_connect, mocked_cm):
        widget = MockWidget()
        dialog = SignIn(widget)
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
        if cm.username:
            del cm.username
        if cm.password:
            del cm.password

    @patch(
        'orangecontrib.bioinformatics.widgets.components.resolwe.get_credential_manager',
        return_value=CredentialManager('resolwe_credentials_test'),
    )
    @patch('orangecontrib.bioinformatics.widgets.components.resolwe.connect', side_effect=ResolweAuthError())
    def test_dialog_fail(self, mocked_connect, mocked_cm):
        widget = MockWidget()
        dialog = SignIn(widget)
        dialog.username_line_edit.setText('foo')
        dialog.password_line_edit.setText('bar')

        self.assertTrue(dialog.error_msg.isHidden())
        dialog.sign_in()
        self.assertFalse(dialog.error_msg.isHidden())

        mocked_connect.assert_called_once()


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


if __name__ == "__main__":
    unittest.main()
