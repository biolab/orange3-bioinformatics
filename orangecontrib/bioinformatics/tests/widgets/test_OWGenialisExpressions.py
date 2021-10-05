import random
import unittest
from datetime import datetime, timezone
from unittest.mock import call, patch

from AnyQt.QtCore import pyqtSignal

from Orange.widgets.widget import OWWidget
from Orange.widgets.settings import SettingProvider
from Orange.widgets.tests.base import WidgetTest
from Orange.widgets.credentials import CredentialManager

from orangecontrib.bioinformatics.resolwe import resapi
from orangecontrib.bioinformatics.widgets.components.resolwe import SignIn
from orangecontrib.bioinformatics.widgets.OWGenialisExpressions import (
    SortBy,
    ItemsPerPage,
    PaginationComponent,
    FilterByDateModified,
    OWGenialisExpressions,
    CollapsibleFilterComponent,
)


class MockWidget(OWWidget):
    name = 'MockWidget'
    want_main_area = False

    filter_component = SettingProvider(CollapsibleFilterComponent)
    pagination_component = SettingProvider(PaginationComponent)

    pagination_availability = pyqtSignal(bool, bool)

    @patch(
        'orangecontrib.bioinformatics.widgets.components.resolwe.get_credential_manager',
        return_value=CredentialManager('resolwe_credentials_test'),
    )
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.filter_component = CollapsibleFilterComponent(self, self.controlArea)
        self.pagination_component = PaginationComponent(self, self.controlArea)
        self.sign_in_dialog = SignIn(self)


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
        mock_resolwe_api.assert_called_once_with(None, None, resapi.DEFAULT_URL)
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
