""" ResolweAPI """
from typing import Dict, List
from urllib.parse import urljoin

import requests
from resdk import Resolwe
from resdk.resources.data import Data, Collection

DEFAULT_URL: str = 'https://app.genialis.com'
RESOLWE_URLS: List[str] = [DEFAULT_URL, 'https://imaps.genialis.com', 'https://bcm.genialis.com']
CREDENTIAL_MANAGER_SERVICE: str = 'resolwe_credentials'


class ResolweAPI:

    COLLECTION_FIELDS = (
        'id',
        'slug',
        'name',
        'created',
        'modified',
        'description',
        'data_count',
        'entity_count',
        'contributor__first_name',
        'contributor__last_name',
        'contributor__username',
        'tags',
    )

    DATA_FIELDS = [
        'output__source',
        'output__species',
        'output__build',
        'output__exp_type',
        'output__feature_type',
        'process__slug',
        'process__name',
    ]

    def __init__(self, username=None, password=None, url=DEFAULT_URL):
        self._res = Resolwe(username, password, url)

    @property
    def url(self):
        return self._res.url

    @property
    def base_url(self) -> str:
        return self._res.api._store['base_url']

    @property
    def collection_url(self):
        return f'{self.url}/expressions/search/collection/'

    @property
    def session(self) -> requests.Session:
        return self._res.session

    def get_currently_logged_user(self):
        url = urljoin(self.base_url, 'user?current_only=1')
        response = self.session.get(url, auth=self._res.auth)
        response.raise_for_status()
        return response.json()

    def get_collections(self, **kwargs):
        kwargs.update({'fields': ','.join(self.COLLECTION_FIELDS)})
        params = '&'.join([f'{key}={value}' for key, value in kwargs.items() if value is not None])
        url = urljoin(self.base_url, f'collection?{params}')

        response = self.session.get(url, auth=self._res.auth)
        response.raise_for_status()
        return response.json()

    def get_species(self, collection_ids: List[str]) -> Dict[str, str]:
        params = '&'.join(f'collection[]={col_id}' for col_id in collection_ids)
        url = urljoin(self.base_url, f'_modules/species/collection?{params}')

        response = self.session.get(url, auth=self._res.auth)
        response.raise_for_status()
        response_data = response.json()

        col_to_species = {
            data['collection_id']: data['species'][0] if len(data['species']) == 1 else '' for data in response_data
        }

        return col_to_species

    def get_collection_by_id(self, collection_id: str) -> Collection:
        return self._res.collection.get(id=collection_id)

    def get_expression_data_objects(self, collection_id: str) -> List[Data]:
        """ Get all expression data objects from a collection. """
        return self._res.data.filter(type='data:expression', collection=collection_id, fields=self.DATA_FIELDS)


if __name__ == '__main__':
    res = ResolweAPI(url=DEFAULT_URL)
    print(res)
