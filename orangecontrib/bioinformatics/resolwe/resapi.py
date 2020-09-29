""" ResolweAPI """
from typing import Dict, List
from datetime import timedelta
from urllib.parse import urljoin

import requests
import requests_cache
from resdk import Resolwe
from requests import Response
from resdk.resources.data import Data

from orangecontrib.bioinformatics.resolwe.utils import CACHE_BACKEND, RESOLWEAPI_CACHE

DEFAULT_URL: str = 'https://app.genialis.com'
RESOLWE_URLS: List[str] = [DEFAULT_URL, 'https://imaps.genialis.com', 'https://bcm.genialis.com']
SAMPLE_DESCRIPTOR_LABELS: List[str] = ['species', 'genotype', 'biosample_source', 'biosample_treatment']
CREDENTIAL_MANAGER_SERVICE: str = 'resolwe_credentials'

expire_after = timedelta(days=365)
requests_cache.install_cache(cache_name=RESOLWEAPI_CACHE, backend=CACHE_BACKEND, expire_after=expire_after)


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

        with self.session.cache_disabled():
            response = self.session.get(url, auth=self._res.auth)
            response.raise_for_status()
            return response.json()

    def get_collections(self, **kwargs):
        kwargs.update({'fields': ','.join(self.COLLECTION_FIELDS)})
        params = '&'.join([f'{key}={value}' for key, value in kwargs.items() if value is not None])
        url = urljoin(self.base_url, f'collection?{params}')

        with self.session.cache_disabled():
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

    def get_expressions(self, data_id, data_file_name) -> Response:
        file_url = urljoin(self.url, f'data/{data_id}/{data_file_name}')
        response = self.session.get(file_url, auth=self._res.auth)
        response.raise_for_status()
        return response

    def get_expression_data_objects(self, collection_id: str) -> List[Data]:
        return self._res.data.filter(type='data:expression', collection=collection_id)


if __name__ == '__main__':
    res = ResolweAPI(url=DEFAULT_URL)
    print(res)
