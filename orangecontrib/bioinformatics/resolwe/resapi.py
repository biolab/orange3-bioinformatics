""" ResolweAPI """
from typing import Dict, List
from pathlib import Path
from tempfile import NamedTemporaryFile
from urllib.parse import urljoin

import requests
from resdk import Resolwe
from requests import Response
from resdk.resources.data import Data, Sample

from Orange.data import Table, Variable

from orangecontrib.bioinformatics.resolwe import sample_metadata

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

    def download(self, data_id, data_file_name) -> Response:
        file_url = urljoin(self.url, f'data/{data_id}/{data_file_name}')
        response = self.session.get(file_url, auth=self._res.auth)
        response.raise_for_status()
        return response

    def get_expression_data_objects(self, collection_id: str) -> List[Data]:
        """ Get all expression data objects from a collection. """
        return self._res.data.filter(type='data:expression', collection=collection_id)

    def get_orange_data_objects(self, collection_id: str) -> List[Data]:
        return self._res.data.filter(type='data:metadata:orange', collection=collection_id, ordering='-modified')

    def get_sample_relations(self, rel_type: str, collection_id: str):
        """ Get all relations (by type) from a collection. """
        return self._res.relation.filter(collection=collection_id).filter(type=rel_type)

    def get_samples_by_id(self, sample_ids: List[str]) -> List[Sample]:
        return self._res.sample.filter(id__in=sample_ids)

    def aggregate_sample_metadata(self, samples: List[Sample], collection_id: str) -> Dict[Variable, List[List[str]]]:
        replicate_groups = list(self.get_sample_relations('group', collection_id))
        replicate_series = list(self.get_sample_relations('series', collection_id))

        descriptors = sample_metadata.descriptors(samples)
        relations = sample_metadata.relations(samples, replicate_groups + replicate_series)

        metadata = {}
        metadata.update(descriptors)
        metadata.update(relations)

        return metadata

    def get_clinical_metadata(self, data_id, data_file_name) -> Table:
        response = self.download(data_id, data_file_name)
        with NamedTemporaryFile(suffix=''.join(Path(data_file_name).suffixes)) as f:
            f.write(response.content)
            f.seek(0)
            return Table.from_file(f.name)


if __name__ == '__main__':
    res = ResolweAPI(url=DEFAULT_URL)
    print(res)
