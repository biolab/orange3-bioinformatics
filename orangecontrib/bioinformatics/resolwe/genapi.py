""" GenAPI """
from collections import defaultdict
import re
from dataclasses import dataclass

import requests
import requests_cache

from orangecontrib.bioinformatics.resolwe.utils import GENAPI_CACHE, CACHE_BACKEND

view_model = ['experiment', 'growth', 'genotype', 'treatment', 'strain', 'time', 'replicate']
DEFAULT_EMAIL = 'anonymous@genialis.com'
DEFAULT_PASSWD = 'anonymous'
DEFAULT_URL = 'https://app.dictyexpress.org'
API_URL = DEFAULT_URL + '/api'
CREDENTIAL_MANAGER_SERVICE: str = 'genapi_credentials'


@dataclass
class GenData:
    """Compatibility wrapper for dictyExpress time series relations."""

    id: int
    name: str
    type: str
    annotation: dict
    var: dict
    static: dict
    relation: dict


class GenAPI:
    """Python module that leverages Genesis PyAPI (Python API for accsess to DictyExpress database).

    It supports connection to the server and data retrieval functionalities.
    """

    def __init__(self, email=DEFAULT_EMAIL, password=DEFAULT_PASSWD, url=DEFAULT_URL):
        self.email = email
        self.url = url.rstrip('/')
        self.api_url = self.url + '/api'
        self._data_endpoint = url + '/data/'
        self._downloaded_ids = set()
        self._relations = {}
        self._session = requests.Session()
        self._session.headers.update(
            {'Accept': 'application/json', 'Referer': self.url + '/bcm/'}
        )

    @staticmethod
    def clear_cache():
        with requests_cache.enabled(cache_name=GENAPI_CACHE, backend=CACHE_BACKEND):
            requests_cache.clear()

    def get_cached_ids(self):
        return list(self._downloaded_ids)

    def fetch_etc_objects(self, *args, **kwargs):
        """Function downloads all available :obj:`GenData` etc objects from DictyExpress database.

        :rtype: list of GenData objects
        """

        with requests_cache.enabled(cache_name=GENAPI_CACHE, backend=CACHE_BACKEND):
            try:
                response = self._session.get(
                    f'{self.api_url}/relation',
                    params={'category': 'Time series', 'tags': 'community:bcm'},
                )
                response.raise_for_status()
                list_of_experiments = response.json()

            except requests.exceptions.ConnectionError as e:
                raise requests.exceptions.ConnectionError('Server not accessible, check your connection.') from e

            store_experiments = [self._relation_to_gen_data(exp) for exp in list_of_experiments]
            self._relations = {exp.id: exp.relation for exp in store_experiments}

            return store_experiments

    def download_etc_data(self, gen_data_id, state=None, *args, **kwargs):
        """Function downloads etc data of a chosen experiment from the server.

        :param gen_data_id: id of GeneData object
        :type gen_data_id: str

        :rtype: data in json like format
        """
        table_name = kwargs.get("table_name", '')

        with requests_cache.enabled(cache_name=GENAPI_CACHE, backend=CACHE_BACKEND):
            try:
                self._set_status(state, 'Downloading time series ...')
                relation = self._relations.get(gen_data_id) or self._fetch_relation(gen_data_id)
                etc_json = self._download_relation_data(relation, state)

            except requests.exceptions.ConnectionError as e:
                raise requests.exceptions.ConnectionError('Server not accessible, check your connection.') from e

            self._downloaded_ids.add(gen_data_id)
            return etc_json, table_name

    def _fetch_relation(self, gen_data_id):
        response = self._session.get(f'{self.api_url}/relation', params={'id': gen_data_id})
        response.raise_for_status()
        relations = response.json()

        if not relations:
            raise ValueError(f'Time series {gen_data_id} does not exist.')

        return relations[0]

    def _relation_to_gen_data(self, relation):
        descriptor = relation.get('descriptor', {})
        citation = descriptor.get('citation') or {}
        citation_name = citation.get('name', '') if isinstance(citation, dict) else citation
        collection = relation.get('collection') or {}
        experiment = collection.get('name', '')

        annotation = {
            'var.project': {'value': descriptor.get('project', '')},
            'static.name': {'value': experiment},
            'static.cite': {'value': [{'name': citation_name}] if citation_name else ''},
            'var.growth': {'value': descriptor.get('growth', '')},
            'var.treatment': {'value': descriptor.get('treatment', '')},
            'var.strain': {'value': descriptor.get('strain', '')},
        }

        return GenData(
            id=relation['id'],
            name=descriptor.get('details') or experiment,
            type='data:etc:',
            annotation=annotation,
            var={'project': descriptor.get('project', '')},
            static={'name': experiment},
            relation=relation,
        )

    def _download_relation_data(self, relation, state=None):
        partitions = sorted(
            relation.get('partitions', []),
            key=lambda part: (part.get('position', 0), part.get('id', 0)),
        )
        entities = [str(part['entity']) for part in partitions]

        if not entities:
            raise ValueError('Selected time series does not contain samples.')

        self._set_progress(state, 5)
        expressions = self._fetch_expressions(entities)
        expression_by_entity = {
            expression['entity']['id']: expression
            for expression in expressions
            if expression.get('entity') and expression.get('output', {}).get('exp_json')
        }
        self._set_progress(state, 10)

        genes_by_time = defaultdict(lambda: defaultdict(list))
        time_points = []
        total_partitions = len(partitions)

        for index, partition in enumerate(partitions, start=1):
            self._check_interruption(state)
            expression = expression_by_entity.get(partition['entity'])
            if expression is None:
                continue

            time_point = self._time_point(partition)
            if time_point not in time_points:
                time_points.append(time_point)

            storage_id = expression['output']['exp_json']
            self._set_status(state, f'Downloading sample {index} of {total_partitions} ...')
            storage = self._session.get(f'{self.api_url}/storage/{storage_id}')
            storage.raise_for_status()

            for gene, value in storage.json()['json']['genes'].items():
                genes_by_time[gene][time_point].append(float(value))

            self._set_progress(state, 10 + 80 * index / total_partitions)

        self._set_status(state, 'Combining replicates ...')
        genes = {
            gene: [
                sum(values_by_time.get(time_point, [])) / len(values_by_time[time_point])
                if values_by_time.get(time_point)
                else None
                for time_point in time_points
            ]
            for gene, values_by_time in genes_by_time.items()
        }

        self._set_progress(state, 100)
        return {'etc': {'genes': genes, 'timePoints': time_points}}

    def _fetch_expressions(self, entities):
        response = self._session.get(
            f'{self.api_url}/data',
            params={'type': 'data:expression', 'entity__in': ','.join(entities)},
        )
        response.raise_for_status()
        data = response.json()
        return data.get('results', data) if isinstance(data, dict) else data

    @staticmethod
    def _time_point(partition):
        label = partition.get('label', '')
        match = re.search(r'\d+', label)
        return int(match.group(0)) if match else partition.get('position', 0)

    @staticmethod
    def _set_status(state, status):
        if state is not None:
            state.set_status(status)

    @staticmethod
    def _set_progress(state, progress):
        if state is not None:
            state.set_progress_value(progress)

    @staticmethod
    def _check_interruption(state):
        if state is not None and state.is_interruption_requested():
            raise InterruptedError
