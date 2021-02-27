""" GenAPI """
import re

import requests
import requests_cache
from genesis import GenData, Genesis

from orangecontrib.bioinformatics.resolwe.utils import GENAPI_CACHE, CACHE_BACKEND, response_to_json

view_model = ['experiment', 'growth', 'genotype', 'treatment', 'strain', 'time', 'replicate']
DEFAULT_EMAIL = 'anonymous@genialis.com'
DEFAULT_PASSWD = 'anonymous'
DEFAULT_URL = 'https://dictyexpress.research.bcm.edu'


class GenAPI:
    """Python module that leverages Genesis PyAPI (Python API for accsess to DictyExpress database).

    It supports connection to the server and data retrieval functionalities.
    """

    def __init__(self, email=DEFAULT_EMAIL, password=DEFAULT_PASSWD, url=DEFAULT_URL):
        self._gen = Genesis(email, password, url)
        self.email = email
        self._data_endpoint = url + '/data/'

    @staticmethod
    def clear_cache():
        with requests_cache.enabled(cache_name=GENAPI_CACHE, backend=CACHE_BACKEND):
            requests_cache.clear()

    def get_cached_ids(self):
        with requests_cache.enabled(cache_name=GENAPI_CACHE, backend=CACHE_BACKEND):
            cached_object = requests_cache.core.get_cache()
            responses = [cached_object.get_response_and_time(response) for response in cached_object.responses]
            gen_ids = []

            for url in [response.url for response, _ in responses]:
                gen_id = re.search(r'{}(.*?)/'.format(self._data_endpoint), url)

                if gen_id is not None:
                    gen_ids.append(gen_id.group(1))

            return gen_ids

    def fetch_etc_objects(self, *args, **kwargs):
        """Function downloads all available :obj:`GenData` etc objects from DictyExpress database.

        :rtype: list of GenData objects
        """

        with requests_cache.enabled(cache_name=GENAPI_CACHE, backend=CACHE_BACKEND):
            try:
                # Note: this is hardcoded for now. When we port this module to Resolwe platform
                #       data retrieval will be different
                list_of_experiments = self._gen.api.data.get(
                    case_ids__contains='5535115cfad58d5e03006217', status='done', type__startswith='data:etc:'
                )['objects']

            except requests.exceptions.ConnectionError as e:
                raise requests.exceptions.ConnectionError('Server not accessible, check your connection.') from e

            store_experiments = [GenData(exp, self._gen) for exp in list_of_experiments]

            return store_experiments

    def download_etc_data(self, gen_data_id, *args, **kwargs):
        """Function downloads etc data of a chosen experiment from the server.

        :param gen_data_id: id of GeneData object
        :type gen_data_id: str

        :rtype: data in json like format
        """
        table_name = kwargs.get("table_name", '')

        with requests_cache.enabled(cache_name=GENAPI_CACHE, backend=CACHE_BACKEND):
            try:
                response = next(self._gen.download([gen_data_id], 'output.etcfile'))

            except requests.exceptions.ConnectionError as e:
                raise requests.exceptions.ConnectionError('Server not accessible, check your connection.') from e

            return response_to_json(response), table_name
