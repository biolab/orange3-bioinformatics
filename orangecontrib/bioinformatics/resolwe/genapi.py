""" GenAPI """
import json
import io
import gzip
import requests


from genesis import Genesis, GenData


view_model = ['experiment', 'growth', 'genotype', 'treatment', 'strain', 'time', 'replicate']
DEFAULT_EMAIL = 'anonymous@genialis.com'
DEFAULT_PASSWD = 'anonymous'
DEFAULT_URL = 'https://dictyexpress.research.bcm.edu'


class GenAPI:

    def __init__(self, email=DEFAULT_EMAIL, password=DEFAULT_PASSWD, url=DEFAULT_URL):
        """ Python module that leverages Genesis PyAPI (Python API for accsess to DictyExpress database).

        It supports connection to the server and data retrieval functionalities.

        Args:
            email (str):
            password (str):
            url (str):
        """

        self._gen = Genesis(email, password, url)
        self.email = email

    def fetch_etc_objects(self, **kwargs):
        """ Function downloads all available :obj:`GenData` etc objects from DictyExpress database.

        Returns:
            :obj:`list`: :obj:`GenData` objects

        """

        callback = kwargs.get("progress_callback", None)

        try:
            # Note: this is hardcoded for now. When we port this module to Resolwe platform
            #       data retrieval will be different
            list_of_experiments = self._gen.api.data.get(case_ids__contains='5535115cfad58d5e03006217', status='done',
                                                         type__startswith='data:etc:')['objects']
            if callback:
                callback.emit()

        except requests.exceptions.ConnectionError as e:
            raise requests.exceptions.ConnectionError('Server not accessible, check your connection.') from e

        store_experiments = [GenData(exp, self._gen) for exp in list_of_experiments]

        if callback:
            callback.emit()

        return store_experiments

    def download_etc_data(self, gen_data_id, **kwargs):
        """ Function downloads etc data of a chosen experiment from the server.

        Args:
            gen_data_id (str): id of :obj:`GenData` object to download.

        Returns:
             :obj:`dict`: data in json like format


        """
        callback = kwargs.get("progress_callback", None)

        try:
            response = next(self._gen.download([gen_data_id], 'output.etcfile'))
            # TODO: maybe edit Genesis module to support callback?
            if callback:
                callback.emit()

        except requests.exceptions.ConnectionError as e:
            raise requests.exceptions.ConnectionError('Server not accessible, check your connection.') from e

        if not response.ok:
            response.raise_for_status()

        response_gzipped = io.BytesIO(response.content)

        response_content = io.TextIOWrapper(gzip.GzipFile(fileobj=response_gzipped), encoding='utf-8')

        if callback:
            callback.emit()

        try:
            json_object = json.load(response_content)
        except ValueError as e:
            raise ValueError('Downloaded data is not a valid JSON') from e

        return json_object
