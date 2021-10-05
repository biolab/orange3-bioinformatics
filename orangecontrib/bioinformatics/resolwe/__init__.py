""" Resolwe module """
from .utils import GENAPI_CACHE, CACHE_BACKEND, RESOLWEAPI_CACHE
from .genapi import GenAPI
from .resapi import ResolweAPI

RESOLWE_PLATFORM = 'resolwe'
GENESIS_PLATFORM = 'genesis'


class ResolweAuthError(Exception):
    """A login error occurred."""


class ResolweServerTypeError(Exception):
    """Unknown server type"""


def connect(username=None, password=None, url=None, server_type='resolwe'):
    """Connect to Resolwe server

    :param username:
    :type username: str

    :param password:
    :type password: str

    :param url:
    :type url: str

    :param server_type: genesis or resolwe
    :type server_type: str

    :return: Instance of GenAPI or ResolweAPI
    """

    if server_type == RESOLWE_PLATFORM:
        _api = ResolweAPI
    elif server_type == GENESIS_PLATFORM:
        _api = GenAPI
    else:
        raise ResolweServerTypeError(f'Unknown server type {server_type}')

    try:
        return _api(username, password, url)
    except Exception as e:
        raise ResolweAuthError(e.args[0]) from e
