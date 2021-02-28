""" Resolwe module """
from .utils import GENAPI_CACHE, CACHE_BACKEND, RESOLWEAPI_CACHE
from .genapi import DEFAULT_EMAIL, DEFAULT_PASSWD, GenAPI
from .resapi import ResolweAPI

__all__ = ('DEFAULT_EMAIL', 'DEFAULT_PASSWD')


class ResolweAuthException(Exception):
    """ A login error occurred. """


class ResolweServerTypeException(Exception):
    """ Unknown server type """


class ResolweDataObjectsNotFound(Exception):
    """ Data Objects not found """


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

    if server_type == 'resolwe':
        _api = ResolweAPI
    elif server_type == 'genesis':
        _api = GenAPI
    else:
        raise ResolweServerTypeException(f'Unknown server type {server_type}')

    try:
        return _api(username, password, url)
    except Exception as e:
        raise ResolweAuthException(e.args[0]) from e
