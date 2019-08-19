"""ServerFiles"""
import os

import serverfiles

from orangecontrib.bioinformatics.utils import local_cache

version = 'v1'
server_url = f'https://download.biolab.si/datasets/bioinformatics/{version}/'


class ServerFiles(serverfiles.ServerFiles):
    def __init__(self, server=server_url):
        serverfiles.ServerFiles.__init__(self, server)


PATH = os.path.join(local_cache, 'serverfiles', version)
LOCALFILES = serverfiles.LocalFiles(PATH, serverfiles=ServerFiles())


def localpath(*args, **kwargs):
    return LOCALFILES.localpath(*args, **kwargs)


def listfiles(*args):
    files = serverfiles.LocalFiles(PATH, serverfiles=ServerFiles())
    return files.listfiles(*args)


def localpath_download(*path, **kwargs):
    files = serverfiles.LocalFiles(PATH, serverfiles=ServerFiles())
    return files.localpath_download(*path, **kwargs)


def download(*path, **kwargs):
    files = serverfiles.LocalFiles(PATH, serverfiles=ServerFiles())
    return files.download(*path, **kwargs)


def allinfo(*args, **kwargs):
    files = serverfiles.LocalFiles(PATH, serverfiles=ServerFiles())
    return files.allinfo(*args, **kwargs)


def info(*args, **kwargs):
    files = serverfiles.LocalFiles(PATH, serverfiles=ServerFiles())
    return files.info(*args, **kwargs)


def need_update(*path):
    files = serverfiles.LocalFiles(PATH, serverfiles=ServerFiles())
    return files.needs_update(*path)


def update(*path, **kwargs):
    files = serverfiles.LocalFiles(PATH, serverfiles=ServerFiles())
    return files.update(*path, **kwargs)


def sizeformat(size):
    return serverfiles.sizeformat(size)
