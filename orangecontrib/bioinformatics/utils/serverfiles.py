"""ServerFiles"""
import serverfiles


from orangecontrib.bioinformatics.utils import serverfile_path

server_url = "http://orange.biolab.si/serverfiles-bio2/"


class ServerFiles(serverfiles.ServerFiles):

    def __init__(self, server=server_url):
        serverfiles.ServerFiles.__init__(self, server)


PATH = serverfile_path()
LOCALFILES = serverfiles.LocalFiles(PATH, serverfiles=ServerFiles())


def localpath(*args, **kwargs):
    return LOCALFILES.localpath(*args, **kwargs)

    
def listfiles(*args, **kwargs):
    return [fname for domain, fname in LOCALFILES.listfiles(*args, **kwargs)]


def localpath_download(*path, **kwargs):
    files = serverfiles.LocalFiles(PATH, serverfiles=ServerFiles())
    return files.localpath_download(*path, **kwargs)


def download(*args, **kwargs):
    files = serverfiles.LocalFiles(PATH, serverfiles=ServerFiles())
    return files.download(*args, **kwargs)


def allinfo(*args, **kwargs):
    files = serverfiles.LocalFiles(PATH, serverfiles=ServerFiles())
    return files.allinfo(*args, **kwargs)


def info(*args, **kwargs):
    files = serverfiles.LocalFiles(PATH, serverfiles=ServerFiles())
    return files.info(*args, **kwargs)


def update(*path, **kwargs):
    files = serverfiles.LocalFiles(PATH, serverfiles=ServerFiles())
    return files.update(*path, **kwargs)


def sizeformat(size):
    return serverfiles.sizeformat(size)
