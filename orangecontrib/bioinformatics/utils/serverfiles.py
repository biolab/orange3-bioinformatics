"""ServerFiles"""
import serverfiles


from orangecontrib.bioinformatics.utils import buffer_folder

server_url = "http://orange.biolab.si/serverfiles-bio2/"


class ServerFiles(serverfiles.ServerFiles):

    def __init__(self, server=server_url):
        serverfiles.ServerFiles.__init__(self, server)


PATH = buffer_folder
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


def update(*path, **kwargs):
    files = serverfiles.LocalFiles(PATH, serverfiles=ServerFiles())
    return files.update(*path, **kwargs)


def sizeformat(size):
    return serverfiles.sizeformat(size)
