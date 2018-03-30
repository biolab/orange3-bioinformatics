"""common functions"""
import json
import os

from collections import OrderedDict
from datetime import datetime

INFO_FILE_SCHEMA = {
    'domain': None,
    'filename': None,
    'source': None,
    'title': None,
    'tags': [],
    'size': None,
    'datetime': None
    # used only if files are compressed
    # 'uncompressed': None,
    # 'compression': None,
}


def file_size_bytes(file_path):
    """ returns file size in bytes """
    return os.stat(file_path).st_size


def create_info_file(file_path, **kwargs):
    info_dict = OrderedDict(INFO_FILE_SCHEMA)

    info_dict.update(**kwargs)
    info_dict['datetime'] = '{0:%Y-%m-%d %H:%M:%S.%f}'.format(datetime.today())
    info_dict['size'] = file_size_bytes(file_path)

    with open(file_path + '.info', 'wt') as f:
        json.dump(info_dict, f)


def create_folder(path):
    try:
        os.makedirs(path)
    except OSError:
        if os.path.exists(path):
            #raise
            pass
        else:
            # There was an error on creation, so make sure we know about it
            raise
