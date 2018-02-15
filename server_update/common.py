"""common functions"""
import json
import os


from datetime import datetime


def file_size_bytes(file_path):
    """ returns file size in bytes"""
    return os.stat(file_path).st_size


def create_info_file(file_path, **kwargs):
    kwargs['datetime'] = '{0:%Y-%m-%d %H:%M:%S.%f}'.format(datetime.today())
    kwargs['size'] = file_size_bytes(file_path)
    print(file_path, file_size_bytes(file_path), os.stat(file_path))
    with open(file_path + '.info', 'wt') as f:
        json.dump(kwargs, f)


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
