import os

from . import environ


def progress_bar_milestones(count, iterations=100):
    return set([int(i*count/float(iterations)) for i in range(iterations)])


buffer_folder = os.path.join(environ.buffer_dir, 'serverfiles-bio/')


def ensure_type(value, types):
    if isinstance(value, types):
        return value
    else:
        raise TypeError('Wrong variable type. {value} is {value_type}, but should be {types}'.format(
            value=value, value_type=type(value), types=types))
