import os

from Orange.misc.environ import data_dir


def progress_bar_milestones(count, iterations=100):
    return {int(i * count / float(iterations)) for i in range(iterations)}


local_cache = os.path.join(data_dir(), 'bioinformatics/')


def ensure_type(value, types):
    if isinstance(value, types):
        return value
    else:
        raise TypeError(
            'Wrong variable type. {value} is {value_type}, but should be {types}'.format(
                value=value, value_type=type(value), types=types
            )
        )
