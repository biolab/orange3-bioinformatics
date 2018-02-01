"""
obiKEGG2 configuration

mostly just caching settings

"""
from __future__ import absolute_import

import os

try:
    import ConfigParser as configparser
except ImportError:
    import configparser

from six import StringIO
from ..utils import serverfiles
kegg_dir = serverfiles.localpath("KEGG2")

default = """
[cache]
# path = %(home)s/.obiKEGG/
path = %(kegg_dir)s/
store = sqlite3
invalidate = weekly

[service]
transport = urllib2
# transport = requests

"""

# Orange kegg files dir

env = dict(os.environ)
env["kegg_dir"] = kegg_dir

parser = configparser.ConfigParser(env)


parser.readfp(StringIO(default), "default")

# TODO: global settings rc file
parser.read([os.path.expanduser("~/.obiKEGG/rc.cfg")])

params = {}

_ALL_PARAMS = [
    "cache.path",
    "cache.store",
    "cache.invalidate",
    "service.transport"
]

for p in _ALL_PARAMS:
    section, option = p.split(".")
    params[p] = parser.get(section, option)
