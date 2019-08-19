"""
KEGG Brite

"""
from __future__ import absolute_import

import io
import os
import re

from orangecontrib.bioinformatics.kegg import conf

try:
    from urllib2 import urlopen
except ImportError:
    from urllib.request import urlopen


class BriteEntry(object):
    _search_re = {
        "ids": re.compile('(?P<ids>\[.*:.*\])'),
        "title": re.compile(r'(<[Bb]>)?(?P<title>\b[a-zA-Z0-9_/\s,;:.+=\-\[\]{}\(\)]+?)(?(1)</[Bb]>)$'),
        "links": re.compile('(?P<links><a href=".+?">.*?</a>)'),
    }

    def __init__(self, line):
        self.entries = []
        self.line = line[1:].strip()
        for name, re in self._search_re.items():
            search = re.search(self.line)
            setattr(self, name, search.group(name) if search else None)

    def __iter__(self):
        return iter(self.entries)


class Brite(BriteEntry):
    VERSION = "v1.0"
    BRITE_URL_FORMAT = "http://www.genome.jp/kegg-bin/download_htext?htext={brite_id}.keg&format=htext&filedir="

    def __init__(self, brite_id, local_cache=None):
        super(Brite, self).__init__("")
        self.brite_id = id
        if local_cache is None:
            local_cache = conf.params["cache.path"]
        self.local_cache = local_cache

        self.load(brite_id)

    def _get_brite(self, brite_id):
        url = self.BRITE_URL_FORMAT.format(brite_id=brite_id)
        local_filename = os.path.join(self.local_cache, brite_id + ".keg")
        if not os.path.exists(local_filename):
            brite = urlopen(url).read()
            with io.open(local_filename, "wb") as f:
                f.write(brite)

        return io.open(local_filename, "r")

    def load(self, brite_id):
        lines = self._get_brite(brite_id).read().split("\n!\n")[1].splitlines()

        # TODO: Implement a proper parser

        def collect(lines, depth, collection):
            while lines:
                line = lines[0]
                if line.startswith("#"):
                    lines.pop(0)
                elif line.startswith(depth) and len(line.strip()) > 1:
                    collection.append(BriteEntry(lines.pop(0)))
                elif line[0] > depth:
                    collect(lines, line[0], collection[-1].entries)
                elif line[0] < depth:
                    return
                else:
                    lines.pop(0)

        collect([line for line in lines if not line.startswith("#") and len(line) > 1], "A", self.entries)
