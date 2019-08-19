"""
Return types from api classes interface for the SOAP kegg api.

"""

from datetime import datetime
from operator import methodcaller
from collections import namedtuple

Definition = namedtuple("Definition", ["entry_id", "definition"])


def _Definition_from_items(items):
    """ Definition 'items' tuple from a list of items
    """
    items_list = []
    for name, val in items:
        if isinstance(name, list):
            name = name[0]
        if isinstance(val, list):
            val = val[0]
        items_list.append((str(name), str(val)))
    return Definition(**dict(items_list))


def _Definition_from_str(text):
    """
    Return a `Definition` item by parsing a tab separated string `text`
    i.e. text must be of the form '<entry_id>\t<definition>'

    """
    return Definition(*text.split("\t", 1))


Definition.from_items = staticmethod(_Definition_from_items)
Definition.from_str = staticmethod(_Definition_from_str)


OrganismSummary = namedtuple("OrganismSummary", ["entry_id", "org_code", "name", "lineage"])


def OrganismSummary_from_str(string):
    # string = string.decode("utf8")
    return OrganismSummary(*string.split("\t"))


OrganismSummary.from_str = staticmethod(OrganismSummary_from_str)


BInfo = namedtuple(
    'BInfo', ["entry_id", "definition", "name", "release", "curator", "contents", "last_update", "supported_formats"]
)


def _BInfo_from_text(text):
    """ Parse the return string from info into a new BInfo instance.
    """
    lines = text.splitlines()
    name, definition = lines[0].split(" ", 1)
    definition = definition.strip()
    entry_id, release = lines[1].split(" ", 1)
    _, release = release.strip().split(" ", 1)
    curator = lines[2].strip()
    contents = "\n".join(map(methodcaller("strip"), lines[3:]))

    return BInfo(entry_id, definition, name, release, curator, contents, None, None)


BInfo.from_text = staticmethod(_BInfo_from_text)


Link = namedtuple("Link", ["entry_id1", "entry_id2"])


SSDBRelation = namedtuple(
    "SSDBRelation",
    [
        "genes_id1",
        "genes_id2",
        "sw_score",
        "bit_score",
        "identity",
        "overlap",
        "start_position1",
        "end_position1",
        "start_position2",
        "end_position2",
        "best_flag_1to2",
        "best_flag_2to1",
        "definition1",
        "definition2",
        "length1",
        "length2",
    ],
)

MotifResult = namedtuple(
    "MotifResult", ["motif_id", "definition", "genes_id", "start_position", "end_position", "score", "evalue"]
)

LinkDBRelation = namedtuple("LinkDBRelation", ["entry_id1", "entry_id2", "type", "path"])

PathwayElement = namedtuple("PathwayElement", ["element_id", "type", "names", "components"])

PathwayElementRelation = namedtuple("PathwayElementRelation", ["element_id1", "element_id2", "type", "subtypes"])

Subtype = namedtuple("Subtype", ["element_id", "relation", "type"])
