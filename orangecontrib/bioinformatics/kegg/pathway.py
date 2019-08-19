"""
============
KEGG Pathway
============

"""
from __future__ import absolute_import

import io
import os
import xml.parsers
from xml.dom import minidom
from functools import reduce
from contextlib import closing

import requests

from orangecontrib.bioinformatics.kegg import api, conf, caching


def cached_method(func, cache_name="_cached_method_cache", store=None):
    def wrapper(self, *args, **kwargs):
        sig = (func.__name__,) + args + tuple(sorted(kwargs.items()))
        if not hasattr(self, cache_name):
            setattr(self, cache_name, store() if store is not None else {})
        if sig not in getattr(self, cache_name):
            getattr(self, cache_name)[sig] = func(self, *args, **kwargs)
        return getattr(self, cache_name)[sig]

    return wrapper


class Pathway(object):
    """
    Class representing a KEGG Pathway (parsed from a "kgml" file)

    :param str pathway_id: A KEGG pathway id (e.g. 'path:hsa05130')

    """

    KGML_URL_FORMAT = "http://rest.kegg.jp/get/{pathway_id}/kgml"

    def __init__(self, pathway_id, local_cache=None, connection=None):
        if pathway_id.startswith("path:"):
            _, pathway_id = pathway_id.split(":", 1)

        self.pathway_id = pathway_id
        if local_cache is None:
            local_cache = conf.params["cache.path"]
        self.local_cache = local_cache
        self.connection = connection

    def cache_store(self):
        caching.touch_path(self.local_cache)
        return caching.Sqlite3Store(os.path.join(self.local_cache, "pathway_store.sqlite3"))

    def _open_last_modified_store(self):
        caching.touch_dir(self.local_cache)
        return caching.Sqlite3Store(os.path.join(self.local_cache, "last_modified.sqlite3"))

    def _get_kgml(self):
        """
        Return an open kgml file for the pathway.
        """
        kegg = api.CachedKeggApi()
        return io.BytesIO(kegg.get(self.pathway_id + "/kgml"))

    def _get_image_filename(self):
        """
        Return a filename of a local copy of the pathway image
        """
        url = str(self.image)

        local_filename = os.path.join(self.local_cache, self.pathway_id + ".png")

        if not os.path.exists(local_filename):
            r = requests.get(url, stream=True)
            modified_since = r.headers['last-modified']
            image = r.raw.read()
        else:
            return local_filename

            # TODO always return the item in cache
            # the whole expiration mechanism would have to be updated
            with closing(self._open_last_modified_store()) as store:
                modified_since = store.get(url, None)

            r = requests.get(url, headers=dict([("If-Modified-Since", modified_since)]), stream=True)
            if r.status_code == 304:
                return local_filename
            modified_since = r.headers["last-modified"]
            image = r.raw.read()

        with open(local_filename, "wb") as f:
            f.write(image)

        with closing(self._open_last_modified_store()) as store:
            store[url] = modified_since

        return local_filename

    def _local_kgml_filename(self):
        """
        Return the local kgml xml filename for the pathway.
        """
        local_filename = os.path.join(self.local_cache, self.pathway_id + ".xml")
        return local_filename

    class entry(object):
        def __init__(self, dom_element):
            self.__dict__.update(dom_element.attributes.items())
            self.graphics = ()
            self.components = []

            graphics = dom_element.getElementsByTagName("graphics")[0]
            self.graphics = dict(graphics.attributes.items())

            components = dom_element.getElementsByTagName("component")
            self.components = [node.getAttribute("id") for node in components]

    class reaction(object):
        def __init__(self, dom_element):
            self.__dict__.update(dom_element.attributes.items())
            self.substrates = [node.getAttribute("name") for node in dom_element.getElementsByTagName("substrate")]
            self.products = [node.getAttribute("name") for node in dom_element.getElementsByTagName("product")]

    class relation(object):
        def __init__(self, dom_element):
            self.__dict__.update(dom_element.attributes.items())
            self.subtypes = [node.attributes.items() for node in dom_element.getElementsByTagName("subtype")]

    @cached_method
    def pathway_attributes(self):
        if self.pathway_dom():
            return dict(self.pathway_dom().attributes.items())
        else:
            return None

    @property
    def name(self):
        """
        Pathway name/id (e.g. "path:hsa05130")
        """
        return self.pathway_attributes().get("name")

    @property
    def org(self):
        """
        Pathway organism code (e.g. 'hsa')
        """
        return self.pathway_attributes().get("org")

    @property
    def number(self):
        """
        Pathway number as a string (e.g. '05130')
        """
        return self.pathway_attributes().get("number")

    @property
    def title(self):
        """
        Pathway title string.
        """
        return self.pathway_attributes().get("title")

    @property
    def image(self):
        """
        URL of the pathway image.
        """
        return self.pathway_attributes().get("image")

    @property
    def link(self):
        """
        URL to a pathway on the KEGG web site.
        """
        return self.pathway_attributes().get("link")

    @cached_method
    def pathway_dom(self):
        with self._get_kgml() as kgml:
            try:
                return minidom.parse(kgml).getElementsByTagName("pathway")[0]
            except xml.parsers.expat.ExpatError:
                # TODO: Should delete the cached xml file.
                return None

    @cached_method
    def entries(self):
        dom = self.pathway_dom()
        if dom:
            return [self.entry(e) for e in dom.getElementsByTagName("entry")]
        else:
            return []

    @cached_method
    def reactions(self):
        dom = self.pathway_dom()
        if dom:
            return [self.reaction(e) for e in dom.getElementsByTagName("reaction")]
        else:
            return []

    @cached_method
    def relations(self):
        dom = self.pathway_dom()
        if dom:
            return [self.relation(e) for e in dom.getElementsByTagName("relation")]
        else:
            return []

    def __iter__(self):
        """
        Iterate over all elements in the pathway.
        """
        return iter(self.all_elements())

    def __contains__(self, element):
        """
        Return ``True`` if element in the pathway.
        """
        return element in self.all_elements()

    @classmethod
    def split_pathway_id(cls, id):
        path, id = id.split(":") if ":" in id else ("path", id)
        org, id = id[:-5], id[-5:]
        return path, org, id

    @cached_method
    def all_elements(self):
        """
        Return all elements
        """
        return reduce(list.__add__, [self.genes(), self.compounds(), self.enzmes(), self.reactions()], [])

    def _get_entries_by_type(self, type1):
        return sorted(
            reduce(set.union, [entry.name.split() for entry in self.entries() if entry.type == type1], set())
        )

    @cached_method
    def genes(self):
        """
        Return all genes on the pathway.
        """
        return self._get_entries_by_type("gene")

    @cached_method
    def compounds(self):
        """
        Return all compounds on the pathway.
        """
        return self._get_entries_by_type("compound")

    @cached_method
    def enzymes(self):
        """
        Return all enzymes on the pathway.
        """
        return self._get_entries_by_type("enzyme")

    @cached_method
    def orthologs(self):
        """
        Return all orthologs on the pathway.
        """
        return self._get_entries_by_type("ortholog")

    @cached_method
    def maps(self):
        """
        Return all linked maps on the pathway.
        """
        return self._get_entries_by_type("map")

    @cached_method
    def groups(self):
        """
        Return all groups on the pathway.
        """
        return self._get_entries_by_type("ortholog")

    def get_image(self):
        """
        Return an local filesystem path to an image of the pathway. The image
        will be downloaded if not already cached.
        """
        return self._get_image_filename()

    @classmethod
    def list(cls, organism):
        """
        List all pathways for KEGG organism code `organism`.
        """
        kegg = api.CachedKeggApi()
        return kegg.list_pathways(organism)
