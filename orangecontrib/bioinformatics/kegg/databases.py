"""
DBGET Database Interface
========================

"""
from __future__ import absolute_import

import re
import sys
from contextlib import closing

from orangecontrib.bioinformatics.kegg import api, entry
from orangecontrib.bioinformatics.kegg.entry import fields


def iter_take(source_iter, n):
    """
    Return a list of the first `n` items in `source_iter`.
    """
    source_iter = iter(source_iter)
    return [item for _, item in zip(range(n), source_iter)]


def batch_iter(source_iter, n):
    """
    Split the `source_iter` into batches of size `n`.
    """
    source_iter = iter(source_iter)
    while True:
        batch = iter_take(source_iter, n)
        if batch:
            yield batch
        else:
            break


def chain_iter(chains_iter):
    for iter in chains_iter:
        for element in iter:
            yield element


# TODO: DBDataBase should be able to be constructed from a flat text
# entry file. The precache etc. should be moved in caching api, that creates
# simple file system hierarchy where the flat database is saved (with db
# release string), e.g.
# genes/hsa.dbget
# genes/hsa.release
# genes/sce.dbget
# path.dbget
# module.dbget
# ligand/compound.dbget


class DBDataBase(object):
    """
    Base class for a DBGET database interface.

    """

    #: ENTRY_TYPE constructor (a :class:`~.entry.DBEntry` subclass). This
    #: should be redefined in subclasses.
    ENTRY_TYPE = entry.DBEntry

    #: A database name/abbreviation (e.g. 'pathway'). Needs to be set in a
    #: subclass or object instance's constructor before calling the base.
    #: __init__
    DB = None

    def __init__(self, **kwargs):
        if not self.DB:
            raise TypeError("Cannot make an instance of abstract base " "class %r." % type(self).__name__)

        self.api = api.CachedKeggApi()
        self._info = None
        # TODO invalidate cache by KEGG release
        # self.api.set_default_release(self.info.release)
        self._keys = []

    @property
    def info(self):  # lazy info loading
        if not self._info:
            self._info = self.api.info(self.DB)
        return self._info

    def iterkeys(self):
        """
        Return an iterator over the `keys`.
        """
        return iter(self._keys)

    def iteritems(self):
        """
        Return an iterator over the `items`.
        """
        batch_size = 100
        iterkeys = self.iterkeys()
        return chain_iter(zip(batch, self.batch_get(batch)) for batch in batch_iter(iterkeys, batch_size))

    def itervalues(self):
        """
        Return an iterator over all :obj:`DBDataBase.ENTRY_TYPE` instances.
        """
        batch_size = 100
        iterkeys = self.iterkeys()
        return chain_iter(self.batch_get(batch) for batch in batch_iter(iterkeys, batch_size))

    if sys.version_info < (3,):

        def keys(self):
            """
            Return a list of database keys. These are unique KEGG identifiers
            that can be used to query the database.
            """
            return list(self._keys)

        def values(self):
            """
            Return a list of all :obj:`DBDataBase.ENTRY_TYPE` instances.
            """
            return self.batch_get(self.keys())

        def items(self):
            """
            Return a list of all (key, :obj:`DBDataBase.ENTRY_TYPE` instance)
            tuples.
            """
            return list(zip(self.keys(), self.batch_get(self.keys())))

    else:

        def keys(self):
            """
            Return an iterator over all database keys. These are unique
            KEGG identifiers that can be used to query the database.
            """
            return iter(self._keys)

        def values(self):
            """
            Return an iterator over all :obj:`DBDataBase.ENTRY_TYPE` instances.
            """
            return self.itervalues()

        def items(self):
            """
            Return an iterator over all (key, :obj:`DBDataBase.ENTRY_TYPE`)
            tuples.
            """
            return self.iteritems()

    def get(self, key, default=None):
        """
        Return an :obj:`DBDataBase.ENTRY_TYPE` instance for the `key`.
        Raises :class:`KeyError` if not found.

        """
        try:
            return self.__getitem__(key)
        except KeyError:
            return default

    def has_key(self, key):
        return self.__contains__(key)

    def __getitem__(self, key):
        e = self.get_entry(key)
        if e is None:
            raise KeyError(key)
        else:
            return e

    def __contains__(self, key):
        return key in set(self.keys())

    def __len__(self):
        return len(self.keys())

    def __iter__(self):
        return iter(self.keys())

    def get_text(self, key):
        """
        Return the database entry for `key` as plain text.
        """
        key = self._add_db(key)
        return self.api.get([key])

    def get_entry(self, key):
        """
        Return the database entry for `key` as an instance of `ENTRY_TYPE`.
        """
        text = self.get_text(key)
        if not text or text == "None":
            return None
        else:
            return self.ENTRY_TYPE(text)

    def find(self, name):
        """
        Find `name` using kegg `find` api.
        """
        res = self.api.find(self.DB, name).splitlines()
        return [r.split(" ", 1)[0] for r in res]

    def pre_cache(self, keys=None, batch_size=10, progress_callback=None):
        """
        Retrieve all the entries for `keys` and cache them locally for faster
        subsequent retrieval. If `keys` is ``None`` then all entries will be
        retrieved.

        """
        if not isinstance(self.api, api.CachedKeggApi):
            raise TypeError("Not an instance of api.CachedKeggApi")

        if batch_size > 10 or batch_size < 1:
            raise ValueError("Invalid batch_size")

        if keys is None:
            keys = self.keys()

        keys = map(self._add_db, keys)

        get = self.api.get

        # drop all keys with a valid cache entry to minimize the number
        # of 'get' requests.
        with closing(get.cache_store()) as store:

            def is_uncached(key):
                cache_key = get.key_from_args((key,))
                return not get.key_has_valid_cache(cache_key, store)

            keys = [key for key in keys if is_uncached(key)]

        start = 0

        while start < len(keys):
            batch = keys[start : start + batch_size]
            self.api.get(batch)

            if progress_callback:
                progress_callback(100.0 * start / len(keys))

            start += batch_size

    def batch_get(self, keys):
        """
        Batch retrieve all entries for keys. This can be significantly
        faster then getting each entry separately especially if entries
        are not yet cached.

        """
        entries = []
        batch_size = 10
        keys = list(map(self._add_db, keys))

        # Precache the entries first
        self.pre_cache(keys)

        start = 0
        while start < len(keys):
            batch = keys[start : start + batch_size]
            batch_entries = self.api.get(batch)
            if batch_entries is not None:
                batch_entries = batch_entries.split("///\n")
                # Remove possible empty last line
                batch_entries = [e for e in batch_entries if e.strip()]
                entries.extend(map(self.ENTRY_TYPE, batch_entries))
            start += batch_size

        return entries

    def _add_db(self, key):
        """
        Prefix the key with '%(DB)s:' string if not already prefixed.
        """
        if not key.startswith(self.DB + ":"):
            return self.DB + ":" + key
        else:
            return key


@entry.entry_decorate
class GenomeEntry(entry.DBEntry):
    """
    Entry for a KEGG Genome database.
    """

    FIELDS = [
        ("ENTRY", fields.DBEntryField),
        ("NAME", fields.DBNameField),
        ("DEFINITION", fields.DBDefinitionField),
        ("ANNOTATION", fields.DBSimpleField),
        ("TAXONOMY", fields.DBTaxonomyField),
        ("DATA_SOURCE", fields.DBSimpleField),
        ("ORIGINAL_DB", fields.DBSimpleField),
        ("KEYWORDS", fields.DBSimpleField),
        ("DISEASE", fields.DBSimpleField),
        ("COMMENT", fields.DBSimpleField),
        ("CHROMOSOME", fields.DBFieldWithSubsections),
        ("PLASMID", fields.DBSimpleField),
        ("STATISTICS", fields.DBSimpleField),
        ("REFERENCE", fields.DBReference),
    ]

    MULTIPLE_FIELDS = ["REFERENCE"]

    def __init__(self, text):
        entry.DBEntry.__init__(self, text)

    @property
    def organism_code(self):
        """
        A three or four letter KEGG organism code (e.g. 'hsa', 'sce', ...)
        """
        return self.name.split(",", 1)[0]

    @property
    def taxid(self):
        """
        Organism NCBI taxonomy id.
        """
        return self.TAXONOMY.taxid

    def org_code(self):
        # for backwards compatibility; return the `organism_code`
        return self.organism_code


class Genome(DBDataBase):
    """
    An interface to the A KEGG GENOME database.
    """

    DB = "genome"
    ENTRY_TYPE = GenomeEntry

    # For obiTaxonomy.common_taxids mapping
    TAXID_MAP = {
        "562": "511145",  # Escherichia coli K-12 MG1655
        "2104": "272634",  # Mycoplasma pneumoniae M129
        "4530": "39947",  # Oryza sativa ssp. japonica cultivar Nipponbare (Japanese rice)
        "4932": "559292",  # Saccharomyces cerevisiae S288C
        "4896": "284812",  # Schizosaccharomyces pombe 972h-
    }

    def __init__(self):
        DBDataBase.__init__(self)
        self._org_list = self.api.list_organisms()
        self._keys = [org.entry_id for org in self._org_list]

    def _key_to_gn_entry_id(self, key):
        res = self.find(key)
        if len(res) == 0:
            raise KeyError("Unknown key")
        elif len(res) > 1:
            raise ValueError("Not a unique key")
        else:
            return res[0]

    @classmethod
    def common_organisms(cls):
        return [
            'ath',
            'bta',
            'cel',
            'cre',
            'dre',
            'ddi',
            'dme',
            'eco',
            'hsa',
            'mmu',
            'mpn',
            'osa',
            'pfa',
            'rno',
            'sce',
            'spo',
            'zma',
            'xla',
        ]

    @classmethod
    def essential_organisms(cls):
        return ['ddi', 'dme', 'hsa', 'mmu', 'sce']

    def org_code_to_entry_key(self, code):
        """
        Map an organism code ('hsa', 'sce', ...) to the corresponding kegg
        identifier (T + 5 digit number).

        """
        for org in self._org_list:
            if org.org_code == code:
                return org.entry_id
        else:
            raise ValueError("Unknown organism code '%s'" % code)

    def search(self, string, relevance=False):
        """
        Search the genome database for string using ``bfind``.
        """
        if relevance:
            raise NotImplementedError("relevance is no longer supported")

        if string in self.TAXID_MAP:
            string = self.TAXID_MAP[string]

        res = self.api.find(self.DB, string).strip()
        if not res:
            return []

        res = res.splitlines()
        res = [r.split(",", 1)[0] for r in res]
        res = [r.split(None, 1)[1] for r in res]
        return res


@entry.entry_decorate
class GeneEntry(entry.DBEntry):
    FIELDS = [
        ("ENTRY", fields.DBEntryField),
        ("NAME", fields.DBNameField),
        ("DEFINITION", fields.DBDefinitionField),
        ("ORTHOLOGY", fields.DBSimpleField),
        ("ORGANISM", fields.DBSimpleField),
        ("PATHWAY", fields.DBPathway),
        ("MODULE", fields.DBSimpleField),
        ("BRITE", fields.DBSimpleField),
        ("DISEASE", fields.DBSimpleField),
        ("DRUG_TARGET", fields.DBSimpleField),
        ("CLASS", fields.DBSimpleField),
        ("MOTIF", fields.DBSimpleField),
        ("DBLINKS", fields.DBDBLinks),
        ("STRUCTURE", fields.DBSimpleField),
        ("POSITION", fields.DBSimpleField),
        ("AASEQ", fields.DBAASeq),
        ("NTSEQ", fields.DBNTSeq),
    ]

    def aliases(self):
        return (
            [self.entry_key]
            + (self.name.split(",") if self.name else [])
            + ([link[1][0] for link in self.dblinks.items()] if self.dblinks else [])
        )

    @property
    def alt_names(self):
        """
        For backwards compatibility.
        """
        return self.aliases()


class Genes(DBDataBase):
    """
    Interface to the KEGG Genes database.

    :param str org_code: KEGG organism code (e.g. 'hsa').

    """

    DB = None  # Needs to be set in __init__
    ENTRY_TYPE = GeneEntry

    def __init__(self, org_code):
        # TODO: Map to org code from kegg id (T + 5 digits)
        self.DB = org_code
        self.org_code = org_code
        DBDataBase.__init__(self)
        self._keys = self.api.get_genes_by_organism(org_code)

    def gene_aliases(self):
        aliases = {}
        for entry in self.itervalues():
            aliases.update(dict.fromkeys(entry.aliases(), self.org_code + ":" + entry.entry_key))

        return aliases


@entry.entry_decorate
class CompoundEntry(entry.DBEntry):
    FIELDS = [
        ("ENTRY", fields.DBEntryField),
        ("NAME", fields.DBNameField),
        ("FORMULA", fields.DBSimpleField),
        ("EXACT_MASS", fields.DBSimpleField),
        ("MOL_WEIGHT", fields.DBSimpleField),
        ("REMARK", fields.DBSimpleField),
        ("COMMENT", fields.DBSimpleField),
        ("REACTION", fields.DBSimpleField),
        ("PATHWAY", fields.DBPathway),
        ("ENZYME", fields.DBSimpleField),
        ("BRITE", fields.DBSimpleField),
        ("REFERENCE", fields.DBSimpleField),
        ("DBLINKS", fields.DBDBLinks),
        ("ATOM", fields.DBSimpleField),
        ("BOND", fields.DBSimpleField),
    ]


class Compound(DBDataBase):
    DB = "cpd"
    ENTRY_TYPE = CompoundEntry

    def __init__(self):
        DBDataBase.__init__(self)
        self._keys = [d.entry_id for d in self.api.list("cpd")]


@entry.entry_decorate
class ReactionEntry(entry.DBEntry):
    FIELDS = [
        ("ENTRY", fields.DBEntryField),
        ("NAME", fields.DBNameField),
        ("DEFINITION", fields.DBDefinitionField),
        ("EQUATION", fields.DBSimpleField),
        ("ENZYME", fields.DBSimpleField),
    ]


class Reaction(DBDataBase):
    DB = "rn"
    ENTRY_TYPE = ReactionEntry

    def __init__(self):
        DBDataBase.__init__(self)
        self._keys = [d.entry_id for d in self.api.list("rn")]


class Brite(DBDataBase):
    DB = "br"


class Disease(DBDataBase):
    DB = "ds"


class Drug(DBDataBase):
    DB = "dr"


@entry.entry_decorate
class EnzymeEntry(entry.DBEntry):
    FIELDS = [
        ("ENTRY", fields.DBEntryField),
        ("NAME", fields.DBNameField),
        ("CLASS", fields.DBSimpleField),
        ("SYSNAME", fields.DBSimpleField),
        ("REACTION", fields.DBSimpleField),
        ("ALL_REAC", fields.DBSimpleField),
        ("SUBSTRATE", fields.DBSimpleField),
        ("PRODUCT", fields.DBSimpleField),
        ("COMMENT", fields.DBSimpleField),
        ("REFERENCE", fields.DBReference),
        ("PATHWAY", fields.DBPathway),
        ("ORTHOLOGY", fields.DBSimpleField),
        ("GENES", fields.DBSimpleField),
        ("DBLINKS", fields.DBDBLinks),
    ]

    MULTIPLE_FIELDS = ["REFERENCE"]


class Enzyme(DBDataBase):
    DB = "ec"
    ENTRY_TYPE = EnzymeEntry

    def __init__(self):
        DBDataBase.__init__(self)
        self._keys = [d.entry_id for d in self.api.list("ec")]


@entry.entry_decorate
class OrthologyEntry(entry.DBEntry):
    FIELDS = [
        ("ENTRY", fields.DBEntryField),
        ("NAME", fields.DBNameField),
        ("CLASS", fields.DBSimpleField),
        ("DBLINKS", fields.DBDBLinks),
        ("GENES", fields.DBSimpleField),
    ]


class Orthology(DBDataBase):
    DB = "ko"
    ENTRY_TYPE = OrthologyEntry

    def __init__(self):
        DBDataBase.__init__(self)
        self._keys = [d.entry_id for d in self.api.list("ko")]


@entry.entry_decorate
class PathwayEntry(entry.DBEntry):
    FIELDS = [
        ("ENTRY", fields.DBEntryField),
        ("NAME", fields.DBNameField),
        ("DESCRIPTION", fields.DBSimpleField),
        ("CLASS", fields.DBSimpleField),
        ("PATHWAY_MAP", fields.DBPathwayMapField),
        ("MODULE", fields.DBSimpleField),
        ("DISEASE", fields.DBSimpleField),
        ("DRUG", fields.DBSimpleField),
        ("DBLINKS", fields.DBDBLinks),
        ("ORGANISM", fields.DBSimpleField),
        ("GENE", fields.DBGeneField),
        ("ENZYME", fields.DBEnzymeField),
        ("COMPOUND", fields.DBCompoundField),
        ("REFERENCE", fields.DBReference),
        ("REL_PATHWAY", fields.DBSimpleField),
        ("KO_PATHWAY", fields.DBSimpleField),
    ]

    MULTIPLE_FIELDS = ["REFERENCE"]

    @property
    def gene(self):
        if hasattr(self, "GENE"):
            genes = self.GENE._convert()
        else:
            return None

        org = self.organism
        org_prefix = ""
        if org:
            match = re.findall(r"\[GN:([a-z]+)\]", org)
            if match:
                org_prefix = match[0] + ":"
        genes = [org_prefix + g for g in genes]
        return genes


class Pathway(DBDataBase):
    """
    KEGG Pathway database

    :param str prefix:
        KEGG Organism code ('hsa', ...) or 'map', 'ko', 'ec' or 'rn'

    """

    DB = "path"
    ENTRY_TYPE = PathwayEntry

    def __init__(self, prefix="map"):
        DBDataBase.__init__(self)
        self.prefix = prefix
        valid = [d.org_code for d in self.api.list_organisms()] + ["map", "ko", "ec", "rn"]

        if prefix not in valid:
            raise ValueError("Invalid prefix %r" % prefix)

        self._keys = [d.entry_id for d in self.api.list("pathway/" + prefix)]
