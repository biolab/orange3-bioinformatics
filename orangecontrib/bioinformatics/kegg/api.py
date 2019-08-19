""" KEGG api interface. """
from __future__ import absolute_import

import os
import warnings
from datetime import datetime
from operator import itemgetter
from contextlib import closing

import six

from orangecontrib.bioinformatics.kegg import caching
from orangecontrib.bioinformatics.kegg.types import Link, BInfo, Definition, OrganismSummary
from orangecontrib.bioinformatics.kegg.caching import touch_dir, cache_entry, cached_method
from orangecontrib.bioinformatics.kegg.service import web_service

# A list of all databases with names, abbreviations
DATABASES = [
    ("KEGG Pathway", "pathway", "path", None),
    ("KEGG Brite", "brite", "br", None),
    ("KEGG Module", "module", "md", "M"),
    ("KEGG Disease", "disease", "ds", "H"),
    ("KEGG Drug", "drug", "dr", "D"),
    ("KEGG Orthology", "orthology", "ko", "K"),
    ("KEGG Genome", "genome", "genome", "T"),
    ("KEGG Genomes", "genomes", "gn", "T"),
    ("KEGG Genes", "genes", None, None),
    ("KEGG Ligand", "ligand", "ligand", None),
    ("KEGG Compound", "compound", "cpd", "C"),
    ("KEGG Glycan", "glycan", "gl", "G"),
    ("KEGG Reaction", "reaction", "rn", "R"),
    ("KEGG RPair", "rpair", "rp", "RP"),
    ("KEGG RClass", "rclass", "rc", "RC"),
    ("KEGG Enzyme", "enzyme", "ec", "E"),
]


def _link_targets(links):
    return sorted(set(map(itemgetter(1), links)))


class KeggApi(object):
    """
    An abstraction of a rest KEGG API.
    """

    def __init__(self):
        self.service = web_service()

    def list_organisms(self):
        """
        Return a list of all available organisms,

        >>> api.list_organisms()  # doctest: +ELLIPSIS
        [OrganismSummary(entry_id='T01001', ...
        """
        return list(map(OrganismSummary.from_str, self.service.list.organism.get().splitlines()))

    def list_pathways(self, organism):
        """
        Return a list of all available pathways for `organism`

        >>> api.list_pathways("hsa")  # doctest: +ELLIPSIS
        [Definition(entry_id='path:hsa00010', ...
        """
        return list(map(Definition.from_str, self.service.list.pathway(organism).get().splitlines()))

    def list(self, db):
        """
        Return a list of all available entries in database `db`.
        """
        return list(map(Definition.from_str, self.service.list(db).get().splitlines()))

    #######
    # DBGET
    #######

    def info(self, db):
        """
        Return info for database `db`

        >>> print(api.info("pathway"))
        BInfo(entry_id='path', definition='KEGG Pathway Database', ...
        """
        result = self.service.info(db).get()
        return BInfo.from_text(result)

    def find(self, db, keywords):
        """
        Search database 'db' for keywords.
        """
        if isinstance(keywords, six.string_types):
            keywords = [keywords]

        return self.service.find(db)("+".join(keywords)).get()

    def get(self, ids):
        """
        Retrieve database entries for `ids` list.
        """
        if not isinstance(ids, six.string_types):
            # Sequence of ids
            ids = "+".join(ids)

        return self.service.get(ids).get()

    def conv(self, target_db, source):
        """
        Return a mapping from source to target_db ids as a list of two
        tuples [(source_id, target_id), ...].

        """
        if not isinstance(source, six.string_types):
            source = "+".join(source)

        res = self.service.conv(target_db)(source).get()
        return [tuple(line.split("\t")) for line in res.splitlines()]

    def link(self, target_db, source_db=None, ids=None):
        if not (source_db or ids):
            raise ValueError("One of 'source_db' or 'ids' must be supplied")
        if source_db and ids:
            raise ValueError("Only one 'source_db' or 'ids' must be supplied")

        if source_db:
            result = self.service.link(target_db)(source_db).get()
        else:
            result = self.service.link(target_db)("+".join(ids)).get()

        return list(map(Link._make, map(str.split, result.splitlines())))

    def get_genes_by_enzyme(self, enzyme_id, org):
        return _link_targets(self.link(org, ids=[enzyme_id]))

    def get_enzymes_by_gene(self, gene_id):
        return _link_targets(self.link("ec", ids=[gene_id]))

    def get_enzymes_by_compound(self, compound_id):
        return _link_targets(self.link("ec", ids=[compound_id]))

    def get_enzymes_by_glycan(self, glycan_id):
        return _link_targets(self.link("ec", ids=[glycan_id]))

    def get_enzymes_by_reaction(self, reaction_id):
        return _link_targets(self.link("ec", ids=[reaction_id]))

    def get_compounds_by_enzyme(self, enzyme_id):
        return _link_targets(self.link("compound", ids=[enzyme_id]))

    def get_compounds_by_reaction(self, reaction_id):
        return _link_targets(self.link("compound", ids=[reaction_id]))

    def get_glycans_by_enzyme(self, enzyme_id):
        return _link_targets(self.link("gl", ids=[enzyme_id]))

    def get_glycans_by_reaction(self, reaction_id):
        return _link_targets(self.link("gl", ids=[reaction_id]))

    def get_reactions_by_enzyme(self, enzyme_id):
        return _link_targets(self.link("rn", ids=[enzyme_id]))

    def get_reactions_by_compound(self, compound_id):
        return _link_targets(self.link("rn", ids=[compound_id]))

    def get_reactions_by_glycan(self, glycan_id):
        return _link_targets(self.link("rn", ids=[glycan_id]))

    ######
    # SSDB
    ######

    # No replacement api in the KEGG REST api.
    def get_best_best_neighbors_by_gene(self, genes_id, offset, limit):
        raise NotImplementedError

    def get_best_neighbors_by_gene(self, genes_id, offset, limit):
        raise NotImplementedError

    def get_reverse_best_neighbors_by_gene(self, genes_id, offset, limit):
        raise NotImplementedError

    def get_paralogs_by_gene(self, genes_id, offset, limit):
        raise NotImplementedError

    #######
    # Motif
    #######

    # No replacement api in KEGG REST api
    def get_motifs_by_gene(self, genes_id, db):
        raise NotImplementedError

    def get_genes_by_motifs(self, motif_id_list, offset, limit):
        raise NotImplementedError

    ####
    # KO
    ####

    def get_ko_by_gene(self, genes_id):
        raise NotImplementedError

    def get_ko_by_ko_class(self, ko_class_id):
        raise NotImplementedError

    def get_genes_by_ko_class(self, ko_class_id, org, offset, limit):
        raise NotImplementedError

    def get_genes_by_ko(self, ko_id, org):
        raise NotImplementedError

    #########
    # Pathway
    #########

    def mark_pathway_by_objects(self, pathway_id, object_id_list):
        raise NotImplementedError

    def color_pathway_by_objects(self, pathway_id, object_id_list, fg_color_list, bg_color_list):
        raise NotImplementedError

    def color_pathway_by_elements(self, pathway_id, element_id_list, fg_color_list, bg_color_list):
        raise NotImplementedError

    def get_html_of_marked_pathway_by_objects(self, pathway_id, object_id_list):
        raise NotImplementedError

    def get_html_of_colored_pathway_by_objects(self, pathway_id, object_id_list, fg_color_list, bg_color_list):
        raise NotImplementedError

    def get_html_of_colored_pathway_by_elements(self, pathway_id, element_id_list, fg_color_list, bg_color_list):
        raise NotImplementedError

    def get_references_by_pathway(self, pathway_id):
        return self.service.get_references_by_pathway(pathway_id)

    def get_element_relations_by_pathway(self, pathway_id):
        return self.service.get_element_relations_by_pathway(pathway_id)

    def get_genes_by_organism(self, organism, offset=None, limit=None):
        if offset is not None:
            raise NotImplementedError("offset is no longer supported")
        if limit is not None:
            raise NotImplementedError("limit is no longer supported.")

        res = self.service.list(organism).get().splitlines()
        return [r.split(None, 1)[0] for r in res]

    def get_number_of_genes_by_organism(self, organism):
        raise NotImplementedError

    ####################
    # Objects by pathway
    ####################

    def get_elements_by_pathway(self, pathway_id):
        raise NotImplementedError

    def get_genes_by_pathway(self, pathway_id):
        return _link_targets(self.link("genes", ids=[pathway_id]))

    def get_enzymes_by_pathway(self, pathway_id):
        return _link_targets(self.link("ec", ids=[pathway_id]))

    def get_compounds_by_pathway(self, pathway_id):
        return _link_targets(self.link("compound", ids=[pathway_id]))

    def get_drugs_by_pathway(self, pathway_id):
        return _link_targets(self.link("drug", ids=[pathway_id]))

    def get_glycans_by_pathway(self, pathway_id):
        return _link_targets(self.link("gl", ids=[pathway_id]))

    def get_reactions_by_pathway(self, pathway_id):
        return _link_targets(self.link("rn", ids=[pathway_id]))

    def get_kos_by_pathway(self, pathway_id):
        return _link_targets(self.link("ko", ids=[pathway_id]))

    #############################################
    # Pathways and genes of a specific organism #
    #############################################

    def get_genes_pathway_organism(self, organism):
        l = self.link("pathway", organism)
        return list(map(tuple, l))

    #####################
    # Pathways by objects
    #####################

    # These functions returned results intersections.
    def get_pathways_by_genes(self, gene_list):
        raise NotImplementedError

    def get_pathways_by_enzymes(self, enzyme_list):
        raise NotImplementedError

    def get_pathways_by_compounds(self, compound_list):
        raise NotImplementedError

    def get_pathways_by_drugs(self, drug_list):
        raise NotImplementedError

    def get_pathways_by_glycans(self, glycan_list):
        raise NotImplementedError

    def get_pathways_by_reactions(self, reaction_list):
        raise NotImplementedError

    def get_pathways_by_kos(self, ko_list):
        raise NotImplementedError

    ##########################
    # Relations among pathways
    ##########################

    def get_linked_pathways(self, pathway_id):
        if not pathway_id.startswith("path:"):
            pathway_id = "path:" + pathway_id
        return _link_targets(self.link("pathway", ids=[pathway_id]))


"""
KEGG api with caching
"""


try:
    from functools import lru_cache
except ImportError:
    # TODO: move a copy of lru_cache in .caching if distributing this as a
    # standalone package
    from Orange.utils import lru_cache


class CachedKeggApi(KeggApi):
    def __init__(self, store=None):
        KeggApi.__init__(self)
        if store is None:
            self.store = {}

    # Needed API for cached decorator.
    def cache_store(self):
        from . import conf

        path = conf.params["cache.path"]
        touch_dir(path)
        return caching.Sqlite3Store(os.path.join(path, "kegg_api_cache_2.sqlite3"))

    def last_modified(self, args, kwargs=None):
        return getattr(self, "default_release", "")

    def set_default_release(self, release):
        self.default_release = release

    @cached_method
    def list_organisms(self):
        return KeggApi.list_organisms(self)

    @cached_method
    def list_pathways(self, organism):
        return KeggApi.list_pathways(self, organism)

    @cached_method
    def list(self, db):
        return KeggApi.list(self, db)

    @lru_cache()  # not persistently cached
    def info(self, db):
        return KeggApi.info(self, db)

    @cached_method
    def find(self, db, keywords):
        return KeggApi.find(self, db, keywords)

    @cached_method
    def get(self, ids):
        if not isinstance(ids, six.string_types):
            return self._batch_get(ids)
        else:
            return KeggApi.get(self, ids)

    @cached_method
    def link(self, target_db, source_db=None, ids=None):
        return KeggApi.link(self, target_db, source_db, ids)

    def _batch_get(self, ids):
        if len(ids) > 10:
            raise ValueError("Can batch at most 10 ids at a time.")

        get = self.get
        uncached = []
        unmatched = set()

        with closing(get.cache_store()) as store:
            # Which ids are already cached
            # TODO: Invalidate entries by release string.
            for id in ids:
                key = get.key_from_args((id,))
                if not get.key_has_valid_cache(key, store):
                    uncached.append(id)

        if uncached:
            # in case there are duplicate ids
            uncached = sorted(set(uncached))

            rval = KeggApi.get(self, uncached)

            if rval is not None:
                entries = rval.split("///\n")
            else:
                entries = []

            if entries and not entries[-1].strip():
                # Delete the last single newline entry if present
                del entries[-1]

            if len(entries) != len(uncached):
                new_uncached, entries = match_by_ids(uncached, entries)
                unmatched = set(uncached) - set(new_uncached)
                uncached = new_uncached
                warnings.warn("Unable to match entries for keys: %s." % ", ".join(map(repr, unmatched)))

            with closing(get.cache_store()) as store:
                for id, entry in zip(uncached, entries):
                    key = get.key_from_args((id,))
                    if entry is not None:
                        entry = entry + "///\n"
                    store[key] = cache_entry(entry, mtime=datetime.now())

        # Finally join all the results, but drop all None objects

        with closing(get.cache_store()):
            keys = [get.key_from_args((id,)) for id in ids]
            entries = [store[key].value for key in keys]

        entries = filter(lambda e: e is not None, entries)

        rval = "".join(entries)
        return rval

    @cached_method
    def conv(self, target_db, source):
        return KeggApi.conv(self, target_db, source)

    ########
    # LinkDB
    ########

    @cached_method
    def get_genes_by_enzyme(self, enzyme_id, org):
        return KeggApi.get_genes_by_enzyme(self, enzyme_id, org)

    @cached_method
    def get_enzymes_by_gene(self, genes_id):
        return KeggApi.get_enzymes_by_gene(self, genes_id)

    @cached_method
    def get_enzymes_by_compound(self, compound_id):
        return KeggApi.get_enzymes_by_compound(self, compound_id)

    @cached_method
    def get_enzymes_by_glycan(self, glycan_id):
        return KeggApi.get_enzymes_by_glycan(self, glycan_id)

    @cached_method
    def get_enzymes_by_reaction(self, reaction_id):
        return KeggApi.get_enzymes_by_reaction(self, reaction_id)

    @cached_method
    def get_compounds_by_enzyme(self, enzyme_id):
        return KeggApi.get_compounds_by_enzyme(self, enzyme_id)

    @cached_method
    def get_compounds_by_reaction(self, reaction_id):
        return KeggApi.get_compounds_by_reaction(self, reaction_id)

    @cached_method
    def get_glycans_by_enzyme(self, enzyme_id):
        return KeggApi.get_glycans_by_enzyme(self, enzyme_id)

    @cached_method
    def get_glycans_by_reaction(self, reaction_id):
        return KeggApi.get_glycans_by_reaction(self, reaction_id)

    @cached_method
    def get_reactions_by_enzyme(self, enzyme_id):
        return KeggApi.get_reactions_by_enzyme(self, enzyme_id)

    @cached_method
    def get_reactions_by_compound(self, compound_id):
        return KeggApi.get_reactions_by_compound(self, compound_id)

    @cached_method
    def get_reactions_by_glycan(self, glycan_id):
        return KeggApi.get_reactions_by_glycan(self, glycan_id)

    ######
    # SSDB
    ######

    @cached_method
    def get_best_best_neighbors_by_gene(self, genes_id, offset, limit):
        return KeggApi.get_best_best_neighbors_by_gene(self, genes_id, offset, limit)

    @cached_method
    def get_best_neighbors_by_gene(self, genes_id, offset, limit):
        return KeggApi.get_best_neighbors_by_gene(self, genes_id, offset, limit)

    @cached_method
    def get_reverse_best_neighbors_by_gene(self, genes_id, offset, limit):
        return KeggApi.get_reverse_best_neighbors_by_gene(self, genes_id, offset, limit)

    @cached_method
    def get_paralogs_by_gene(self, genes_id, offset, limit):
        return KeggApi.get_paralogs_by_gene(self, genes_id, offset, limit)

    #######
    # Motif
    #######

    @cached_method
    def get_motifs_by_gene(self, genes_id, db):
        return KeggApi.get_motifs_by_gene(self, genes_id, db)

    @cached_method
    def get_genes_by_motifs(self, motif_id_list, offset, limit):
        return KeggApi.get_genes_by_motifs(self, motif_id_list, offset, limit)

    ####
    # KO
    ####

    @cached_method
    def get_ko_by_gene(self, genes_id):
        return KeggApi.get_ko_by_gene(self, genes_id)

    @cached_method
    def get_ko_by_ko_class(self, ko_class_id):
        return KeggApi.service.get_ko_by_ko_class(self, ko_class_id)

    @cached_method
    def get_genes_by_ko_class(self, ko_class_id, org, offset, limit):
        return KeggApi.get_genes_by_ko_class(self, ko_class_id, org, offset, limit)

    @cached_method
    def get_genes_by_ko(self, ko_id, org):
        return KeggApi.get_genes_by_ko(self, ko_id, org)

    #########
    # Pathway
    #########

    @cached_method
    def get_genes_by_organism(self, organism, offset=None, limit=None):
        return KeggApi.get_genes_by_organism(self, organism, offset=offset, limit=limit)

    @cached_method
    def get_number_of_genes_by_organism(self, organism):
        return KeggApi.get_number_of_genes_by_organism(self, organism)

    @cached_method
    def get_pathways_by_genes(self, gene_list):
        return KeggApi.get_pathways_by_genes(self, gene_list)

    @cached_method
    def get_pathways_by_enzymes(self, enzyme_list):
        return KeggApi.get_pathways_by_enzymes(self, enzyme_list)

    @cached_method
    def get_pathways_by_compounds(self, compound_list):
        return KeggApi.get_pathways_by_compounds(self, compound_list)

    @cached_method
    def get_pathways_by_drugs(self, drug_list):
        return KeggApi.get_pathways_by_drugs(self, drug_list)

    @cached_method
    def get_pathways_by_glycans(self, glycan_list):
        return KeggApi.get_pathways_by_glycans(self, glycan_list)

    @cached_method
    def get_pathways_by_reactions(self, reaction_list):
        return KeggApi.get_pathways_by_reactions(self, reaction_list)

    @cached_method
    def get_pathways_by_kos(self, ko_list):
        return KeggApi.get_pathways_by_kos(self, ko_list)

    @cached_method
    def get_elements_by_pathway(self, pathway_id):
        return KeggApi.get_elements_by_pathway(self, pathway_id)

    @cached_method
    def get_genes_by_pathway(self, pathway_id):
        return KeggApi.get_genes_by_pathway(self, pathway_id)

    @cached_method
    def get_enzymes_by_pathway(self, pathway_id):
        return KeggApi.get_enzymes_by_pathway(self, pathway_id)

    @cached_method
    def get_compounds_by_pathway(self, pathway_id):
        return KeggApi.get_compounds_by_pathway(self, pathway_id)

    @cached_method
    def get_drugs_by_pathway(self, pathway_id):
        return KeggApi.get_drugs_by_pathway(self, pathway_id)

    @cached_method
    def get_glycans_by_pathway(self, pathway_id):
        return KeggApi.get_glycans_by_pathway(self, pathway_id)

    @cached_method
    def get_reactions_by_pathway(self, pathway_id):
        return KeggApi.get_reactions_by_pathway(self, pathway_id)

    @cached_method
    def get_kos_by_pathway(self, pathway_id):
        return KeggApi.get_kos_by_pathway(self, pathway_id)

    @cached_method
    def get_genes_pathway_organism(self, org):
        return KeggApi.get_genes_pathway_organism(self, org)


def match_by_ids(ids, entries):
    """

    """

    unmatched_ids = set(ids)
    unmatched_entries = set(entries)

    matched_ids = []
    matched_entries = []

    def match_add(search_id, entry):
        """
        Move search_id and entry to the matched lists.
        """
        matched_ids.append(search_id)
        matched_entries.append(entry)

        # Remove from the unmatched set
        unmatched_ids.remove(search_id)
        unmatched_entries.remove(entry)

    def entry_split(entry_text):
        line, _ = entry_text.split("\n", 1)
        return line.split(None, 2)

    entries_by_id = {}

    for entry in entries:
        _, eid, _ = entry_split(entry)
        entries_by_id[eid] = entry

    # First match full search ids
    for search_id in list(unmatched_ids):
        if search_id in entries_by_id:
            entry = entries_by_id.pop(search_id)
            match_add(search_id, entry)

    # Second pass, split the search ids by ':' to db and identifier part,
    # match by identifier
    for search_id in list(unmatched_ids):
        if ":" in search_id:
            db_id, rest = search_id.split(":", 1)
            if rest in entries_by_id:
                entry = entries_by_id.pop(rest)
                match_add(search_id, entry)

    return matched_ids, matched_entries
