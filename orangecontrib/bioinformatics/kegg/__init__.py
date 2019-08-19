"""
==============================================
KEGG - Kyoto Encyclopedia of Genes and Genomes
==============================================

:mod:`kegg` is a python module for accessing `KEGG (Kyoto Encyclopedia
of Genes and Genomes) <http://www.genome.jp/kegg/>`_ using its web services.

.. note:: This module requires `slumber`_ and `requests`_ packages.

.. _`slumber`: https://pypi.python.org/pypi/slumber/

.. _`requests`: https://pypi.python.org/pypi/requests

>>> # Create a KEGG Genes database interface
>>> genome = KEGGGenome()
>>> # List all available entry ids
>>> keys = list(genome.keys())
>>> print(keys[0])
T01001
>>> # Retrieve the entry for the key.
>>> entry = genome[keys[0]]
>>> print(entry.entry_key)
T01001
>>> print(entry.definition)
Homo sapiens (human)
>>> print(entry)  # doctest: +SKIP
ENTRY       T01001            Complete  Genome
NAME        hsa, HUMAN, 9606
DEFINITION  Homo sapiens (human)
...

The :class:`Organism` class can be a convenient starting point
for organism specific databases.

>>> organism = Organism("Homo sapiens")  # searches for the organism by name
>>> print(organism.org_code)  # prints the KEGG organism code
hsa
>>> genes = organism.genes  # get the genes database for the organism
>>> gene_ids = list(genes.keys()) # KEGG gene identifiers
>>> entry = genes["hsa:672"]
>>> print(entry.definition) # doctest: +SKIP
(RefSeq) BRCA1, DNA repair associated
>>> # print the entry in DBGET database format.
>>> print(entry) # doctest: +SKIP
ENTRY       672               CDS       T01001
NAME        BRCA1, BRCAI, BRCC1, BROVCA1, FANCS, IRIS, PNCA4, PPP1R53, PSCP, RNF53
DEFINITION  ...
"""
from __future__ import absolute_import

import os
import sys
import threading
from datetime import datetime
from functools import reduce
from itertools import chain
from contextlib import contextmanager
from collections import defaultdict

from orangecontrib.bioinformatics.kegg import api, conf, entry, pathway, databases
from orangecontrib.bioinformatics.ncbi import taxonomy
from orangecontrib.bioinformatics.utils import statistics, progress_bar_milestones
from orangecontrib.bioinformatics.kegg.brite import Brite, BriteEntry

KEGGGenome = databases.Genome
KEGGGenes = databases.Genes
KEGGEnzyme = databases.Enzyme
KEGGReaction = databases.Reaction
KEGGPathways = databases.Pathway
KEGGCompound = databases.Compound

KEGGBrite = Brite
KEGGBriteEntry = BriteEntry

KEGGPathway = pathway.Pathway

DEFAULT_CACHE_DIR = conf.params["cache.path"]


class Organism(object):
    """
    A convenience class for retrieving information regarding an
    organism in the KEGG Genes database.

    :param str org: KEGG organism code (e.g. "hsa", "sce"). Can also be a
        descriptive name (e.g. 'yeast', "homo sapiens") in which case the
        organism code will be searched for by using KEGG `find` api.

    .. seealso::

        :func:`organism_name_search`
            Search KEGG for an organism code

    """

    def __init__(self, org):
        self.org_code = self.organism_name_search(org)
        self.api = api.CachedKeggApi()

    @property
    def org(self):
        """
        KEGG organism code.
        """
        return self.org_code

    @property
    def genes(self):
        """
        An :class:`~.databases.Genes` database instance for this organism.
        """
        # TODO: This should not be a property but a method.
        # I think it was only put here as back compatibility with old obiKEGG.
        if not hasattr(self, "_genes"):
            genes = KEGGGenes(self.org_code)
            self._genes = genes
        return self._genes

    def kegg_to_ncbi_mapper(self):
        return {kegg.upper(): ncbi.split(':')[-1] for ncbi, kegg in self.api.conv(self.org_code, "ncbi-geneid")}

    def get_ncbi_ids(self):
        # return [ncbi.split(':')[-1] for ncbi, _ in self.api.conv(self.org_code, "ncbi-geneid")]
        return [ncbi for ncbi in self.kegg_to_ncbi_mapper().values()]

    def gene_aliases(self):
        """
        Return a list of sets of equal genes (synonyms) in KEGG for
        this organism.

        .. note::

            This only includes 'ncbi-geneid' and 'ncbi-proteinid' records
            from the KEGG Genes DBLINKS entries.

        """
        definitions = self.api.list(self.org_code)
        ncbi_geneid = self.api.conv(self.org_code, "ncbi-geneid")
        ncbi_gi = self.api.conv(self.org_code, "ncbi-proteinid")

        aliases = defaultdict(set)

        for entry_id, definition in definitions:
            # genes entry id without the organism code
            aliases[entry_id].add(entry_id.split(":", 1)[1])
            # all names in the NAME field (KEGG API list returns
            # 'NAME; DEFINITION') fields for genes
            names = definition.split(";")[0].split(",")
            aliases[entry_id].update([name.strip() for name in names])

        for source_id, target_id in chain(ncbi_geneid, ncbi_gi):
            aliases[target_id].add(source_id.split(":", 1)[1])

        return [set([entry_id]).union(names) for entry_id, names in aliases.items()]

    def pathways(self, with_ids=None):
        """
        Return a list of all pathways for this organism.
        """
        if with_ids is not None:
            return self.api.get_pathways_by_genes(with_ids)
        else:
            return [p.entry_id for p in self.api.list_pathways(self.org_code)]

    def list_pathways(self):
        """
        List all pathways for this organism.

        .. deprecated: 2.5
            Use :func:`pathways` instead.

        """
        # NOTE: remove/deprecate and use pathways()
        return self.pathways()

    def get_linked_pathways(self, pathway_id):
        self.api.get_linked_pathways(pathway_id)

    def enzymes(self, genes=None):
        raise NotImplementedError()

    def get_enriched_pathways(self, genes, reference=None, prob=statistics.Binomial(), callback=None):
        """
        Return a dictionary with enriched pathways ids as keys
        and (list_of_genes, p_value, num_of_reference_genes) tuples
        as items.

        """
        if reference is None:
            reference = self.genes.keys()
        reference = set(reference)

        allPathways = defaultdict(lambda: [[], 1.0, []])
        milestones = progress_bar_milestones(len(genes), 100)
        pathways_db = KEGGPathways()

        pathways_for_gene = []
        for i, gene in enumerate(genes):
            pathways_for_gene.append(self.pathways([gene]))
            if callback and i in milestones:
                callback(i * 50.0 / len(genes))

        # pre-cache for speed
        pathways_db.pre_cache([pid for pfg in pathways_for_gene for pid in pfg])
        for i, (gene, pathways) in enumerate(zip(genes, pathways_for_gene)):
            for pathway in pathways:
                if pathways_db.get_entry(pathway).gene:
                    allPathways[pathway][0].append(gene)
            if callback and i in milestones:
                callback(50.0 + i * 50.0 / len(genes))

        pItems = allPathways.items()

        for i, (p_id, entry) in enumerate(pItems):
            pathway = pathways_db.get_entry(p_id)
            entry[2].extend(reference.intersection(pathway.gene or []))
            entry[1] = prob.p_value(len(entry[0]), len(reference), len(entry[2]), len(genes))
        return dict([(pid, (genes, p, len(ref))) for pid, (genes, p, ref) in allPathways.items()])

    def get_genes_by_enzyme(self, enzyme):
        enzyme = KEGGEnzyme().get_entry(enzyme)
        return enzyme.genes.get(self.org_code, []) if enzyme.genes else []

    def get_genes_by_pathway(self, pathway_id):
        # print(len(KEGGPathway(pathway_id).genes()), KEGGPathway(pathway_id).genes()[1])
        return KEGGPathway(pathway_id).genes()

    def get_enzymes_by_pathway(self, pathway_id):
        return KEGGPathway(pathway_id).enzymes()

    def get_compounds_by_pathway(self, pathway_id):
        return KEGGPathway(pathway_id).compounds()

    def get_pathways_by_genes(self, gene_ids):
        """ Pathways that include all genes in gene_ids. """
        l = self.api.get_genes_pathway_organism(self.org_code)
        gene_ids = set(gene_ids)
        gtp = defaultdict(set)
        for a, b in l:
            gtp[a].add(b)
        pathways = [gtp[g] for g in gene_ids]
        pathways = reduce(set.intersection, pathways)
        return sorted(pathways)

    def get_pathways_by_enzymes(self, enzyme_ids):
        enzyme_ids = set(enzyme_ids)
        pathways = [KEGGEnzyme()[id].pathway for id in enzyme_ids]
        pathways = reduce(set.union, pathways, set())
        return [id for id in pathways if enzyme_ids.issubset(KEGGPathway(id).enzymes())]

    def get_pathways_by_compounds(self, compound_ids):
        compound_ids = set(compound_ids)
        pathways = [KEGGCompound()[id].pathway for id in compound_ids]
        pathways = reduce(set.union, pathways, set())
        return [id for id in pathways if compound_ids.issubset(KEGGPathway(id).compounds())]

    def get_enzymes_by_compound(self, compound_id):
        return KEGGCompound()[compound_id].enzyme

    def get_enzymes_by_gene(self, gene_id):
        return self.genes[gene_id].enzymes

    def get_compounds_by_enzyme(self, enzyme_id):
        return self._enzymes_to_compounds.get(enzyme_id)

    def get_genes(self):
        return self.genes

    @classmethod
    def organism_name_search(cls, name):
        return organism_name_search(name)

    @classmethod
    def organism_version(cls, name):
        name = cls.organism_name_search(name)
        with _global_genome_instance() as genome:
            info = genome.api.info(name)
            return info.releas


KEGGOrganism = Organism


def organism_name_search(name):
    """
    Search for a organism by `name` and return it's KEGG organism code.
    """
    with _global_genome_instance() as genome:
        try:
            name = genome.org_code_to_entry_key(name)
        except ValueError:
            pass

        if name not in genome:
            ids = genome.search(name)
            if ids:
                name = ids.pop(0) if ids else name
            else:
                raise taxonomy.UnknownSpeciesIdentifier(name)

        try:
            return genome[name].organism_code
        except Exception as e:
            raise taxonomy.UnknownSpeciesIdentifier(name)


def pathways(org):
    """
    Return a list of all KEGG pathways for an KEGG organism code `org`.
    """
    return KEGGPathway.list(org)


def from_taxid(taxid):
    """
    Return a KEGG organism code for a an NCBI Taxonomy id string `taxid`.
    """
    with _global_genome_instance() as genome:
        res = genome.search(taxid)
        for r in res:
            e = genome[r]

            if e.taxid in [taxid, genome.TAXID_MAP.get(taxid, taxid)]:
                return e.org_code()

        return None


def to_taxid(name):
    """
    Return a NCBI Taxonomy id for a given KEGG Organism name
    """
    with _global_genome_instance() as genome:

        if name in genome:  # a T string
            return genome[name].taxid

        name2 = genome.org_code_to_entry_key(name)
        if name2 in genome:
            return genome[name2].taxid

        keys = genome.search(name)
        if keys:
            return genome[keys[0]].taxid
        else:
            return None


def _global_genome_instance():
    genome = getattr(_global_genome_instance, "_genome", None)
    if genome is None:
        genome = KEGGGenome()
        genome._lock = threading.RLock()
        _global_genome_instance._genome = genome

    @contextmanager
    def instance_locked(instance, lock):
        lock.acquire()
        try:
            yield instance
        finally:
            lock.release()

    return instance_locked(genome, genome._lock)
