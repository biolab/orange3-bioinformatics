""" Mutant Phenotypes """
import json
from functools import reduce
from collections import defaultdict

from orangecontrib.bioinformatics.utils import serverfiles

DOMAIN = 'dictybase'
PHENOTYPES_FILENAME = 'mutant_phenotypes.json'


class DictyMutant:
    __slots__ = ('name', 'descriptor', 'genes', 'phenotypes', 'mutant_types')

    def __init__(self, mutant_entry: dict) -> None:
        """ A single Dictyostelium discoideum mutant from the Dictybase.

        :param mutant_entry: A single mutant entry from `curated mutants file
            <http://dictybase.org/db/cgi-bin/dictyBase/download/download.pl?area=mutant_phenotypes&ID=all-mutants.txt>`_.

        """
        # dictyBase ID
        self.name = mutant_entry.get('systematic_name', None)

        # dictyBase strain descriptor
        self.descriptor = mutant_entry.get('strain_descriptor', None)

        # all associated genes
        self.genes = mutant_entry.get('associated_genes', [])

        # all associated phenotypes
        self.phenotypes = mutant_entry.get('phenotypes', [])

        self.mutant_types = mutant_entry.get('mutant_types', None)

    def __repr__(self):
        return f'<dictyBase ID={self.name}, strain descriptor={self.descriptor}>'


class DictyMutants:
    DEFAULT_DATABASE_PATH = serverfiles.localpath(DOMAIN)  # use a default local folder for storing the genesets

    def __init__(self, file_path=None):
        """  A collection of Dictybase mutants as a dictionary of :obj:`DictyMutant` objects.
        """
        if file_path is None:
            file_path = serverfiles.localpath_download(DOMAIN, PHENOTYPES_FILENAME)

        with open(file_path, 'r') as fp:
            _mutants = [DictyMutant(mutant) for mutant in json.load(fp)]
            self._mutants = {m: m for m in _mutants}

    @classmethod
    def get_instance(cls):
        if not hasattr(cls, "_shared_dict"):
            dicty = DictyMutants()
            cls._shared_dict = dicty.__dict__
        instance = DictyMutants.__new__(DictyMutants)
        instance.__dict__ = cls._shared_dict
        return instance

    def mutants(self):
        return list(self._mutants.keys())

    def genes(self):
        return sorted(set(reduce(list.__add__, [self.mutant_genes(mutant) for mutant in self.mutants()], [])))

    def phenotypes(self):
        return sorted(set(reduce(list.__add__, [self.mutant_phenotypes(mutant) for mutant in self.mutants()], [])))

    def mutant_genes(self, mutant):
        return self._mutants[mutant].genes

    def mutant_phenotypes(self, mutant):
        return self._mutants[mutant].phenotypes

    def gene_mutants(self):
        dgm = defaultdict(set)
        for mutant, genes in [(mutant, self.mutant_genes(mutant)) for mutant in self.mutants()]:
            for gene in genes:
                dgm[gene].add(mutant)
        return dgm

    def phenotype_mutants(self):
        dpm = defaultdict(set)
        for mutant, phenotypes in [(mutant, self.mutant_phenotypes(mutant)) for mutant in self.mutants()]:
            for phenotype in phenotypes:
                dpm[phenotype].add(mutant)
        return dpm


def mutants():
    """ Return all :obj:`DictyMutant` objects.
    """
    return DictyMutants.get_instance().mutants()


def genes():
    """ Return a set of all genes referenced in the Dictybase.
    """
    return DictyMutants.get_instance().genes()


def phenotypes():
    """ Return a set of all phenotypes referenced in Dictybase.
    """
    return DictyMutants.get_instance().phenotypes()


def mutant_genes(mutant):
    """ Return a set of all genes referenced by a `mutant` in Dictybase.
    """
    return DictyMutants.get_instance().mutant_genes(mutant)


def mutant_phenotypes(mutant):
    """ Return a set of all phenotypes referenced by a `mutant` in Dictybase.
    """
    return DictyMutants.get_instance().mutant_phenotypes(mutant)


def gene_mutants():
    """ Return a dictionary { gene: set(mutant_objects for mutant), ... }.
    """
    return DictyMutants.get_instance().gene_mutants()


def phenotype_mutants():
    """ Return a dictionary { phenotype: set(mutant_objects for mutant), ... }.
    """
    return DictyMutants.get_instance().phenotype_mutants()


if __name__ == "__main__":
    dicty_mutants = mutants()
    mutant = list(dicty_mutants)[0]
    print(mutant.name, mutant.descriptor, mutant.genes)
    print(gene_mutants())
