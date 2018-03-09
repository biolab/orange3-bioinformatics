""" Mutant Phenotypes """
import os
import shutil
import pickle

from urllib.request import urlopen
from collections import defaultdict
from functools import reduce

from orangecontrib.bioinformatics.dicty.config import DOMAIN, PHENOTYPES_FILENAME
from orangecontrib.bioinformatics.utils import serverfiles


class DictyMutant:

    def __init__(self, mutant_entry):  # type: (str) -> None
        """ A single Dictyostelium discoideum mutant from the Dictybase.

        Args:
            mutant_entry: A single mutant entry from `curated mutants file
            <http://dictybase.org/db/cgi-bin/dictyBase/download/download.pl?area=mutant_phenotypes&ID=all-mutants.txt>`_.
        """

        mutant = mutant_entry.split("\t")

        # dictyBase ID
        self.name = mutant[0]

        # dictyBase strain descriptor
        self.descriptor = mutant[1]

        # all associated genes
        self.genes = mutant[2].split(" | ")

        # all associated phenotypes
        self.phenotypes = mutant[3].split(" | ")

        self.null = False
        self.overexp = False
        self.multiple = False
        self.develop = False
        self.other = False


class DictyMutants:
    DEFAULT_DATABASE_PATH = serverfiles.localpath(DOMAIN)  # use a default local folder for storing the genesets

    def __init__(self, local_database_path=None):
        """ A collection of Dictybase mutants as a dictionary of :obj:`DictyMutant` objects.

        Args:
            local_database_path: A path for storing D. dictyostelium mutants objects. If `None` then
                                a default database path is used.
        """

        self.local_database_path = local_database_path \
            if local_database_path is not None else self.DEFAULT_DATABASE_PATH

        if not os.path.exists(self.local_database_path):
            os.mkdir(self.local_database_path)

        self._mutants = pickle.load(open(serverfiles.localpath_download(DOMAIN, PHENOTYPES_FILENAME), "rb"))

    def update_file(self, name):
        url = "http://dictybase.org/db/cgi-bin/dictyBase/download/download.pl?area=mutant_phenotypes&ID="
        filename = os.path.join(self.local_database_path, name)
        temp_file = os.path.join(self.local_database_path, name + "_temp")
        stream = urlopen(url + name)

        with open(temp_file, "wb") as file:
            shutil.copyfileobj(stream, file)

        os.rename(temp_file, filename)
        return filename

    def load_mutants(self, file):
        data = open(file)
        data.readline()  # remove data_header
        data = data.read()
        return data.splitlines()

    def download_mutants(self):
        all_mutants = self.load_mutants(self.update_file("all-mutants.txt"))
        null_mutants = self.load_mutants(
            self.update_file("null-mutants.txt"))
        overexp_mutants = self.load_mutants(
            self.update_file("overexpression-mutants.txt"))
        multiple_mutants = self.load_mutants(
            self.update_file("multiple-mutants.txt"))
        develop_mutants = self.load_mutants(
            self.update_file("developmental-mutants.txt"))
        other_mutants = self.load_mutants(
            self.update_file("other-mutants.txt"))

        _mutants = [DictyMutant(mutant) for mutant in all_mutants]

        the_nulls = set([DictyMutant(line).name for line in null_mutants])
        the_overexps = set([DictyMutant(line).name for line in overexp_mutants])
        the_multiples = set([DictyMutant(line).name for line in multiple_mutants])
        the_develops = set([DictyMutant(line).name for line in develop_mutants])
        the_others = set([DictyMutant(line).name for line in other_mutants])

        for mutant in _mutants:
            if mutant.name in the_nulls: mutant.null = True
            if mutant.name in the_overexps: mutant.overexp = True
            if mutant.name in the_multiples: mutant.multiple = True
            if mutant.name in the_develops: mutant.develop = True
            if mutant.name in the_others: mutant.other = True

        final_mutants = {x: x for x in _mutants}
        return final_mutants

    def pickle_data(self):
        return pickle.dumps(self.download_mutants(), -1)

    @classmethod
    def get_instance(cls):
        if not hasattr(cls, "_shared_dict"):
            dicty = DictyMutants()
            cls._shared_dict = dicty.__dict__
        instance = DictyMutants.__new__(DictyMutants)
        instance.__dict__ = cls._shared_dict
        return instance

    def mutants(self):
        return self._mutants.keys()

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
        for mutant, genes in [(mutant, self.mutant_genes(mutant))
                              for mutant in self.mutants()]:
            for gene in genes:
                dgm[gene].add(mutant)
        return dgm

    def phenotype_mutants(self):
        dpm = defaultdict(set)
        for mutant, phenotypes in [(mutant, self.mutant_phenotypes(mutant))
                                   for mutant in self.mutants()]:
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


def download_mutants():
    return DictyMutants.get_instance().pickle_data()


if __name__ == "__main__":
    dicty_mutants = mutants()
    mutant = list(dicty_mutants)[0]
    print(mutant.name, mutant.descriptor, mutant.genes)
    print(gene_mutants())
