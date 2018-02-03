""" GeneSets utility functions """
import os
import re
import pickle
import tempfile

from orangecontrib.bioinformatics.geneset.config import DOMAIN
from orangecontrib.bioinformatics.utils import environ, serverfiles


def filename(hierarchy, organism):
    """ Obtain a filename for given hierarchy and organism.
    """
    return "gs_" + "_._".join(hierarchy + (organism if organism is not None else '',)) + '.pck'


def pickle_temp(obj):
    """ Pickle a file to a temporary file and returns its name
    """
    file_descriptor, temp_file_name = tempfile.mkstemp()
    os.close(file_descriptor)

    with open(temp_file_name, 'wb') as f:
        pickle.dump(obj, f)
        return temp_file_name


def local_path():
    """ Returns local path for gene sets. Creates it if it does not exists yet.
    """
    return os.makedirs(os.path.join(environ.buffer_dir, DOMAIN), exist_ok=True)


def filename_parse(fn):
    """ Returns a hierarchy and the organism from the gene set filename format.
    """
    fn = fn[3:-4]
    parts = fn.split('_._')
    hierarchy = tuple(parts[:-1])
    org = parts[-1] if parts[-1] != '' else None
    return hierarchy, org


def list_gene_sets():
    """ Returns available gene sets from the server files repository: a list of (hierarchy, organism, on_local)
    """

    return [filename_parse(f_name) + (False,) for domain, f_name in serverfiles.ServerFiles().listfiles(DOMAIN)]


def only_option(a):
    if len(a) == 1:
        return list(a)[0]
    else:
        raise Exception()


def is_sequencens(x):
    """ Is x a sequence and not string ?
    We say it is if it has a __getitem__ method and it is not an instance of string.
    """
    return hasattr(x, '__getitem__') and not isinstance(x, str)


def str_or_none(x):
    return str(x) if x is not None else x


def gmt_file_loader(contents, name):
    """ Each line consists of tab separated elements. First is
    the geneset name, next is it's description.

    For now the description is skipped.

    Example Gene Set (.gmt) file:
        anti-liver_sw   anti-liver_sw   BMPR1A  APOD    WSB1    BMI1    SLC2A1  ...
        B-cells_sw  B-cells_sw  E2F5    NCF1    PALM2-AKAP2 IRF4    SLC2A1  ...
        Bladder_sw  Bladder_sw  PLCD4   ANGPTL1 LOC286191   ST0N1   LOC283904   ...
        cerebellum_sw   cerebellum_sw   C19orf43    LOC653464   KI110802    ...
        Cervix_sw   Cervix_sw   LAMA4   GSTM5   SNX19   DKK1    NT5E    ...
    """
    linkre = re.compile("(.*?)\s?(?:\[(https?://[^[^\]]*)\])?$")

    def _hline(line):
        tabs = [tab.strip() for tab in line.split("\t")]
        groups = linkre.match(tabs[1]).groups()
        return GeneSet(id=tabs[0],
                       description=groups[0],
                       link=groups[1],
                       hierarchy=('Custom', name),
                       genes=tabs[2:])

    def _handle_ne_lines(s, fn):
        """ Run function on nonempty lines of a string.
        Return a list of results for each line.
        """
        lines = (l.strip() for l in s.splitlines())
        return [fn(l) for l in lines if l]

    return GeneSets(_handle_ne_lines(contents, _hline))


class GeneSet:
    """ A single set of genes.
    """

    def __init__(self, genes=[], name=None, id=None, description=None, link=None,
                 organism=None, hierarchy=None, pair=None):
        """
        :param pair: Only for backward compatibility: convert a tuple (name, genes)
            into this object.
        """

        self.hierarchy = hierarchy
        """ Hierarchy should be formated as a tuple, for example ``("GO", "biological_process")``"""

        self.genes = set(genes)
        """ A set of genes. Genes are strings. """

        self.name = name
        """ Gene set name. """

        self.id = id
        """ Short gene set ID. """

        self.description = description
        """ Gene set description. """

        self.link = link
        """ Link to further information about this gene set. """

        self.organism = organism
        """ Organism as a NCBI taxonomy ID. """

        if pair:
            self.id, self.genes = pair[0], set(pair[1])

    """ The following functions are needed for sets of gene sets to be able to assess equality
    """

    def __hash__(self):
        return self.id.__hash__() + self.name.__hash__()

    def __eq__(self, other):
        if isinstance(other, self.__class__):
            return self.__dict__ == other.__dict__
        else:
            return False

    def __ne__(self, other):
        return not self.__eq__(other)

    def size(self):
        return len(self.genes)

    def cname(self, source=True, name=True):
        """ Return a gene set name with hieararchy. """
        oname = self.id
        if source and self.hierarchy:
            oname = "[ " + ", ".join(self.hierarchy) + " ] " + oname
        if name and self.name:
            oname = oname + " " + self.name
        return oname

    def to_odict(self, source=True, name=True):
        """ For backward compatibility. Return a gene set as a tuple
        (id, list of genes).
        """
        return self.cname(source=source, name=name), self.genes

    def __repr__(self):
        return "GeneSet(" + ", ".join([
            "id=" + str(self.id),
            "genes=" + str(self.genes),
            "name=" + str(self.name),
            "link=" + str(self.link),
            "hierarchy=" + str(self.hierarchy)]) + ")"


class GeneSets(set):
    """ A collection of gene sets: contains :class:`GeneSet` objects.
    """

    def __init__(self, input=None):
        """ If `input` is a dictionary, the gene sets are converted to the current format.
        """
        if input is not None and len(input) > 0:
            self.update(input)

    def update(self, input):
        if input.__class__.__name__ == "GeneSets":  # HACK: because location can appear different
            super(GeneSets, self).update(input)
        else:
            prepared_genesets = []  # parse them all before adding, so that it fails on error
            if hasattr(input, "items"):
                for i, g in input.items():
                    prepared_genesets.append(GeneSet(pair=(i, g)))
            else:
                for i in input:
                    if isinstance(i, GeneSet):
                        prepared_genesets.append(i)
                    else:
                        i, g = i
                        prepared_genesets.append(GeneSet(pair=(i, g)))

            for g in prepared_genesets:
                self.add(g)

    def to_odict(self):
        """ Return gene sets in old dictionary format. """
        return dict(gs.to_odict() for gs in self)

    def set_hierarchy(self, hierarchy):
        """ Sets hierarchy for all gene sets. """
        for gs in self:
            gs.hierarchy = hierarchy

    def __repr__(self):
        return "GeneSets(" + set.__repr__(self) + ")"

    def common_org(self):
        """ Return a common organism. """
        if len(self) == 0:
            raise GeneSetRegException("Empty gene sets.")

        organisms = set(a.organism for a in self)
        try:
            return only_option(organisms)
        except:
            raise GeneSetRegException("multiple organisms: " + str(organisms))

    def hierarchies(self):
        """ Return all hierarchies. """
        if len(self) == 0:
            raise GeneSetRegException("Empty gene sets.")
        return set(a.hierarchy for a in self)

    def common_hierarchy(self):
        """ Return a common hierarchy. """
        hierarchies = self.hierarchies()

        def common_hierarchy1(hierarchies):
            def hier(l): return set(map(lambda x: x[:currentl], hierarchies))

            currentl = max(map(len, hierarchies))
            while len(hier(currentl)) > 1:
                currentl -= 1
            return only_option(hier(currentl))

        return common_hierarchy1(hierarchies)

    def split_by_hierarchy(self):
        """ Split gene sets by hierarchies. Return a list of :class:`GeneSets` objects. """
        hd = dict((h, GeneSets()) for h in self.hierarchies())
        for gs in self:
            hd[gs.hierarchy].add(gs)
        return hd.values()


class GMTFileFormatException(Exception):
    pass


class NoGeneSetsException(Exception):
    """ Raised when provided taxonomy is not in orangecontrib.bio.ncbi.taxonomy.common_taxids """


class GeneSetRegException(Exception):
    pass


class GeneSetIDException(Exception):
    pass


if __name__ == '__main__':
    print(list_gene_sets())
