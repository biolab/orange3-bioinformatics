""" GeneSet module """
import pickle


from pathlib import Path
from collections import defaultdict


from orangecontrib.bioinformatics.geneset.config import *
from orangecontrib.bioinformatics.geneset.utils import (
    pickle_temp, filename, local_path, list_gene_sets,
    GeneSets, GeneSet, is_sequencens, gmt_file_loader, str_or_none,
    GMTFileFormatException, NoGeneSetsException, GeneSetRegException)

from orangecontrib.bioinformatics.utils import serverfiles
from orangecontrib.bioinformatics.ncbi import taxonomy


def list_all(organism=None, local=None):
    """ Return gene sets available in the local and ServerFiles repositories.
    It returns a list of tuples of (hierarchy, organism, available_locally)

    Results can be filtered with the following parameters.

    :param str organism: Organism tax id.
    :param bool local: Available locally.
    """
    d = {}
    for hier, org, is_local in list_gene_sets():
        d[hier, org] = min(is_local, d.get((hier, org), True))

    return [(hier, org, is_local) for (hier, org), is_local in d.items()
            if (local is None or is_local == local) and (organism is None or org == str(organism))]


def load_gene_sets(hierarchy):
    file_path = serverfiles.localpath_download(DOMAIN, filename(*hierarchy))
    with open(file_path, 'rb') as file:
        return pickle.load(file, encoding="latin1")


def collections(*args):
    """ Load gene sets from various sources: GMT file, GO, KEGG, and others.
    Return an instance of :class:`GeneSets`.

    Each arguments specifies a gene set and can be either:

    * a filename of a GMT file,
    * a tuple (hierarchy, organism) (for example ``(("KEGG",), "10090")``), or
    * an instance of :class:`GeneSets`
    """

    result = GeneSets()

    for collection in args:
        try:
            result.update(collection)
        except (ValueError, TypeError) as e:
            if is_sequencens(collection):  # have a hierarchy, organism specification
                new = load(*collection)
                result.update(new)
            else:
                if Path(collection).suffix == '.gmt':  # supported file format
                    with open(collection, 'rt') as col:
                        result.update(gmt_file_loader(col.read(), collection))
                else:
                    raise GMTFileFormatException('collection() accepts files in .gmt format only.')

    return result


def load(hierarchy, organism):
    """ load gene_set file from given hierarchy and organism.
    If the file is not available, load it from the server files.
    """
    if organism is not None and organism not in taxonomy.common_taxids():
            raise NoGeneSetsException('Unknown organism id {}. '
                                      'Orange supports {}'.format(str(organism), str(taxonomy.common_taxids())))

    def _build_hierarchy_dict(files):
        hierd = defaultdict(list)
        for ind, f in enumerate(files):
            hier, org = f
            for i in range(len(hier) + 1):
                hierd[(hier[:i], org)].append(ind)
        return hierd

    def _load_gene_sets(hierarchy, organism):
        files = list(map(lambda x: x[:2], list_gene_sets()))
        hierd = _build_hierarchy_dict(files)
        out = GeneSets()
        matches = hierd[(hierarchy, organism)]
        if not matches:
            raise NoGeneSetsException('No gene sets for {} (org {} )'.format(str(hierarchy), str(organism)))

        for (h, o) in [files[i] for i in hierd[(hierarchy, organism)]]:
            fname = serverfiles.localpath_download(DOMAIN, filename(h, o))
            out.update(pickle.load(open(fname, 'rb'), encoding="latin1"))

        return out

    return _load_gene_sets(hierarchy, str_or_none(organism))
