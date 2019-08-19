""" GeneSet module """
from typing import Tuple

from orangecontrib.bioinformatics.utils import serverfiles
from orangecontrib.bioinformatics.geneset.utils import (
    GeneSet,
    GeneSets,
    GeneSetException,
    NoGeneSetsException,
    filename,
    filename_parse,
)

__all__ = (GeneSet, GeneSets, GeneSetException, NoGeneSetsException)
DOMAIN = 'gene_sets'


def list_all(**kwargs):
    """ Returns available gene sets from the server files repository.

    :param kwargs:
        * *organism* (``str``) -- Taxonomy id (NCBI taxonomy database)


    :rtype: :obj:`list` of (hierarchy, organism)


    Example
    --------
    The available gene set collection can be listed with:
        >>> list_all(organism='10090')

    """

    organism = kwargs.get("organism", None)

    all_available = {
        filename_parse(f_name)
        for domain, f_name in serverfiles.ServerFiles().listfiles(DOMAIN) + serverfiles.listfiles(DOMAIN)
    }

    if organism:
        return [hier for hier, org in all_available if org == organism]
    else:
        return all_available


def load_gene_sets(hierarchy, tax_id):
    # type: (Tuple[Tuple(str, str), str]) -> GeneSets
    """ Initialize gene sets from a given hierarchy.

    :param tuple hierarchy: gene set hierarchy.
    :rtype: :obj:`GeneSets`

    Example
    --------
    Gene sets provided with Orange are organized hierarchically:
        >>> list_of_genesets= list_all(organism='10090')
            [(('KEGG', 'Pathways'), '10090'),
             (('KEGG', 'pathways'), '10090'),
             (('GO', 'biological_process'), '10090'),
             (('GO', 'molecular_function'), '10090'),
             (('GO', 'cellular_component'), '10090')]
        >>> load_gene_sets(list_of_genesets[0])

    """
    file_path = serverfiles.localpath_download(DOMAIN, filename(hierarchy, tax_id))
    return GeneSets.from_gmt_file_format(file_path)
