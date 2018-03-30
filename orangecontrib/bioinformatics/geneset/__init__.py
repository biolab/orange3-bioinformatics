""" GeneSet module """
from orangecontrib.bioinformatics.geneset.config import *
from orangecontrib.bioinformatics.geneset.utils import (
    filename, GeneSets, GeneSet, NoGeneSetsException, GeneSetException, filename_parse
)
from orangecontrib.bioinformatics.utils import serverfiles


def list_all(**kwargs):
    """ Returns available gene sets from the server files repository: a list of (hierarchy, organism)
    """
    organism = kwargs.get("organism", None)

    all_available = [filename_parse(f_name) for domain, f_name
                     in serverfiles.ServerFiles().listfiles(DOMAIN) + serverfiles.listfiles(DOMAIN)]
    if organism:
        return [(hier, org) for hier, org in all_available if org == organism]
    else:
        return all_available


def load_gene_sets(hierarchy):
    file_path = serverfiles.localpath_download(DOMAIN, filename(*hierarchy))
    return GeneSets.from_gmt_file_format(file_path)
