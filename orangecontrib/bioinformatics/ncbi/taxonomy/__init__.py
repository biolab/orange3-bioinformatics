""" NCBI Taxonomy browser module """
import warnings

from orangecontrib.bioinformatics.utils import serverfiles
from .utils import UnknownSpeciesIdentifier, MultipleSpeciesException, Taxonomy
from .config import *


def common_taxids():
    """ Return taxonomy IDs for most common organisms. """
    return [tax_id for tax_id, _ in COMMON_NAMES]


def common_taxid_to_name(tax_id):
    """ Return a name for a common organism taxonomy id. """
    return COMMON_NAMES_MAPPING[str(tax_id)]


def taxname_to_taxid(name):
    """ Return taxonomy ID for a taxonomy name. """
    name_to_taxid = dict(map(reversed, COMMON_NAMES))

    if name in name_to_taxid:
        return name_to_taxid[name]
    return None


def shortname(taxid):
    """ Short names for common_taxids organisms """
    if taxid in SHORT_NAMES:
        return SHORT_NAMES[taxid]
    return []


def name(tax_id):
    """ Return the scientific name for organism with provided taxonomy id.
    """
    # Most of the lookups will be for the common names, so in most
    # situations we can avoid loading the taxonomy.
    if tax_id in COMMON_NAMES_MAPPING:
        return COMMON_NAMES_MAPPING[tax_id]
    else:
        return Taxonomy()[tax_id]


def other_names(tax_id):
    """ Return a list of (name, name_type) tuples excluding the scientific name.
    Use :func:`name` to retrieve the scientific name.
    """
    return Taxonomy().other_names(tax_id)


def search(string, onlySpecies=True, exact=False):
    """ Search the NCBI taxonomy database for an organism.

    :param string: Search string.
    :param onlySpecies: Return only taxids of species (and subspecies).
    :param exact:  Return only taxids of organism that exactly match the string.
    """
    ids = Taxonomy().search(string, onlySpecies, exact)
    return list(ids)


def lineage(taxid):
    """ Return a list of taxids ordered from the topmost node (root) to taxid.
    """
    return Taxonomy().lineage(taxid)


def to_taxid(code, mapTo=None):
    """ See if the code is a valid code in GO or KEGG and return a set of its taxids.
    """
    warnings.warn("'to_taxid' is deprecated", DeprecationWarning)

    try:
        name(code)
        results = set([code])
    except UnknownSpeciesIdentifier:
        results = set()
        from orangecontrib.bio import kegg, go
        for test in [kegg.to_taxid, go.to_taxid]:
            try:
                r = test(code)
                if type(r) == set:
                    results.update(r)
                elif r:
                    results.add(r)
            except Exception:
                pass

    if mapTo:
        mapTo = set(mapTo)
        if mapTo.intersection(results):
            return mapTo.intersection(results)

        mapped = set()
        for r in results:
            r_lin = lineage(r)
            if mapTo.intersection(r_lin):
                mapped.update(mapTo.intersection(r_lin))

        if not mapped:
            for m in mapTo:
                m_lin = lineage(m)
                if results.intersection(m_lin):
                    mapped.add(m)
        results = mapped

    return results


def ensure_downloaded(callback=None, verbose=True):
    """ Retrieve the taxonomy database if not already downloaded.
    """
    serverfiles.localpath_download(DOMAIN, FILENAME, callback=callback, verbose=verbose)
