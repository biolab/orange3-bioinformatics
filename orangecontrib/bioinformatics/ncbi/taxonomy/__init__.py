""" NCBI Taxonomy browser module """

from orangecontrib.bioinformatics.ncbi.taxonomy.utils import Taxonomy, UnknownSpeciesIdentifier

COMMON_NAMES = (
    ("6500", "Aplysia californica"),
    ("3702", "Arabidopsis thaliana"),
    ("9913", "Bos taurus"),
    ("6239", "Caenorhabditis elegans"),
    ("5476", "Candida albicans"),
    ("3055", "Chlamydomonas reinhardtii"),
    ("7955", "Danio rerio"),
    ("44689", "Dictyostelium discoideum"),
    ("7227", "Drosophila melanogaster"),
    ("562", "Escherichia coli"),
    ("11103", "Hepatitis C virus"),
    ("9606", "Homo sapiens"),
    ("10090", "Mus musculus"),
    ("2104", "Mycoplasma pneumoniae"),
    ("4530", "Oryza sativa"),
    ("5833", "Plasmodium falciparum"),
    ("4754", "Pneumocystis carinii"),
    ("10116", "Rattus norvegicus"),
    ("4932", "Saccharomyces cerevisiae"),
    ("4896", "Schizosaccharomyces pombe"),
    ("31033", "Takifugu rubripes"),
    ("8355", "Xenopus laevis"),
    ("4577", "Zea mays"),
)

SHORT_NAMES = {
    "6500": ["aplysia"],
    "3702": ["arabidopsis", "thaliana", "plant"],
    "9913": ["cattle", "cow"],
    "6239": ["nematode", "roundworm"],
    "5476": ["thrush", "candidiasis", "candida"],
    "3055": ["algae"],
    "7955": ["zebrafish"],
    "44689": ["dicty", "amoeba", "slime mold"],
    "7227": ["fly", "fruit fly", "vinegar fly"],
    "562": ["ecoli", "coli", "bacterium"],
    "11103": ["virus, hepatitis"],
    "9606": ["human"],
    "10090": ["mouse", "mus"],
    "2104": ["bacterium", "mycoplasma"],
    "4530": ["asian rice", "rice", "cereal", "plant"],
    "5833": ["plasmodium", "malaria", "parasite"],
    "4754": ["pneumonia", "fungus"],
    "10116": ["rat", "laboratory rat"],
    "4932": ["yeast", "baker yeast", "brewer yeast"],
    "4896": ["yeast", "fission yeast"],
    "31033": ["fish", "pufferfish"],
    "8355": ["frog", "african clawed frog"],
    "4577": ["corn", "cereal grain", "plant"],
}


COMMON_NAMES_MAPPING = dict(COMMON_NAMES)


def common_taxids():
    """ Return taxonomy IDs for most common organisms.

    :rtype: :class:`list` of common organisms
    """
    return [tax_id for tax_id, _ in COMMON_NAMES]


def common_taxid_to_name(tax_id):
    """ Return a name for a common organism taxonomy id.

    :param tax_id: Taxonomy if (NCBI taxonomy database)
    :type tax_id: str

    :rtype: :class:`str`
    """
    if str(tax_id) == '352472':
        tax_id = '44689'
    return COMMON_NAMES_MAPPING[str(tax_id)]


def species_name_to_taxid(organism_name):
    """ Return taxonomy ID for a taxonomy name.

    :param organism_name: Official organism name, e.g., 'Homo sapiens'
    :type organism_name: str

    :rtype: :class:`str`
    """
    name_to_taxid = dict(map(reversed, COMMON_NAMES))

    if organism_name in name_to_taxid:
        return name_to_taxid[organism_name]
    return None


def shortname(tax_id):
    """ Short names for common_taxids organisms

    :param tax_id: Taxonomy if (NCBI taxonomy database)
    :type tax_id: str

    """
    if tax_id in SHORT_NAMES:
        return SHORT_NAMES[tax_id]
    return []


def name(tax_id):
    """ Return the scientific name for organism with provided taxonomy id.

    :param tax_id: Taxonomy if (NCBI taxonomy database)
    :type tax_id: str

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

    :param tax_id: Taxonomy if (NCBI taxonomy database)
    :type tax_id: str

    """
    return Taxonomy().other_names(tax_id)


def search(string, only_species=True, exact=False):
    """ Search the NCBI taxonomy database for an organism.

    :param string: Search string.
    :type string: str

    :param only_species: Return only taxids of species (and subspecies).
    :type only_species: bool

    :param exact:  Return only taxids of organism that exactly match the string.
    :type exact: bool
    """
    ids = Taxonomy().search(string, only_species, exact)
    return list(ids)


def lineage(tax_id):
    """ Return a list of taxids ordered from the topmost node (root) to taxid.

    :param tax_id: Taxonomy if (NCBI taxonomy database)
    :type tax_id: str
    """

    return Taxonomy().lineage(tax_id)
