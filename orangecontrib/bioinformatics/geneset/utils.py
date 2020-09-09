""" GeneSets utility functions """
from typing import List, Tuple, NamedTuple

import numpy as np

from orangecontrib.bioinformatics.utils import ensure_type
from orangecontrib.bioinformatics.utils.statistics import Hypergeometric


def filename(hierarchy, organism):  # type: (Tuple[str, str], str) -> str
    """Obtain a filename for given hierarchy and organism.

    :param hierarchy: GeneSet hierarchy, example: ('GO', 'biological_process')
    :param organism: Taxonomy ID

    :return: Filename for given hierarchy and organism

    Example
    --------
        >>> filename(('CustomSet', 'subsets'), '6500')
        'CustomSet-subsets-6500.gmt'
    """

    return '-'.join(hierarchy + (organism if organism is not None else '',)) + '.gmt'


def filename_parse(fn):  # type: (str) -> (Tuple[Tuple[str, str], str])
    """Returns a hierarchy and the organism from the gene set filename format.

    :param fn: GeneSets file name (.gmt)

    :return: A hierarchy and taxonomy id for given filename

    Example
    --------
        >>> filename_parse('Custom-set-6500.gmt')
        (('Custom', 'set'), '6500')
    """

    file_name = fn[:-4]  # removes .gmt extension
    parts = file_name.split('-')
    hierarchy = tuple(parts[:-1])
    org = parts[-1] if parts[-1] != '' else None
    return hierarchy, org


GENE_SET_ATTRIBUTES = ('gs_id', 'hierarchy', 'organism', 'name', 'genes', 'description', 'link')
HYPERGEOMETRIC = Hypergeometric()

# change this when python 3.4 is not supported anymore
enrichment_result = NamedTuple(
    'enrichment_result', [('query', set), ('reference', set), ('p_value', float), ('enrichment_score', float)]
)


class GeneSet:
    __slots__ = GENE_SET_ATTRIBUTES

    def __init__(self, gs_id=None, hierarchy=None, organism=None, name=None, genes=None, description=None, link=None):
        """Object representing a single set of genes

        :param gs_id: Short gene set ID.
        :param hierarchy: Hierarchy should be formated as a tuple, for example ``("GO", "biological_process")``
        :param organism: Organism as a NCBI taxonomy ID.
        :param name: Gene set name.
        :param genes: A set of genes. Genes are strings.
        :param description: Gene set description.
        :param link: Link to further information about this gene set.
        """

        self.gs_id = gs_id
        self.hierarchy = hierarchy
        self.organism = organism
        self.name = name
        self.genes = genes
        self.description = description
        self.link = link

    def __hash__(self):
        return self.gs_id.__hash__() + self.name.__hash__()

    def __eq__(self, other):
        if isinstance(other, self.__class__):
            if self.__slots__ == other.__slots__:
                return all(getattr(self, attr) == getattr(other, attr) for attr in self.__slots__)

        return False

    def __lt__(self, other):
        return (self.name + self.gs_id) < (other.name + other.gs_id)

    def set_enrichment(self, reference, query):  # type: (List, List) -> enrichment_result
        """
        Args:
            reference:
            query:
        """

        assert len(reference) > 0
        query_mapped = self.genes.intersection(query)
        reference_mapped = self.genes.intersection(reference)

        query_p = len(query_mapped) / len(query) if query else np.nan
        ref_p = len(reference_mapped) / len(reference) if reference else np.nan
        enrichment = query_p / ref_p if ref_p else np.nan

        return enrichment_result(
            set(query_mapped),
            set(reference_mapped),
            HYPERGEOMETRIC.p_value(len(query_mapped), len(reference), len(reference_mapped), len(query)),
            enrichment,
        )

    def gmt_description(self):
        """Represent GeneSet as line in GMT file format

        :return: Comma-separated GeneSet attributes.
        """

        # this is just so that we preserve order from GENE_SET_ATTRIBUTES for easy access
        empty_field = '_'
        hierarchy = '-'.join((hier for hier in self.hierarchy)) if self.hierarchy else empty_field

        return ','.join(
            (
                self.gs_id,
                hierarchy,
                self.organism,
                self.name if self.name else empty_field,
                empty_field,
                self.description if self.description else empty_field,
                self.link if self.link else empty_field,
            )
        )


class GeneSets(set):
    """A collection of gene sets: contains :obj:`GeneSet` objects."""

    def __init__(self, sets=None):
        # type: (List[GeneSet]) -> None
        super().__init__()

        if sets:
            self.update(ensure_type(sets, list))

    def update(self, sets):
        # type: (List[GeneSet]) -> None

        for g_set in sets:
            self.add(ensure_type(g_set, GeneSet))

    def common_org(self):
        """ Return a common organism. """
        if len(self) == 0:
            raise GeneSetException("Empty gene sets.")

        organisms = {g_set.organism for g_set in self}
        if not len(organisms) == 1:
            raise GeneSetException("multiple organisms in a set " + str(organisms))
        else:
            for org in organisms:
                return org

    def hierarchies(self):
        """ Return all hierarchies. """
        return {g_set.hierarchy for g_set in self}

    def common_hierarchy(self):
        """ Return a common hierarchy. """
        hierarchies = self.hierarchies()

        if not len(hierarchies) == 1:
            raise GeneSetException('multiple hierarchies in a set', str(hierarchies))
        else:
            for org in hierarchies:
                return org

    def delete_sets_by_hierarchy(self, hier):
        selected_sets = self.map_hierarchy_to_sets().get(hier, None)
        if selected_sets:
            for gene_set in selected_sets:
                self.remove(gene_set)

    def map_hierarchy_to_sets(self):
        try:
            split_by_hier = {hier: GeneSets() for hier in self.hierarchies()}
            for gs in self:
                split_by_hier[gs.hierarchy].update([gs])
            return split_by_hier

        except GeneSetException:
            return {}

    def split_by_hierarchy(self):
        """ Split gene sets by hierarchies. Return a list of :class:`GeneSets` objects. """

        try:
            split_by_hier = {hier: GeneSets() for hier in self.hierarchies()}
            [split_by_hier[gs.hierarchy].update([gs]) for gs in self]
            return list(split_by_hier.values())

        except GeneSetException:
            return []

    def genes(self):
        """
        Returns:
            All genes from GeneSets
        """
        genes = set()
        [genes.update(gene_set.genes) for gene_set in self]
        return genes

    def to_gmt_file_format(self, file_path):  # type: (str) -> None
        """The GMT file format is a tab delimited file format that describes gene sets.

        In the GMT format, each row represents a gene set.
        Columns: gs_id    gmt_description    Gene    Gene    Gene    ...
        gmt_description: 'gs_id','hierarchy','organism','name','genes','description','link'

        :param file_path: Path to where file will be created

        """

        with open(file_path, 'w') as gmt_file:
            for gene_set in sorted(self):
                genes = sorted(map(str, gene_set.genes), key=lambda x: int(x))
                line = '\t'.join([gene_set.gs_id, gene_set.gmt_description()] + genes)
                gmt_file.write(line + '\n')

    @staticmethod
    def from_gmt_file_format(file_path):  # type: (str) -> GeneSets
        """Load GeneSets object from GMT file.

        :param file_path: path to a file on local disk
        :rtype: :obj:`GeneSets`
        """
        index = {label: index for index, label in enumerate(GENE_SET_ATTRIBUTES)}

        with open(file_path, 'r', encoding='utf-8') as gmt_file:
            gene_sets = []

            for line in gmt_file:
                columns = [column.strip() for column in line.split('\t')]
                gs_info = columns[1].split(',')
                hierarchy = tuple(gs_info[index['hierarchy']].split('-'))
                genes = {str(gene) for gene in columns[2:]}

                gene_set = GeneSet(
                    gs_id=columns[0],
                    genes=genes,
                    hierarchy=hierarchy,
                    name=gs_info[index['name']],
                    organism=gs_info[index['organism']],
                    description=gs_info[index['description']],
                    link=gs_info[index['link']],
                )

                gene_sets.append(gene_set)

            return GeneSets(gene_sets)


class NoGeneSetsException(Exception):
    """ Raised when provided taxonomy is not in orangecontrib.bio.ncbi.taxonomy.common_taxids """


class GeneSetException(Exception):
    pass


if __name__ == '__main__':
    pass
