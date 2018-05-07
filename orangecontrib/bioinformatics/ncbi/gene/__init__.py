""" NCBI GeneInformation module """
import pickle

from collections import defaultdict
from typing import List

from orangecontrib.bioinformatics.ncbi.gene.config import *
from orangecontrib.bioinformatics.ncbi.gene.utils import GeneInfoDB, parse_sources, parse_synonyms
from orangecontrib.bioinformatics.utils import serverfiles, ensure_type


_no_hits, _single_hit, _multiple_hits = 0, 1, 2
_source, _symbol, _synonym, _locus, _gene_id, _nom_symbol = \
    'External reference', 'Symbol', 'Synonym', 'Locus tag', 'NCBI ID', 'Nomenclature symbol'


class NoGeneNamesException(Exception):
    """ Can not extract gene names from the table. Table can't have both row and column gene variables """


class Gene:
    """ Basic representation of a gene. Class holds NCBI and GeneMatcher data created upon class initialization.
    GeneMatcher tags are used only for matching purposes.

    We deny __dict__ creation by defining __slots__ for memory optimization. We expect a lot of instances
    of this class to be created.
    """
    __slots__ = GENE_INFO_TAGS + GENE_MATCHER_TAGS  # Note: slots are ordered.

    def __init__(self, input_name=None):
        self.input_name = input_name
        self.type_of_match = None
        self.ncbi_id = None
        self._possible_hits = []

    @property
    def possible_hits(self):
        return self._possible_hits

    @possible_hits.setter
    def possible_hits(self, genes):
        for gene in genes:
            possible_match = Gene()
            possible_match.ncbi_id = gene
            self._possible_hits.append(possible_match)

    def load_ncbi_info(self):
        if not self.ncbi_id:
            return

        info = GeneInfoDB().select_gene_info(self.ncbi_id)
        for attr, value in zip(self.__slots__, info):
            if attr == 'db_refs':
                value = parse_sources(value)
            elif attr == 'synonyms':
                value = parse_synonyms(value)

            setattr(self, attr, value)


class GeneInfo(dict):

    def __init__(self, organism):
        super().__init__()

        self._init_gene_info(organism)

    def _init_gene_info(self, organism):
        for gene in GeneInfoDB().select_genes_by_organism(organism):
            gene_obj = Gene()

            for attr, value in zip(gene_obj.__slots__, gene):
                if attr == 'db_refs':
                    value = parse_sources(value)
                elif attr == 'synonyms':
                    value = parse_synonyms(value)

                setattr(gene_obj, attr, value)

            self[gene_obj.gene_id] = gene_obj

    def get_gene_by_id(self, gene_id):
        """ Search and return the Gene object for gene_id
        """
        try:
            return self[int(gene_id)]
        except KeyError as e:
            # TODO: handle this properly?
            # raise e
            pass


class GeneMatcher:

    def __init__(self,  tax_id, **kwargs):
        self._organism = ensure_type(tax_id, str)  # type: str
        self._genes = []                           # type: (List[str])

        self._case_insensitive = kwargs.get("case_insensitive", False)
        self._matcher = self.load_matcher_file(DOMAIN, MATCHER_FILENAME.format(tax_id))

    @property
    def organism(self):
        return self.organism

    @organism.setter
    def organism(self, tax_id):
        self._organism = ensure_type(tax_id, str)
        self._matcher = self.load_matcher_file(DOMAIN, MATCHER_FILENAME.format(tax_id))

    @property
    def genes(self):
        return self._genes

    @genes.setter
    def genes(self, gene_data):  # type: (List[str]) -> None
        self._genes = [Gene(input_name=gene) for gene in gene_data]

    def get_known_genes(self):
        return [gene for gene in self.genes if gene.ncbi_id]

    def map_input_to_ncbi(self):
        return {gene.input_name: gene.ncbi_id for gene in self.genes if gene.ncbi_id}

    def match_table_attributes(self, data_table):
        input_gene_names = [var.name for var in data_table.domain.attributes]

        if input_gene_names:
            self.genes = input_gene_names
            self.run_matcher()

            for gene in self.genes:
                if gene.ncbi_id:
                    data_table.domain[gene.input_name].attributes[NCBI_ID] = gene.ncbi_id

    def run_matcher(self, progress_callback=None):
        """ This will try to match genes, with ncbi ids, based on provided input of genes.


        Args:
            progress_callback: Used for progress bar in widgets. It emits the signal back to main thread

        """
        if progress_callback:
            self._match(callback=progress_callback.emit)
        else:
            self._match()

    def _match(self, **kwargs):
        callback = kwargs.get("callback", None)

        def match_input(mapper, input_name):
            return mapper[input_name]

        for gene in self.genes:
            if callback:
                callback()

            try:
                ncbi_match = match_input(self._matcher[MAP_GENE_ID], int(gene.input_name))
                if ncbi_match:
                    gene.ncbi_id = ncbi_match[0][MAP_GENE_ID]
                    gene.type_of_match = _gene_id
                    continue
            except ValueError:
                # NCBI ids are stored as Integers. If ValueError is raised, probably not NCBI ID.
                # We expect unique match here.
                pass

            source_match = match_input(self._matcher[MAP_SOURCES], gene.input_name)
            # ids from different sources are unique. We do not expect to get multiple hits here.
            # There is exceptions with organism 3702. It has same source id from Araport or TAIR databases.
            if source_match:
                gene.ncbi_id = source_match[0][MAP_GENE_ID]
                gene.type_of_match = _source
                continue

            symbol_match = match_input(self._matcher[MAP_SYMBOL], gene.input_name)
            if len(symbol_match) == _single_hit:
                gene.ncbi_id = symbol_match[0][MAP_GENE_ID]
                gene.type_of_match = _symbol
                continue
            elif len(symbol_match) >= _multiple_hits:
                gene.possible_hits = symbol_match

            locus_match = match_input(self._matcher[MAP_LOCUS], gene.input_name)
            if len(locus_match) == _single_hit:
                gene.ncbi_id = locus_match[0][MAP_GENE_ID]
                gene.type_of_match = _locus
                continue
            elif len(symbol_match) >= _multiple_hits:
                gene.possible_hits = locus_match

            synonym_match = match_input(self._matcher[MAP_SYNONYMS], gene.input_name)
            if len(synonym_match) == _single_hit:
                gene.ncbi_id = synonym_match[0][MAP_GENE_ID]
                gene.type_of_match = _synonym
                continue
            elif len(synonym_match) >= _multiple_hits:
                gene.possible_hits = synonym_match

            nomenclature_match = match_input(self._matcher[MAP_NOMENCLATURE], gene.input_name)
            if len(nomenclature_match) == _single_hit:
                gene.ncbi_id = nomenclature_match[0][MAP_GENE_ID]
                gene.type_of_match = _nom_symbol
                continue
            elif len(nomenclature_match) >= _multiple_hits:
                gene.possible_hits = nomenclature_match

    def load_matcher_file(self, domain, filename):
        # this starts download if files are not on local machine
        file_path = serverfiles.localpath_download(domain, filename)
        # download new version before using this file for gene name matching
        serverfiles.update(domain, filename)

        def case_insensitive_keys(matcher_dict):
            updated_dict = {MAP_SOURCES:  matcher_dict[MAP_SOURCES],
                            MAP_GENE_ID: matcher_dict[MAP_GENE_ID],
                            MAP_LOCUS:    matcher_dict[MAP_LOCUS],
                            MAP_SYNONYMS: defaultdict(list),
                            MAP_SYMBOL:  defaultdict(list),
                            MAP_NOMENCLATURE: defaultdict(list)}

            for key, value in matcher_dict[MAP_SYMBOL].items():
                # ensure string, we are using string methods (upper, lower)
                key = ensure_type(str(key), str)
                updated_dict[MAP_SYMBOL][key] = value
                updated_dict[MAP_SYMBOL][key.upper()] = value
                updated_dict[MAP_SYMBOL][key.lower()] = value

            for key, value in matcher_dict[MAP_SYNONYMS].items():
                key = ensure_type(str(key), str)
                updated_dict[MAP_SYNONYMS][key] = value
                updated_dict[MAP_SYNONYMS][key.upper()] = value
                updated_dict[MAP_SYNONYMS][key.lower()] = value

            for key, value in matcher_dict[MAP_NOMENCLATURE].items():
                key = ensure_type(str(key), str)
                updated_dict[MAP_NOMENCLATURE][key] = value
                updated_dict[MAP_NOMENCLATURE][key.upper()] = value
                updated_dict[MAP_NOMENCLATURE][key.lower()] = value

            return updated_dict

        with open(file_path, 'rb') as pickle_file:
            if not self._case_insensitive:
                return pickle.load(pickle_file)
            else:
                return case_insensitive_keys(pickle.load(pickle_file))


if __name__ == "__main__":
    g = Gene()
    g.ncbi_id = 1
    g.load_ncbi_info()
    print(g.__slots__)
    print(g.ncbi_id, g.description, g.db_refs, g.synonyms)
