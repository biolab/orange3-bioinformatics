""" NCBI GeneInformation module """
import pickle

from collections import defaultdict
from functools import lru_cache
from typing import List, Union

from Orange.data import StringVariable, DiscreteVariable, Domain, Table


from orangecontrib.bioinformatics.ncbi.gene.config import *
from orangecontrib.bioinformatics.ncbi.gene.utils import GeneInfoDB, parse_sources, parse_synonyms
from orangecontrib.bioinformatics.utils import serverfiles, ensure_type
from orangecontrib.bioinformatics.widgets.utils.data import TAX_ID, GENE_ID_COLUMN, GENE_AS_ATTRIBUTE_NAME

_no_hits, _single_hit, _multiple_hits = 0, 1, 2
_source, _symbol, _synonym, _locus, _gene_id, _nom_symbol = \
    'External reference', 'Symbol', 'Synonym', 'Locus tag', 'NCBI ID', 'Nomenclature symbol'


class NoGeneNamesException(Exception):
    """ Can not extract gene names from the table. Table can't have both row and column gene variables """


class Gene:
    """ Basic representation of a gene. Class holds NCBI and GeneMatcher data created upon class initialization.

    GeneMatcher tags are used only for matching purposes.
    """

    # We deny __dict__ creation by defining __slots__ for memory optimization.
    # We expect a lot of instances of this class to be created.
    __slots__ = GENE_INFO_TAGS + GENE_MATCHER_TAGS  # Note: slots are ordered.

    def __init__(self, input_name=None):
        """
        :param input_name: Name of the gene (symbol, synonym, external reference)
        :type input_name: :class:`str`
        """
        self.input_name = input_name
        self.type_of_match = None
        self.ncbi_id = None
        self._possible_hits = []

    def __getattr__(self, attribute):
        if attribute not in self.__slots__:
            return ''

    @property
    def possible_hits(self):
        return self._possible_hits

    @possible_hits.setter
    def possible_hits(self, genes):
        for gene in genes:
            possible_match = Gene()
            possible_match.ncbi_id = gene
            self._possible_hits.append(possible_match)

    @lru_cache(10000)
    def load_ncbi_info(self):
        """ Populate :class:`Gene` with data from NCBI gene database
        """
        if not self.ncbi_id:
            return

        info = GeneInfoDB().select_gene_info(self.ncbi_id)
        for attr, value in zip(self.__slots__, info):
            if attr == 'db_refs':
                value = parse_sources(value)
            elif attr == 'synonyms':
                value = parse_synonyms(value)

            setattr(self, attr, value)

    def to_list(self):
        _, header_tags = GENE_MATCHER_HEADER

        def parse_attribute(tag):
            gene_attr = getattr(self, '{}'.format(tag))

            if isinstance(gene_attr, dict):
                # note: db_refs are stored as dicts
                gene_attr = ', '.join('{}: {}'.format(key, val)
                                      for (key, val) in gene_attr.items()) if gene_attr else ' '
            elif isinstance(gene_attr, list):
                # note: synonyms are stored as lists
                gene_attr = ', '.join(gene_attr) if gene_attr else ' '

            return gene_attr

        return [parse_attribute(tag) for tag in header_tags]

    def to_html(self):
        """ Return a basic HTML representation of :class:`Gene` class.

        :rtype: :class:`str`
        """
        self.load_ncbi_info()
        db_refs = getattr(self, 'db_refs')
        external_links = []
        if db_refs:
            external_links = ['<dd>- <b>{}</b>: {}<dd>'.format(ref, ref_id) for ref, ref_id in db_refs.items()]

        html_string = """<span>
                      <h3><b>Gene info summary</b></h3>
                      <dl>
                      <dt><b>Gene ID:</b></dt>
                      <dd>- {}</dd>
                      
                      <dt><b>Symbol:</b></dt>
                      <dd>- {}</dd>
   
                      <dt><b>Synonyms:</b></dt>
                      <dd>- {}</dd> 
                      
                      <dt><b>External references:</b></dt>
                      {}
                      
                      <dt><b>Description:</b></dt>
                      <dd>{}</dd>
                      
                      <dt><b>Type of gene:</b></dt>
                      <dd>{}</dd>
                      </dl>
                      </span>""".format(self.ncbi_id,
                                        getattr(self, 'symbol'),
                                        ', '.join([synonym for synonym in getattr(self, 'synonyms')]),
                                        '\n'.join([link for link in external_links]),
                                        getattr(self, 'description'),
                                        getattr(self, 'type_of_gene'))

        return html_string


class GeneInfo(dict):
    def __init__(self, organism):
        """ Load genes for given organism.

        Each instance of :class:`Gene` is mapped to corresponding Entrez ID

        :param organism: Taxonomy if (NCBI taxonomy database)
        :type organism:  class:`str` or :class:`int`
        """
        super().__init__()
        self.__init_gene_info(organism)

    def __init_gene_info(self, organism):
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
        """ Return the :class:`Gene` object for given gene_id

        :param gene_id: Entrez ID (NCBI gene database)
        :type gene_id: str or int

        :rtype: :class:`Gene`
        """

        try:
            return self[int(gene_id)]
        except KeyError as e:
            # TODO: handle this properly?
            # raise e
            pass


class GeneMatcher:

    def __init__(self,  tax_id, **kwargs):
        """ Gene name matching interface

        :param tax_id: Taxonomy if (NCBI taxonomy database)
        :type tax_id: str

        """
        self._organism = ensure_type(tax_id, str)  # type: str
        self._genes = []                           # type: (List[Union[str, Gene]])

        self._case_insensitive = kwargs.get("case_insensitive", False)
        self._matcher = self.load_matcher_file(DOMAIN, MATCHER_FILENAME.format(tax_id))

    @property
    def organism(self):
        return self._organism

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
        """ Return genes with known Entrez Id from NCBI gene database

        Unknown genes are not included.

        :rtype: :class:`list` of :class:`Gene` instances
        """
        return [gene for gene in self.genes if gene.ncbi_id]

    def map_input_to_ncbi(self):
        """ Return a map where key is a given input (gene name) and value is gene id (Entrez Id)

        Unknown genes are not included.

        :rtype: :class:`dict`
        """
        return {gene.input_name: gene.ncbi_id for gene in self.genes if gene.ncbi_id}

    def match_table_attributes(self, data_table):
        """ Helper function for gene name matching in data table.

        Match table attributes and if a unique match is found create a new column attribute for Entrez Id.
        Attribute name is defined here: `orangecontrib.bioinformatics.ncbi.gene.config.NCBI_ID`


        :param data_table: data table
        :type data_table: :class:`Orange.data.Table`

        """
        input_gene_names = [var.name for var in data_table.domain.attributes]

        if input_gene_names:
            self.genes = input_gene_names
            self.run_matcher()

            for gene in self.genes:
                if gene.ncbi_id:
                    data_table.domain[gene.input_name].attributes[NCBI_ID] = gene.ncbi_id

    @staticmethod
    def gene_match_status(gene):
        # type: (Gene) -> str

        if gene.ncbi_id:
            return 'Matched'
        elif gene.possible_hits:
            return 'Match Conflict'
        else:
            return 'Unmatched'

    def to_data_table(self, selected_genes=None):
        tax_id = set()
        data_x = []
        metas = [
            StringVariable('Input gene ID'),
            DiscreteVariable('Match result', values=['Matched', 'Match Conflict', 'Unmatched']),
            StringVariable(NCBI_ID),
            StringVariable('Symbol'),
            StringVariable('Synonyms'),
            StringVariable('Description'),
            StringVariable('Other IDs'),
            StringVariable('Type of gene'),
            StringVariable('Chromosome'),
            StringVariable('Map location'),
            StringVariable('Locus tag'),
            StringVariable('Symbol from nomenclature authority'),
            StringVariable('Full name from nomenclature authority'),
            StringVariable('Nomenclature status'),
            StringVariable('Other designations'),
            StringVariable('Taxonomy ID'),
        ]
        domain = Domain([], metas=metas)

        genes = self.genes
        if selected_genes is not None:
            genes = [gene for gene in self.genes if str(gene.ncbi_id) in selected_genes]

        for gene in genes:
            gene.load_ncbi_info()
            tax_id.add(gene.tax_id)
            match_status = self.gene_match_status(gene)

            db_refs = ', '.join('{}: {}'.format(key, val)
                                for (key, val) in gene.db_refs.items()) if gene.db_refs else ''
            synonyms = ', '.join(gene.synonyms) if gene.synonyms else ''

            line = [
                gene.input_name, match_status, gene.ncbi_id, gene.symbol, synonyms, gene.description, db_refs,
                gene.type_of_gene, gene.chromosome, gene.map_location, gene.locus_tag,
                gene.symbol_from_nomenclature_authority, gene.full_name_from_nomenclature_authority,
                gene.nomenclature_status, gene.other_designations, gene.tax_id
            ]

            data_x.append(line)

        tax_id = filter(None.__ne__, tax_id)

        table = Table(domain, data_x)
        table.name = 'Gene Matcher Results'
        table.attributes[TAX_ID] = next(tax_id)
        table.attributes[GENE_AS_ATTRIBUTE_NAME] = False
        table.attributes[GENE_ID_COLUMN] = NCBI_ID
        return table

    def run_matcher(self, progress_callback=None):
        """ This will try to match genes, with ncbi ids, based on provided input of genes.

        :param progress_callback: Used for progress bar in widgets. It emits the signal back to main thread

        """

        if progress_callback:
            self._match(callback=progress_callback.emit)
        else:
            self._match()

    def _match(self, **kwargs):
        # TODO: refactor this madness ...
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
                    gene.load_ncbi_info()
                    continue
            except ValueError:
                # NCBI ids are stored as Integers. If ValueError is raised, probably not NCBI ID.
                # We expect unique match here.
                pass

            input_name = gene.input_name
            if self._case_insensitive:
                input_name = gene.input_name.lower()

            source_match = match_input(self._matcher[MAP_SOURCES], gene.input_name)
            # ids from different sources are unique. We do not expect to get multiple hits here.
            # There is exceptions with organism 3702. It has same source id from Araport or TAIR databases.
            if source_match:
                gene.ncbi_id = source_match[0][MAP_GENE_ID]
                gene.type_of_match = _source
                gene.load_ncbi_info()
                continue

            symbol_match = match_input(self._matcher[MAP_SYMBOL], input_name)
            if len(symbol_match) == _single_hit:
                gene.ncbi_id = symbol_match[0][MAP_GENE_ID]
                gene.type_of_match = _symbol
                gene.load_ncbi_info()
                continue
            elif len(symbol_match) >= _multiple_hits:
                gene.possible_hits = symbol_match
                continue

            locus_match = match_input(self._matcher[MAP_LOCUS], gene.input_name)
            if len(locus_match) == _single_hit:
                gene.ncbi_id = locus_match[0][MAP_GENE_ID]
                gene.type_of_match = _locus
                gene.load_ncbi_info()
                continue
            elif len(symbol_match) >= _multiple_hits:
                gene.possible_hits = locus_match
                continue

            synonym_match = match_input(self._matcher[MAP_SYNONYMS], input_name)
            if len(synonym_match) == _single_hit:
                gene.ncbi_id = synonym_match[0][MAP_GENE_ID]
                gene.type_of_match = _synonym
                gene.load_ncbi_info()
                continue
            elif len(synonym_match) >= _multiple_hits:
                gene.possible_hits = synonym_match
                continue

            nomenclature_match = match_input(self._matcher[MAP_NOMENCLATURE], input_name)
            if len(nomenclature_match) == _single_hit:
                gene.ncbi_id = nomenclature_match[0][MAP_GENE_ID]
                gene.type_of_match = _nom_symbol
                gene.load_ncbi_info()
                continue
            elif len(nomenclature_match) >= _multiple_hits:
                gene.possible_hits = nomenclature_match
                continue

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
                updated_dict[MAP_SYMBOL][key.lower()] = value

            for key, value in matcher_dict[MAP_SYNONYMS].items():
                key = ensure_type(str(key), str)
                updated_dict[MAP_SYNONYMS][key] = value
                updated_dict[MAP_SYNONYMS][key.lower()] = value

            for key, value in matcher_dict[MAP_NOMENCLATURE].items():
                key = ensure_type(str(key), str)
                updated_dict[MAP_NOMENCLATURE][key] = value
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
