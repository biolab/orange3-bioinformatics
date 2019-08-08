""" NCBI GeneInformation module """
import json
import sqlite3
import contextlib

from typing import List, Optional, Tuple, Dict, Union

from Orange.data import StringVariable, Domain, Table

from orangecontrib.bioinformatics.utils import serverfiles
from orangecontrib.bioinformatics.ncbi.taxonomy import species_name_to_taxid
from orangecontrib.bioinformatics.widgets.utils.data import TableAnnotation
from orangecontrib.bioinformatics.ncbi.gene.config import (
    gene_info_attributes,
    _select_gene_info_columns,
    owgenes_header,
    ENTREZ_ID,
    DOMAIN,
    NCBI_DETAIL_LINK,
)


class Gene:
    """ Base Gene class. """

    __slots__ = gene_info_attributes + ('input_identifier',)

    def __init__(self, input_identifier=None):
        self.input_identifier = input_identifier

    def __getattr__(self, attribute):
        if attribute not in self.__slots__:
            return None

    def __repr__(self):
        return f'<Gene symbol={self.symbol}, tax_id={self.tax_id}, gene_id={self.gene_id}>'

    def load_attributes(self, values: Tuple[str, ...], attributes: Tuple[str, ...] = gene_info_attributes):
        for attr, val in zip(attributes, values):
            setattr(self, attr, json.loads(val) if attr in ('synonyms', 'db_refs', 'homologs') else val)

    def to_list(self) -> List[str]:
        _, header_tags = owgenes_header

        def parse_attribute(tag):
            gene_attr = getattr(self, '{}'.format(tag))

            if isinstance(gene_attr, dict):
                # note: db_refs are stored as dicts
                gene_attr = (
                    ', '.join('{}: {}'.format(key, val) for (key, val) in gene_attr.items()) if gene_attr else ' '
                )
            elif isinstance(gene_attr, list):
                # note: synonyms are stored as lists
                gene_attr = ', '.join(gene_attr) if gene_attr else ' '

            return gene_attr

        return [parse_attribute(tag) for tag in header_tags]

    def homolog_gene(self, taxonomy_id: str) -> Union[Dict[str, str], str, None]:
        """ Returns gene homologs.

        :param taxonomy_id: target organism.
        :return:
        """
        return self.homologs.get(taxonomy_id, None)


class GeneMatcher:
    def __init__(self, tax_id: str, progress_callback=None, auto_start=True):
        """ Gene name matching interface

        :param tax_id: Taxonomy id (from NCBI taxonomy database)
        :type tax_id: str

        """
        self._tax_id: str = tax_id
        self._genes: List[Gene] = []
        self._progress_callback = progress_callback
        self._auto_start = auto_start
        self.gene_db_path = self._gene_db_path()

    @property
    def tax_id(self):
        return self._tax_id

    @tax_id.setter
    def tax_id(self, tax_id: str) -> None:
        self._tax_id = tax_id
        self.gene_db_path = self._gene_db_path()

    @property
    def genes(self):
        return self._genes

    @genes.setter
    def genes(self, genes: List[str]) -> None:
        self._genes = [Gene(input_identifier=gene) for gene in genes]
        if self._auto_start:
            self._match()

    def get_known_genes(self) -> List[Gene]:
        """ Return genes with known Entrez Id from NCBI gene database

        Unknown genes are not included.

        :rtype: :class:`list` of :class:`Gene` instances
        """
        return [gene for gene in self.genes if gene.gene_id]

    def to_data_table(self, selected_genes=None):
        data_x = []
        metas = [
            StringVariable('Input gene ID'),
            StringVariable(ENTREZ_ID),
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
            StringVariable('Species'),
            StringVariable('Taxonomy ID'),
        ]
        domain = Domain([], metas=metas)

        genes = self.genes
        if selected_genes is not None:
            genes = [gene for gene in self.genes if str(gene.gene_id) in selected_genes]

        for gene in genes:
            db_refs = (
                ', '.join('{}: {}'.format(key, val) for (key, val) in gene.db_refs.items()) if gene.db_refs else ''
            )
            synonyms = ', '.join(gene.synonyms) if gene.synonyms else ''

            line = [
                gene.input_identifier,
                gene.gene_id,
                gene.symbol,
                synonyms,
                gene.description,
                db_refs,
                gene.type_of_gene,
                gene.chromosome,
                gene.map_location,
                gene.locus_tag,
                gene.symbol_from_nomenclature_authority,
                gene.full_name_from_nomenclature_authority,
                gene.nomenclature_status,
                gene.other_designations,
                species_name_to_taxid(gene.species),
                gene.tax_id,
            ]

            data_x.append(line)

        table = Table(domain, data_x)
        table.name = 'Gene Matcher Results'
        table.attributes[TableAnnotation.tax_id] = self.tax_id
        table.attributes[TableAnnotation.gene_as_attr_name] = False
        table.attributes[TableAnnotation.gene_id_column] = ENTREZ_ID
        return table

    def match_table_column(
        self, data_table: Table, column_name: str, target_column: Optional[StringVariable] = None
    ) -> Table:
        """ Helper function for gene name matching in data table.

        :param data_table: data table
        :param column_name: Name of the column where gene symbols are located
        :param target_column: Column
        :type data_table: :class:`Orange.data.Table`

        """

        if column_name in data_table.domain:
            self.genes = data_table.get_column_view(column_name)[0]

            if target_column is None:
                target_column = StringVariable(ENTREZ_ID)

            domain = Domain([], metas=[target_column])
            data = [[str(gene.gene_id) if gene.gene_id else '?'] for gene in self.genes]
            table = Table(domain, data)

            return Table.concatenate([data_table, table])

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

            for gene in self.genes:
                if gene.gene_id:
                    data_table.domain[gene.input_identifier].attributes[ENTREZ_ID] = gene.gene_id

    def match_genes(self):
        self._match()

    def _gene_db_path(self):
        return serverfiles.localpath_download(DOMAIN, f'{self.tax_id}.sqlite')

    def _query_exact(self):
        return f""" 
            SELECT  {_select_gene_info_columns}
            FROM gene_info
            JOIN gene_info_fts on gene_info.rowid = gene_info_fts.rowid
            WHERE gene_info_fts 
                MATCH ?
                    AND (lower(gene_info_fts.gene_id)=?
                     OR  lower(gene_info_fts.symbol)=?
                     OR  lower(gene_info_fts.locus_tag)=?
                     OR  lower(gene_info_fts.symbol_from_nomenclature_authority)=?)
        """

    def _query(self):
        return f""" 
            SELECT {_select_gene_info_columns}
            FROM gene_info
            JOIN gene_info_fts on gene_info.rowid = gene_info_fts.rowid
            WHERE gene_info_fts MATCH ?
        """

    def _match(self):
        synonyms, db_refs = 4, 5

        with contextlib.closing(sqlite3.connect(self.gene_db_path)) as con:
            with con as cursor:
                for gene in self.genes:

                    if self._progress_callback:
                        self._progress_callback()

                    search_param = gene.input_identifier.lower()

                    if search_param:
                        match_statement = (
                            '{gene_id symbol locus_tag symbol_from_nomenclature_authority}:^"' + search_param + '"'
                        )
                        match = cursor.execute(
                            self._query_exact(), (match_statement,) + tuple([search_param] * 4)
                        ).fetchall()
                        # if unique match
                        if len(match) == 1:
                            gene.load_attributes(match[0])
                            continue

                        match = cursor.execute(self._query(), (f'synonyms:"{search_param}"',)).fetchall()
                        synonym_matched_rows = [
                            m for m in match if search_param in [x.lower() for x in json.loads(m[synonyms])]
                        ]
                        # if unique match
                        if len(synonym_matched_rows) == 1:
                            gene.load_attributes(synonym_matched_rows[0])
                            continue

                        match = cursor.execute(self._query(), (f'db_refs:"{search_param}"',)).fetchall()
                        db_ref_matched_rows = [
                            m for m in match if search_param in [x.lower() for x in json.loads(m[db_refs]).values()]
                        ]
                        # if unique match
                        if len(db_ref_matched_rows) == 1:
                            gene.load_attributes(db_ref_matched_rows[0])
                            continue


class GeneInfo(dict):
    def __init__(self, tax_id: str):
        """ Load genes for given organism.

        Each instance of :class:`Gene` is mapped to corresponding Entrez ID

        :param tax_id: Taxonomy id (from NCBI taxonomy database)
        :type tax_id: str
        """
        super().__init__()
        self.tax_id: str = tax_id
        self.gene_db_path: str = self._gene_db_path()

        connection = sqlite3.connect(self.gene_db_path)
        cursor = connection.cursor()

        for gene_info in cursor.execute('SELECT * FROM gene_info').fetchall():
            gene = Gene()
            gene.load_attributes(gene_info)
            self[gene.gene_id] = gene

        cursor.close()
        connection.close()

    def _gene_db_path(self):
        return serverfiles.localpath_download(DOMAIN, f'{self.tax_id}.sqlite')


def load_gene_summary(tax_d: str, genes: List[Union[str, Gene]]) -> List[Gene]:
    gene_db_path = serverfiles.localpath_download(DOMAIN, f'{tax_d}.sqlite')

    # filter NoneTypes
    _genes = [g for g in genes if g]

    with contextlib.closing(sqlite3.connect(gene_db_path)) as con:
        with con as cur:
            if all(isinstance(g, Gene) for g in _genes):
                gene_ids = [g.gene_id for g in _genes]
            else:
                gene_ids = _genes

            gene_map = {}
            for gene_info in cur.execute(f'SELECT * FROM gene_info WHERE gene_id in ({",".join(gene_ids)})').fetchall():
                gene = Gene()
                gene.load_attributes(gene_info)
                gene_map[gene.gene_id] = gene

            return [gene_map.get(gid, None) for gid in genes]


if __name__ == "__main__":
    gm = GeneMatcher('9606')
    gm.genes = ['CD4', '614535', 'ENSG00000205426', "2'-PDE", 'HB-1Y']
    print(list(zip(gm.genes, [g.input_identifier for g in gm.genes])))

    homologs = [g.homolog_gene(taxonomy_id='10090') for g in gm.genes]
    homologs = load_gene_summary('10090', homologs)
    print(homologs)
