""" NCBI GeneInformation module """
import json
import sqlite3
import contextlib
from typing import Dict, List, Tuple, Optional

from Orange.data import Table, Domain, StringVariable

from orangecontrib.bioinformatics.utils import serverfiles
from orangecontrib.bioinformatics.ncbi.taxonomy import species_name_to_taxid
from orangecontrib.bioinformatics.ncbi.gene.config import DOMAIN, ENTREZ_ID, query, query_exact, gene_info_attributes
from orangecontrib.bioinformatics.widgets.utils.data import TableAnnotation


class Gene:
    """ Representation of gene summary. """

    __slots__ = gene_info_attributes + ('input_identifier',)

    def __init__(self, input_identifier: Optional[str] = None):
        """
        If we want to match gene to it's corresponding Entrez ID we must,
        upon class initialization, provide some `input identifier`. This way
        :class:`GeneMatcher` will know what to match it against in Gene Database.

        Parameters
        ----------
        input_identifier : str
            This can be any of the following: symbol, synonym, locus tag, other database id, ...
        """
        self.input_identifier = input_identifier

    def __getattr__(self, attribute):
        if attribute not in self.__slots__:
            return None

    def __repr__(self):
        return f'<Gene symbol={self.symbol}, tax_id={self.tax_id}, gene_id={self.gene_id}>'

    def load_attributes(self, values: Tuple[str, ...], attributes: Tuple[str, ...] = gene_info_attributes):
        for attr, val in zip(attributes, values):
            setattr(self, attr, json.loads(val) if attr in ('synonyms', 'db_refs', 'homologs') else val)

    def homolog_gene(self, taxonomy_id: str) -> Optional[str]:
        """Returns gene homolog for given organism.

        Parameters
        ----------
        taxonomy_id: str
            Taxonomy id of target organism.

        Returns
        -------
        str
            Entrez ID (if available).
        """
        return self.homologs.get(taxonomy_id, None)


class GeneMatcher:
    """ Gene name matching interface. """

    def __init__(self, tax_id: str, progress_callback=None, auto_start=True):
        """

        Parameters
        ----------
        tax_id:: str
            Taxonomy id of target organism.

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
    def genes(self) -> List[Gene]:
        return self._genes

    @genes.setter
    def genes(self, genes: List[str]) -> None:
        self._genes = [Gene(input_identifier=gene) for gene in genes]
        if self._auto_start:
            self._match()

    def get_known_genes(self) -> List[Gene]:
        """Return Genes with known Entrez ID

        Returns
        -------
        :class:`list` of :class:`Gene` instances
            Genes with unique match

        """
        return [gene for gene in self.genes if gene.gene_id]

    def to_data_table(self, selected_genes: Optional[List[str]] = None) -> Table:
        """Transform GeneMatcher results to Orange data table.

        Optionally we can provide a list of genes (Entrez Ids).
        The table on the output will be populated only with provided genes.

        Parameters
        ----------
        selected_genes: list
            List of Entrez Ids

        Returns
        -------
        Orange.data.Table
            Summary of Gene info in tabular format
        """
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

        genes: List[Gene] = self.genes
        if selected_genes is not None:
            selected_genes_set = set(selected_genes)
            genes = [gene for gene in self.genes if str(gene.gene_id) in selected_genes_set]

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
        """Helper function for gene name matching with :class:`Orange.data.Table`.

        Give a column of genes, GeneMatcher will try to map genes to their
        corresponding Entrez Ids.


        Parameters
        ----------
        data_table: :class:`Orange.data.Table`
            Data table

        column_name: str
            Name of the column where gene symbols are located

        target_column: :class:`StringVariable`
            Column where we store Entrez Ids.
            Defaults to StringVariable(ncbi.gene.config.NCBI_ID)

        Returns
        -------
        :class:`Orange.data.Table`
            Data table with a column of Gene Ids
        """

        if column_name in data_table.domain:
            self.genes = data_table.get_column_view(column_name)[0]

            if target_column is None:
                target_column = StringVariable(ENTREZ_ID)

            new_domain = Domain(
                data_table.domain.attributes, data_table.domain.class_vars, data_table.domain.metas + (target_column,)
            )

            new_data = data_table.transform(new_domain)
            new_data[:, target_column] = [[str(gene.gene_id) if gene.gene_id else '?'] for gene in self.genes]

            return new_data

    def match_table_attributes(self, data_table, rename=False, source_name='Source ID') -> Table:
        """Helper function for gene name matching with :class:`Orange.data.Table`.

        Match table attributes and if a unique match is found create a new column attribute
        for Entrez Id. Attribute name is defined here: `orangecontrib.bioinformatics.ncbi.gene.config.NCBI_ID`


        Parameters
        ----------
        data_table: :class:`Orange.data.Table`
            Data table

        Returns
        -------

        :class:`Orange.data.Table`
            Data table column attributes are populated with Entrez Ids

        """

        # run gene matcher
        self.genes = [var.name for var in data_table.domain.attributes]

        def helper(gene, attribute):
            if gene.gene_id:
                if rename:
                    attribute = attribute.renamed(gene.symbol)
                    attribute.attributes[source_name] = gene.input_identifier

                attribute.attributes[ENTREZ_ID] = gene.gene_id
            return attribute

        attributes = [helper(gene, attr) for gene, attr in zip(self.genes, data_table.domain.attributes)]
        domain = Domain(attributes, data_table.domain.class_vars, data_table.domain.metas)

        return data_table.transform(domain)

    def match_genes(self):
        self._match()

    def _gene_db_path(self):
        return serverfiles.localpath_download(DOMAIN, f'{self.tax_id}.sqlite')

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
                        match = cursor.execute(query_exact, (match_statement,) + tuple([search_param] * 4)).fetchall()
                        # if unique match
                        if len(match) == 1:
                            gene.load_attributes(match[0])
                            continue

                        match = cursor.execute(query, (f'synonyms:"{search_param}"',)).fetchall()
                        synonym_matched_rows = [
                            m for m in match if search_param in (x.lower() for x in json.loads(m[synonyms]))
                        ]
                        # if unique match
                        if len(synonym_matched_rows) == 1:
                            gene.load_attributes(synonym_matched_rows[0])
                            continue

                        match = cursor.execute(query, (f'db_refs:"{search_param}"',)).fetchall()
                        db_ref_matched_rows = [
                            m for m in match if search_param in (x.lower() for x in json.loads(m[db_refs]).values())
                        ]
                        # if unique match
                        if len(db_ref_matched_rows) == 1:
                            gene.load_attributes(db_ref_matched_rows[0])
                            continue


class GeneInfo(dict):
    def __init__(self, tax_id: str):
        """Loads genes for given organism in a dict.

        Each instance of :class:`Gene` is mapped to corresponding Entrez ID

        Parameters
        ----------
        tax_id: str
            Taxonomy id of target organism.

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


def load_gene_summary(tax_d: str, genes: List[Optional[str]]) -> List[Optional[Gene]]:
    gene_db_path = serverfiles.localpath_download(DOMAIN, f'{tax_d}.sqlite')

    # filter NoneTypes
    _genes = [g for g in genes if g]

    with contextlib.closing(sqlite3.connect(gene_db_path)) as con:
        with con as cur:

            gene_map: Dict[str, Gene] = {}
            for gene_info in cur.execute(f'SELECT * FROM gene_info WHERE gene_id in ({",".join(_genes)})').fetchall():
                gene = Gene()
                gene.load_attributes(gene_info)
                gene_map[gene.gene_id] = gene

            return [gene_map.get(gid, None) if gid else None for gid in genes]


if __name__ == "__main__":
    gm = GeneMatcher('9606')
    gm.genes = ['CD4', '614535', 'ENSG00000205426', "2'-PDE", 'HB-1Y']
    print(list(zip(gm.genes, [g.input_identifier for g in gm.genes])))
    _homologs = load_gene_summary('10090', [g.homolog_gene(taxonomy_id='10090') for g in gm.genes])
    print(_homologs)
