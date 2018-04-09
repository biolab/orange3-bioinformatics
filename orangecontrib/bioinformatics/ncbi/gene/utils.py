""" GeneInfo utils """
import os
import sqlite3

from contextlib import closing

from orangecontrib.bioinformatics.ncbi.gene import DOMAIN, FILENAME
from orangecontrib.bioinformatics.utils import serverfiles


class GeneInfoFileNotFound(Exception):
    pass


def parse_synonyms(values):
    return [val for val in values.split('|') if val]


def parse_sources(values):
    """ Parse source string from NCBI gene info.

    ftp://ftp.ncbi.nlm.nih.gov/gene/DATA/README:

    bar-delimited set of identifiers in other databases
    for this gene.  The unit of the set is database:value.
    Note that HGNC and MGI include 'HGNC' and 'MGI', respectively,
    in the value part of their identifier.

    Consequently, dbXrefs for these databases will appear like: HGNC:HGNC:1100
    This would be interpreted as database='HGNC', value='HGNC:1100'


    Args:
        values (str): string of gene sources

    Returns:
        :obj:`dict`: Keys are source names, values are source ids

    """
    external_ids = values.split('|')
    out_dict = {}

    for ref in external_ids:
        if ref != '-':

            source_dest, source_id = ref.split(':', 1)
            out_dict[source_dest] = source_id

    return out_dict


class GeneInfoDB:

    def __init__(self):
        db_path = serverfiles.localpath_download(DOMAIN, FILENAME)

        if os.path.isfile(db_path):
            self._db_con = sqlite3.connect(db_path)
        else:
            raise GeneInfoFileNotFound(db_path)

    def __len__(self):
        with closing(self._db_con.cursor()) as cursor:
            return cursor.execute('SELECT COUNT(gene_id) FROM gene_info').fetchone()[0]

    def select_gene_info(self, gene_id):
        with closing(self._db_con.cursor()) as cursor:
            return cursor.execute('SELECT tax_id, gene_id, symbol, synonyms, db_refs, description, locus_tag,'
                                  'chromosome, map_location, type_of_gene, symbol_from_nomenclature_authority,'
                                  'full_name_from_nomenclature_authority, nomenclature_status, other_designations,'
                                  'modification_date FROM gene_info '
                                  'WHERE gene_id = ?', (gene_id,)).fetchone()

    def select_genes_by_organism(self, organism):
        with closing(self._db_con.cursor()) as cursor:
            return cursor.execute('SELECT tax_id, gene_id, symbol, synonyms, db_refs, description, locus_tag,'
                                  'chromosome, map_location, type_of_gene, symbol_from_nomenclature_authority,'
                                  'full_name_from_nomenclature_authority, nomenclature_status, other_designations,'
                                  'modification_date FROM gene_info '
                                  'WHERE species = ? or tax_id = ?', (organism, organism)).fetchall()

    def select_gene_matcher_data(self, organism):
        with closing(self._db_con.cursor()) as cursor:
            return cursor.execute('SELECT tax_id, gene_id, symbol, synonyms, db_refs, locus_tag, '
                                  'symbol_from_nomenclature_authority FROM gene_info '
                                  'WHERE species = ? or tax_id = ?', (organism, organism)).fetchall()

    def __del__(self):
        self._db_con.close()


if __name__ == "__main__":

    def main():
        test = GeneInfoDB()
        genes = test.select_gene_matcher_data("562")
        print(type(genes), len(genes), genes[0])
        sel_genes = test.select_genes_by_organism('562')
        print(type(sel_genes), len(sel_genes), sel_genes[0])

    main()
