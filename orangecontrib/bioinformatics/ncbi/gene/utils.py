""" GeneInfo utils """
import os
import sqlite3

from contextlib import closing
from itertools import islice

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


def iter_slice(it, length=100):
    return [_ for _ in islice(it, length)]


def _fetch_organism_strains(organism):
    """  Return a list of taxonomy ids for given organism and all its strains.
    """
    from orangecontrib.bioinformatics.ncbi.taxonomy.utils import Taxonomy
    tax_obj = Taxonomy()
    strains = tax_obj.get_all_strains(organism)
    # append parent tax_id
    strains.append(organism)
    return strains


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
            return cursor.execute('SELECT * FROM gene_info WHERE gene_id = ?', (gene_id,)).fetchone()

    def select_genes_by_organism(self, organism):
        strains = _fetch_organism_strains(organism)
        with closing(self._db_con.cursor()) as cursor:
            return cursor.execute('SELECT * FROM gene_info '
                                  'WHERE tax_id in ({})'.format(', '.join('?' for _ in strains)), strains).fetchall()

    def select_gene_matcher_data(self, organism):
        strains = _fetch_organism_strains(organism)
        with closing(self._db_con.cursor()) as cursor:
            return cursor.execute('SELECT tax_id, gene_id, symbol, synonyms, db_refs, locus_tag FROM gene_info '
                                  'WHERE tax_id in ({})'.format(', '.join('?' for _ in strains)),
                                  strains).fetchall()
    '''def has_external_reference(self, gene_id):
        with closing(self._db_con.cursor()) as cursor:
            if not cursor.execute('SELECT gene_id FROM gene_match WHERE gene_id = ?', (gene_id,)).fetchall():
                return False
            return True

    def match_symbol(self, gene, tax_id):
        with closing(self._db_con.cursor()) as cursor:
            return cursor.execute('SELECT gene_id, symbol FROM gene_info WHERE symbol = ? AND tax_id = ?',
                                  (gene, tax_id)).fetchall()

    def match_synonym(self, gene, tax_id):
        with closing(self._db_con.cursor()) as cursor:
            synonym_tag = '%|' + gene + '|%'
            return cursor.execute('SELECT gene_id FROM gene_info WHERE synonyms LIKE ? AND tax_id = ?',
                                  (synonym_tag, tax_id)).fetchall()

    def match_source_id(self, gene, tax_id):
        with closing(self._db_con.cursor()) as cursor:
            return cursor.execute('SELECT gene_id, source_id FROM gene_match WHERE source_id = ? AND tax_id = ?',
                                  (gene, tax_id)).fetchall()'''

    def __del__(self):
        self._db_con.close()


if __name__ == "__main__":
    def main():
        test = GeneInfoDB()
        genes = test.select_gene_matcher_data("9606")
        print(type(genes), len(genes), genes[0])


    main()
