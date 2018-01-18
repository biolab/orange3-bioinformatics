""" GeneInfo utils """
import os
import sqlite3


from contextlib import closing
from orangecontrib.bioinformatics.ncbi.gene import DOMAIN, FILENAME
from orangecontrib.bioinformatics.utils import serverfiles


class GeneInfoFileNotFound(Exception):
    pass


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
        with closing(self._db_con.cursor()) as cursor:
            return cursor.execute('SELECT * FROM gene_info WHERE tax_id = ?', (organism,)).fetchall()

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
        genes = test.select_genes_by_organism("9606")
        print(type(genes), len(genes), genes[0])

    main()

