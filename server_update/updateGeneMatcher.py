import sqlite3
import bz2
import pickle


from collections import defaultdict
from server_update import *
from server_update.tests.test_GeneInfo import GeneInfo
from orangecontrib.bioinformatics.ncbi.gene import (
    DOMAIN, FILENAME, gene_matcher_tuple, MATCHER_FILENAME, MATCHER_TITLE, MATCHER_TAGS,
    MAP_GENE_IDS, MAP_SOURCES, MAP_SYMBOLS, MAP_SYNONYMS, MAP_LOCUS)

from orangecontrib.bioinformatics.ncbi.gene.utils import parse_sources, parse_synonyms, GeneInfoDB
from orangecontrib.bioinformatics.ncbi.taxonomy import common_taxids, common_taxid_to_name
from orangecontrib.bioinformatics.utils import serverfiles


serverfiles.localpath_download(DOMAIN, FILENAME)

tax_id, gene_id, symbol, synonyms, source, locus_tag = 0, 1, 2, 3, 4, 5

domain_path = sf_local.localpath(DOMAIN)
temp_path = os.path.join(domain_path, sf_temp)

db_path = os.path.join(domain_path, FILENAME)


create_folder(temp_path)
create_folder(domain_path)


def parse_gene_record(parent_tax, mapper, gene_record):

    gene = gene_matcher_tuple(parent_tax,
                              gene_record[gene_id],
                              gene_record[symbol],
                              parse_synonyms(gene_record[synonyms]),
                              parse_sources(gene_record[source]),
                              gene_record[locus_tag])

    # construct gene mapper
    mapper[MAP_SYMBOLS][gene.symbol].append(gene)
    mapper[MAP_LOCUS][gene.locus_tag].append(gene)
    mapper[MAP_GENE_IDS][gene.gene_id].append(gene)

    for gene_synonym in gene.synonyms:
        mapper[MAP_SYNONYMS][gene_synonym].append(gene)

    for source_id in gene.sources.values():
        mapper[MAP_SOURCES][source_id].append(gene)


print("Creating gene name mapper ...")

con = sqlite3.connect(db_path, timeout=15)
cursor = con.cursor()

for taxonomy_id in common_taxids():
    g_db = GeneInfoDB()
    gene_mapper = {MAP_GENE_IDS: defaultdict(list),
                   MAP_SOURCES: defaultdict(list),
                   MAP_SYMBOLS: defaultdict(list),
                   MAP_SYNONYMS: defaultdict(list),
                   MAP_LOCUS: defaultdict(list)}

    for record in g_db.select_gene_matcher_data(taxonomy_id):
        parse_gene_record(taxonomy_id, gene_mapper, record)

    with open(os.path.join(domain_path, MATCHER_FILENAME.format(taxonomy_id)), 'wb') as file:
        pickle.dump(gene_mapper, file, protocol=pickle.HIGHEST_PROTOCOL)
        uncompressed_size = os.stat(os.path.join(domain_path, MATCHER_FILENAME.format(taxonomy_id))).st_size

    with bz2.BZ2File(os.path.join(temp_path, MATCHER_FILENAME.format(taxonomy_id)), mode='w', compresslevel=9) as f:
        shutil.copyfileobj(open(os.path.join(domain_path, MATCHER_FILENAME.format(taxonomy_id)), "rb"), f)

    create_info_file(os.path.join(temp_path, MATCHER_FILENAME.format(taxonomy_id)),
                     title=MATCHER_TITLE + ' for ' + common_taxid_to_name(taxonomy_id),
                     tags=MATCHER_TAGS + [taxonomy_id],
                     uncompressed=uncompressed_size, compression='bz2')

con.close()

helper = SyncHelper(DOMAIN, GeneInfo)

# sync files with remote server
helper.run_tests()
helper.sync_files()

helper.remove_update_folder()
