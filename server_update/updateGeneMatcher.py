import sqlite3
import bz2
import pickle


from collections import defaultdict
from server_update import *
from server_update.tests.test_GeneInfo import GeneInfo
from orangecontrib.bioinformatics.ncbi.gene import (
    DOMAIN, FILENAME, MATCHER_FILENAME, MATCHER_TITLE, MATCHER_TAGS,
    MAP_GENE_ID, MAP_SOURCES, MAP_SYMBOL, MAP_SYNONYMS, MAP_LOCUS, MAP_NOMENCLATURE, MAP_TAX_ID)

from orangecontrib.bioinformatics.ncbi.gene.utils import parse_sources, parse_synonyms, GeneInfoDB
from orangecontrib.bioinformatics.ncbi.taxonomy import common_taxids, common_taxid_to_name
from orangecontrib.bioinformatics.utils import serverfiles


serverfiles.localpath_download(DOMAIN, FILENAME)

# for indexing results from gene info database
tax_id, gene_id, symbol, synonyms, sources, locus_tag, nomenclature = range(7)

domain_path = sf_local.localpath(DOMAIN)
temp_path = os.path.join(domain_path, sf_temp)

db_path = os.path.join(domain_path, FILENAME)


create_folder(temp_path)
create_folder(domain_path)


def parse_gene_record(parent_tax, mapper, gene_record):
    gene = {MAP_TAX_ID:        parent_tax,
            MAP_GENE_ID:       gene_record[gene_id],
            MAP_SYMBOL:        gene_record[symbol],
            MAP_SYNONYMS:      parse_synonyms(gene_record[synonyms]),
            MAP_SOURCES:       parse_sources(gene_record[sources]),
            MAP_LOCUS:         gene_record[locus_tag],
            MAP_NOMENCLATURE:  gene_record[nomenclature]
            }

    # construct gene mapper
    mapper[MAP_SYMBOL][gene[MAP_SYMBOL]].append(gene)
    mapper[MAP_LOCUS][gene[MAP_LOCUS]].append(gene)
    mapper[MAP_GENE_ID][gene[MAP_GENE_ID]].append(gene)
    mapper[MAP_NOMENCLATURE][gene[MAP_NOMENCLATURE]].append(gene)

    for gene_synonym in gene[MAP_SYNONYMS]:
        mapper[MAP_SYNONYMS][gene_synonym].append(gene)

    for source_id in gene[MAP_SOURCES].values():
        mapper[MAP_SOURCES][source_id].append(gene)


print("Creating gene name mapper ...")

con = sqlite3.connect(db_path, timeout=15)
cursor = con.cursor()

for taxonomy_id in common_taxids():
    g_db = GeneInfoDB()
    gene_mapper = {MAP_GENE_ID: defaultdict(list),
                   MAP_SOURCES: defaultdict(list),
                   MAP_SYMBOL: defaultdict(list),
                   MAP_SYNONYMS: defaultdict(list),
                   MAP_LOCUS: defaultdict(list),
                   MAP_NOMENCLATURE: defaultdict(list)}

    for record in g_db.select_gene_matcher_data(taxonomy_id):
        parse_gene_record(taxonomy_id, gene_mapper, record)

    with open(os.path.join(domain_path, MATCHER_FILENAME.format(taxonomy_id)), 'wb') as file:
        pickle.dump(gene_mapper, file, protocol=pickle.HIGHEST_PROTOCOL)

    uncompressed_size = os.stat(os.path.join(domain_path, MATCHER_FILENAME.format(taxonomy_id))).st_size

    with bz2.BZ2File(os.path.join(temp_path, MATCHER_FILENAME.format(taxonomy_id)), mode='w', compresslevel=9) as f:
        shutil.copyfileobj(open(os.path.join(domain_path, MATCHER_FILENAME.format(taxonomy_id)), "rb"), f)

    create_info_file(os.path.join(temp_path, MATCHER_FILENAME.format(taxonomy_id)),
                     domain=DOMAIN,
                     filename=MATCHER_FILENAME.format(taxonomy_id),
                     source=SOURCE_SERVER,
                     title=MATCHER_TITLE + ' for ' + common_taxid_to_name(taxonomy_id),
                     tags=MATCHER_TAGS + [taxonomy_id],
                     uncompressed=uncompressed_size,
                     compression='bz2')


con.close()

helper = SyncHelper(DOMAIN, GeneInfo)

# sync files with remote server
helper.run_tests()
helper.sync_files()

helper.remove_update_folder()
