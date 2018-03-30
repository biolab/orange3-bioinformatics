""" NCBI Gene Information data update """
import sqlite3
import gzip
import bz2


from server_update import *
from server_update.tests.test_GeneInfo import GeneInfo
from orangecontrib.bioinformatics.ncbi.gene import (
    DOMAIN, FILENAME, FTP_FILENAME, FTP_URL, TITLE, TAGS
)
from orangecontrib.bioinformatics.dicty.config import DICTY_INFO_URL, DICTY_MAPPING_URL
from orangecontrib.bioinformatics.ncbi.taxonomy import common_taxids
from orangecontrib.bioinformatics.ncbi.taxonomy.utils import Taxonomy


# columns indexes
# ftp://ftp.ncbi.nlm.nih.gov/gene/README under "gene_info" section
tax_id, gene_id, symbol, synonyms, db_refs, description = 0, 1, 2, 4, 5, 8
locus_tag, chromosome, map_location, type_of_gene, modification_date = 3, 6, 7, 9, 14
symbol_from_nomenclature_authority, full_name_from_nomenclature_authority = 10, 11
nomenclature_status, other_designations = 12, 13

domain_path = sf_local.localpath(DOMAIN)
temp_path = os.path.join(domain_path, sf_temp)
file_path = os.path.join(domain_path, FILENAME)

create_folder(domain_path)
create_folder(temp_path)
parent_tax_ids = common_taxids()

# we must include all strains of organism. Genes refer to specific strain of an organism
tax_ids = []
map_strain_to_species = {}
tax_obj = Taxonomy()
for parent_id in parent_tax_ids:
    strains = tax_obj.get_all_strains(parent_id)
    tax_ids.append(parent_id)

    if strains:
        [tax_ids.append(strain_id) for strain_id in strains]


init_table = """ 
 CREATE TABLE "gene_info" ( 
    species INTEGER NOT NULL,
    tax_id INTEGER NOT NULL, 
    gene_id INTEGER NOT NULL UNIQUE, 
    symbol INTEGER NOT NULL,
    synonyms TEXT,
    db_refs TEXT,
    description TEXT,
    locus_tag TEXT,
    chromosome TEXT,
    map_location TEXT,
    type_of_gene TEXT,
    symbol_from_nomenclature_authority TEXT,
    full_name_from_nomenclature_authority TEXT,
    nomenclature_status TEXT,
    other_designations TEXT,
    modification_date TEXT,
    PRIMARY KEY(`tax_id`,`gene_id`));

 CREATE TABLE "gene_source" ( 
    tax_id INTEGER NOT NULL, 
    gene_id INTEGER NOT NULL, 
    source TEXT, 
    source_id TEXT,
    FOREIGN KEY (tax_id, gene_id) REFERENCES gene_info(tax_id, gene_id) ON DELETE CASCADE ON UPDATE NO ACTION

  )
"""

gene_info_lines = []
gene_source_lines = []
gene_names_dictyBase = []

# DOWNLOAD FILES
print("Downloading ncbi file ...")
stream = urlopen(FTP_URL, timeout=30)
with open(os.path.join(domain_path, FTP_FILENAME), 'wb') as f:
    shutil.copyfileobj(stream, f)

print("Downloading dictyBase file ...")
stream = urlopen(DICTY_INFO_URL, timeout=30)
with open(os.path.join(domain_path, 'DictyBase'), 'wb') as f:
    shutil.copyfileobj(stream, f)

print("Downloading dictyMappings file ...")
stream = urlopen(DICTY_MAPPING_URL, timeout=30)
with open(os.path.join(domain_path, 'DictyMapping'), 'wb') as f:
    shutil.copyfileobj(stream, f)


# PARSE DICTY BASE FILE
dicty_id, dicty_name = 0, 1
dictyBase_map = dict()

with open(os.path.join(domain_path, 'DictyBase'), 'r') as f:
    f.readline()  # skip header
    for line in f:
        split_line = line.strip().split('\t')
        dictyBase_map[split_line[dicty_id]] = split_line[dicty_name]


# PARSE DICTY MAPPINGS FILE
ddb_id,	ddb_g_id, name, uniProt_id = 0, 1, 2, 3
dictyMapping_map = dict()

with open(os.path.join(domain_path, 'DictyMapping'), 'r') as f:
    f.readline()  # skip header
    for line in f:
        split_line = line.strip().split('\t')
        try:
            dictyMapping_map[split_line[ddb_g_id]] = split_line[uniProt_id]
        except IndexError:
            # no uniprot id for this gene
            pass


# PARSE GENE_INFO FILE
def parse_refs(gene):
    external_ids = gene[db_refs].split('|')

    for ref in external_ids:
        source, source_id = ref.split(':', 1)
        # gene_source_lines.append([int(gene[tax_id]), int(gene[gene_id]), source, source_id])

        if source == 'dictyBase':
            try:
                # Update gene symbol column, we want to have gene symbols from source database for dictyBase
                gene_names_dictyBase.append([int(gene[gene_id]), dictyBase_map[source_id]])

                # add uniProt ID as a source to dicty genes
                if dictyMapping_map[source_id] is not None:
                    gene[db_refs] = gene[db_refs] + '|UniProt:' + dictyMapping_map[source_id]

            except KeyError:
                # dictyBase map is constructed from official database source for dicty.
                # if source_id in ncbi database is not found in that map, it is probably deprecated.
                pass


print("Parsing NCBI gene info file ...")
with gzip.open(os.path.join(domain_path, FTP_FILENAME), 'rb') as info_file:
    info_file.readline()

    for line in info_file:
        split_line = line.decode().split('\t')

        if split_line[0] in tax_ids:
            gene_synonyms = split_line[synonyms]
            split_line[synonyms] = '|' + gene_synonyms + '|' if gene_synonyms != '-' else ''
            gene_info_lines.append(split_line)

            if split_line[db_refs] != '-':
                parse_refs(split_line)


# CREATE SQLITE DB
con = sqlite3.connect(file_path, timeout=15)
cursor = con.cursor()

for table in ["gene_info", "gene_source"]:
    cursor.execute("DROP TABLE IF EXISTS %s" % table)

cursor.executescript(init_table)

cursor.executemany("INSERT INTO gene_info VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)",
                   ((int(tax_obj.get_species(gene[tax_id])),
                     int(gene[tax_id]),
                     int(gene[gene_id]),
                     gene[symbol],
                     gene[synonyms],
                     gene[db_refs],
                     gene[description],
                     gene[locus_tag],
                     gene[chromosome],
                     gene[map_location],
                     gene[type_of_gene],
                     gene[symbol_from_nomenclature_authority],
                     gene[full_name_from_nomenclature_authority],
                     gene[nomenclature_status],
                     gene[other_designations],
                     gene[modification_date])
                    for gene in gene_info_lines))

# cursor.executemany("INSERT INTO gene_source VALUES (?, ?, ?, ?)",
# ([column for column in match] for match in gene_source_lines))

con.commit()

cursor.executemany("UPDATE gene_info SET symbol = ? WHERE gene_id = ?",
                   ([dicty[1], dicty[0]] for dicty in gene_names_dictyBase))

con.commit()
con.close()


# PREPARE SERVERFILES
print("Creating serverfiles ...")
db_size = os.stat(file_path).st_size  # store uncompressed database size

with bz2.BZ2File(os.path.join(temp_path, FILENAME), mode='w', compresslevel=9) as f:
    shutil.copyfileobj(open(os.path.join(domain_path, FILENAME), 'rb'), f)

create_info_file(os.path.join(temp_path, FILENAME),
                 domain=DOMAIN,
                 filename=FILENAME,
                 source=SOURCE_SERVER,
                 title=TITLE,
                 tags=TAGS,
                 uncompressed=db_size,
                 compression='bz2')

helper = SyncHelper(DOMAIN, GeneInfo)

# sync files with remote server
helper.run_tests()
helper.sync_files()

helper.remove_update_folder()
