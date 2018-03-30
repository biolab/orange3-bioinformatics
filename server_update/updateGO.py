""" GO update """
import gzip
import bz2

from collections import defaultdict
from server_update import *
from urllib.request import urlopen

from server_update.tests.test_GO import GOTest

from orangecontrib.bioinformatics.go.config import (
    FTP_URL_ANNOTATIONS, FTP_FILE_ANNOTATIONS, DOMAIN, FILENAME_ANNOTATION, FILENAME_ONTOLOGY, FTP_URL_ONTOLOGY,
    ONTOLOGY_TITLE, ONTOLOGY_TAGS
)

from orangecontrib.bioinformatics.ncbi.taxonomy import common_taxids, Taxonomy, common_taxid_to_name


domain_path = sf_local.localpath(DOMAIN)
temp_path = os.path.join(domain_path, sf_temp)

create_folder(domain_path)
create_folder(temp_path)

# GENE ONTOLOGY
stream = urlopen(FTP_URL_ONTOLOGY, timeout=30)

with open(os.path.join(domain_path, FILENAME_ONTOLOGY), 'wb') as f:
    shutil.copyfileobj(stream, f)
    db_size = os.stat(os.path.join(domain_path, FILENAME_ONTOLOGY)).st_size  # store uncompressed database size

with bz2.BZ2File(os.path.join(temp_path, FILENAME_ONTOLOGY), mode='w', compresslevel=9) as f_compressed:
    shutil.copyfileobj(open(os.path.join(domain_path, FILENAME_ONTOLOGY), 'rb'), f_compressed)

create_info_file(os.path.join(temp_path, FILENAME_ONTOLOGY),
                 domain=DOMAIN,
                 filename=FILENAME_ONTOLOGY,
                 source=SOURCE_SERVER,
                 title=ONTOLOGY_TITLE,
                 tags=ONTOLOGY_TAGS,
                 uncompressed=db_size,
                 compression='bz2')

# GENE ANNOTATIONS

tax_ids = common_taxids()
taxonomy = Taxonomy()
store_lines_by_taxid = defaultdict(list)


stream = urlopen(FTP_URL_ANNOTATIONS, timeout=30)
with open(os.path.join(domain_path, FTP_FILE_ANNOTATIONS), 'wb') as f:
    shutil.copyfileobj(stream, f)


with gzip.open(os.path.join(domain_path, FTP_FILE_ANNOTATIONS), 'rb') as gene2go:
    header = gene2go.readline()

    for line in gene2go:
        split_line = line.decode().split('\t')

        if split_line[0] in tax_ids:
            store_lines_by_taxid[split_line[0]].append(line)
        #else:
            #parent = taxonomy.parent(split_line[0])
            #if parent in tax_ids:
                #store_lines_by_taxid[parent].append(line)


for org, lines in store_lines_by_taxid.items():
    filename = FILENAME_ANNOTATION.format(org)
    FILE_PATH = os.path.join(domain_path, filename)
    TITLE = "GO Annotations for " + common_taxid_to_name(org)
    TAGS = ["gene", "annotation", "ontology", "GO", org]

    with open(FILE_PATH, 'wb') as f:
        f.write(header)
        f.writelines(lines)

    db_size = os.stat(FILE_PATH).st_size  # store uncompressed database size

    with bz2.BZ2File(os.path.join(temp_path, filename), mode='w', compresslevel=9) as f_compressed:
        shutil.copyfileobj(open(os.path.join(domain_path, filename), 'rb'), f_compressed)

    create_info_file(os.path.join(temp_path, filename),
                     domain=DOMAIN,
                     filename=filename,
                     source=SOURCE_SERVER,
                     title=TITLE,
                     tags=TAGS,
                     uncompressed=db_size,
                     compression='bz2')


helper = SyncHelper(DOMAIN, GOTest)

# sync files with remote server
helper.run_tests()
helper.sync_files()

helper.remove_update_folder()
