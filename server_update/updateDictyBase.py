""" Dictybase update """
import bz2

from server_update import *
from server_update.tests.test_Dicty import DictyMutantsTest
from orangecontrib.bioinformatics.dicty import phenotypes
from orangecontrib.bioinformatics.dicty.config import PHENOTYPES_FILENAME, PHENOTYPES_TAGS, DOMAIN, PHENOTYPES_TITLE


""" Update DictyMutants
"""
domain_path = sf_local.localpath(DOMAIN)
temp_path = os.path.join(domain_path, sf_temp)

create_folder(domain_path)
create_folder(temp_path)
localfile = sf_local.localpath(temp_path, PHENOTYPES_FILENAME)

mutants = phenotypes.DictyMutants()
base_mutants = phenotypes.download_mutants()

print("updating dictBase mutants")
with open(localfile, 'wb') as f:
    f.write(base_mutants)
    f.close()

uncompressed = os.stat(os.path.join(domain_path, PHENOTYPES_FILENAME)).st_size  # store uncompressed database size

with bz2.BZ2File(os.path.join(temp_path, PHENOTYPES_FILENAME), mode='w', compresslevel=9) as f_compressed:
    shutil.copyfileobj(open(os.path.join(domain_path, PHENOTYPES_FILENAME), 'rb'), f_compressed)

create_info_file(localfile,
                 title=PHENOTYPES_TITLE,
                 tags=PHENOTYPES_TAGS,
                 source=SOURCE_SERVER,
                 domain=DOMAIN,
                 filename=PHENOTYPES_FILENAME,
                 uncompressed=uncompressed,
                 compression='bz2')

helper = SyncHelper(DOMAIN, DictyMutantsTest)
helper.run_tests()
helper.sync_files()
print("DictyBase mutant updated")

helper.remove_update_folder()
