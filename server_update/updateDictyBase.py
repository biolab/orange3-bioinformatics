""" Dictybase update """

from server_update import *
from server_update.tests.test_Dicty import DictyMutantsTest
from orangecontrib.bioinformatics.dicty import phenotypes
from orangecontrib.bioinformatics.dicty.config import PHENOTYPES_FILENAME, PHENOTYPES_TAGS, DOMAIN, PHENOTYPES_TITLE

""" Update DictyMutants
"""
mutants = phenotypes.DictyMutants()
base_mutants = phenotypes.download_mutants()

create_folder(sf_local.localpath(DOMAIN))
localfile = sf_local.localpath(DOMAIN, PHENOTYPES_FILENAME)


print("updating dictBase mutants")
with open(localfile, 'wb') as f:
    f.write(base_mutants)
    f.close()

create_info_file(localfile, title=PHENOTYPES_TITLE, tags=PHENOTYPES_TAGS)

helper = SyncHelper(DOMAIN, DictyMutantsTest)
helper.run_tests()
helper.sync_files()
print("DictyBase mutant updated")

helper.remove_update_folder()
