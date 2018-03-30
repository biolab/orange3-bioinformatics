""" GEO update """
import ftplib
import pickle
import re
import time

from server_update import *
from orangecontrib.bioinformatics.ncbi import taxonomy
from orangecontrib.bioinformatics import geo
from server_update.tests.test_GEO import GEOTest

force_update = False
DOMAIN = 'GEO'
GDS_INFO = 'gds_info.pickled'
TITLE = 'Gene Expression Omnibus data sets information'
TAGS = ['Gene Expression Omnibus', 'data sets', 'GEO', 'GDS']

FTP_NCBI = 'ftp.ncbi.nih.gov'
NCBI_DIR = 'pub/geo/DATA/SOFT/GDS'
domain_path = sf_local.localpath(DOMAIN)
localfile = sf_local.localpath(DOMAIN, GDS_INFO)
create_folder(domain_path)

try:
    sf_local.localpath_download(DOMAIN, GDS_INFO)
    # read the information from the local file
    with open(localfile, 'rb') as f:
        gds_info, excluded = pickle.load(f, encoding='latin1')
        f.close()

except FileNotFoundError as e:
    print('{} file on the server not found!'.format(GDS_INFO))
    force_update = True


# if needed to refresh the data base
if force_update:
    gds_info, excluded = ({}, {})

# list of common organisms may have changed, rescan excluded list
excluded = dict([(id, taxid) for id, taxid in excluded.items()
                 if taxid not in taxonomy.common_taxids()])
excluded.update([(id, info["taxid"]) for id, info in gds_info.items()
                 if info["taxid"] not in taxonomy.common_taxids()])
gds_info = dict([(id, info) for id, info in gds_info.items()
                 if info["taxid"] in taxonomy.common_taxids()])

# get the list of GDS files from NCBI directory
print("Retrieving ftp directory ...")
ftp = ftplib.FTP(FTP_NCBI)
ftp.login()
ftp.cwd(NCBI_DIR)
dirlist = []
ftp.dir(dirlist.append)

m = re.compile("GDS[0-9]*")
gds_names = [m.search(d).group(0) for d in dirlist if m.search(d)]
gds_names = [name for name in gds_names if not (name in gds_info or name in excluded)]
print('{} new files will be added!'.format(len(gds_names)))
skipped = []

helper = SyncHelper(DOMAIN, GEOTest)

if len(gds_names):
    for count, gds_name in enumerate(gds_names):
        print("%3d of %3d -- Adding %s ..." % (count + 1, len(gds_names), gds_name))
        try:
            time.sleep(1)
            gds = geo.GDS(gds_name)
            if gds.info["taxid"] not in taxonomy.common_taxids():
                excluded[gds_name] = gds.info["taxid"]
                print("... excluded (%s)." % gds.info["sample_organism"])
            else:
                gds_info.update({gds_name: gds.info})
                with open(localfile, 'wb') as f:
                    pickle.dump((gds_info, excluded), f, True)
                    f.close()
                print("... added.")
        except Exception as ex:
            print("... skipped (error):", str(ex))
            skipped.append(gds_name)

    # update .info file
    create_info_file(localfile,
                     domain=DOMAIN,
                     filename=GDS_INFO,
                     source=SOURCE_SERVER,
                     title=TITLE,
                     tags=TAGS)

    print("GDS data sets: %d" % len(gds_info))
    print("Organisms:")

    organisms = [info["sample_organism"] for info in gds_info.values()]
    for org in set(organisms):
        print("  %s (%d)" % (org, organisms.count(org)))

    # sync files with remote server
    helper.run_tests()
    helper.sync_files()
    helper.remove_update_folder()

else:
    helper.remove_update_folder()
    sys.exit(up_to_date)

