""" server update module """
import subprocess
import shutil
import sys


from urllib.request import urlopen
from server_update.common import *
from orangecontrib.bioinformatics.utils import local_cache
from unittest import TextTestRunner, makeSuite

# create_folder(update_folder)  # before calling utils.serverfiles
from orangecontrib.bioinformatics.utils.serverfiles import ServerFiles, LOCALFILES, PATH

print(PATH)
info_date_fmt = '%Y-%m-%d %H:%M:%S.%f'
http_date_fmt = '%a, %d %b %Y %H:%M:%S %Z'

# End script with exit code 10 if files on the server is up-to-date!
up_to_date = 10

sf_server = ServerFiles()
sf_local = LOCALFILES
sf_temp = 'temp'
_sync_path = 'serverfiles-bio@orange.biolab.si:/serverfiles-bio2'

# File sources
SOURCE_SERVER = 'server_file'  # files on the serverfiles-bio repository
SOURCE_USER = 'user_file'      # user defined files


def sf_last_modified(domain, filename):
    return datetime.strptime(sf_server.info(domain, filename)['datetime'], info_date_fmt)


def http_last_modified(url):
    return datetime.strptime(urlopen(url).headers.get('Last-Modified'), http_date_fmt)


class SyncHelper:

    def __init__(self, domain, test_case):
        self.domain = domain
        self.domain_path = sf_local.localpath(self.domain)
        self.test_case = test_case
        self.runner = TextTestRunner(verbosity=2)
        self.suit = makeSuite(self.test_case)
        self.results = None

    def _listfiles(self):
        """ return files that are ready to be synced to the main repository
        """
        return [fname for domain, fname in sf_local.listfiles(self.domain)]

    def _allinfo(self):
        # delete old __INFO__.txt from the server
        subprocess.call(['rsync', '-av', '--delete', '--include=__INFO__', '--exclude=*', local_cache, _sync_path])

        # create new allinfo file
        info_path = os.path.join(local_cache, '__INFO__')
        with open(info_path, 'wt') as f:
            # we must initialize ServerFiles object again because old one has __INFO__ cached
            json.dump(list(ServerFiles().allinfo().items()), f)

        subprocess.call(["rsync", info_path, _sync_path])

    def _prepare_domain_path(self):
        """ Some tests require uncompressed files to run. So we move ready(compressed) files
        to temporary folder (domain_path/temp). When/if all tests pass, we delete uncompressed files
        and move "ready to sync" files back to domain_path folder.
        """
        # temporary folder name is standardized across all update scripts
        temp_folder = os.path.join(self.domain_path, sf_temp)

        if os.path.exists(temp_folder):
            entries = os.scandir(self.domain_path)
            [os.remove(entry.path) for entry in entries if entry.is_file()]

            prepared_files = os.scandir(temp_folder)
            [shutil.move(file.path, self.domain_path) for file in prepared_files]

            try:
                os.rmdir(temp_folder)
            except OSError as e:
                # Directory not empty, copy failed?
                raise e

    def run_tests(self):
        self.results = self.runner.run(self.suit)

    def check_results(self):
        if self.results:
            return self.results.wasSuccessful()

    @staticmethod
    def remove_update_folder():
        shutil.rmtree(local_cache)

    def sync_files(self):
        """ use rsync to upload file on serverfiles-bio repo (info file included).
        """
        if not self.check_results():
            # self.remove_update_folder()
            print('Tests failed!?')
            sys.exit(1)

        self._prepare_domain_path()

        for file in self._listfiles():
            source_file = os.path.join(self.domain_path, file)
            source_info = os.path.join(self.domain_path, file + '.info')

            dest_files = os.path.join(_sync_path, self.domain)
            subprocess.call(["rsync", source_file, source_info, dest_files])
            print("Files {} and {} synced".format(file, file + '.info'))

        self._allinfo()


if __name__ == "__main__":
    # helper = SyncHelper('temp', str)
    # helper._allinfo()
    pass
