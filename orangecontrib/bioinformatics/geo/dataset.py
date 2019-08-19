""" GEO Dataset (GDS) """
from collections import OrderedDict

from orangecontrib.bioinformatics.geo import dataset_all_info, dataset_download


class GDSInfo(OrderedDict):
    __slots__ = ()

    def __init__(self):
        """ Retrieve information about `GEO DataSets <http://www.ncbi.nlm.nih.gov/sites/GDSbrowser>`_.

        The class accesses the Orange server file that either resides on the local computer or
        is automatically retrieved from Orange server. Calls to this class do not access any NCBI's servers.

        The constructor will download information on GEO DataSets that are curated and stored on Orange database
        servers.

        An instance behaves like a dictionary: the keys are GEO DataSets IDs, and the dictionary values
        for is a dictionary providing various information about the particular data set.

        Example
        --------
            >>> info = GDSInfo()
            >>> list(info.keys())[:5]
            ['GDS10', 'GDS100', 'GDS1001', 'GDS1002', 'GDS1003']
            >>> info['GDS10']['title']
            'Type 1 diabetes gene expression profiling'
            >>> info['GDS10']['platform_organism']
            'Mus musculus'

        """
        super().__init__()
        self.update({gds_info['name']: gds_info for _, gds_info in dataset_all_info()})


def get_samples(gds_info: dict):
    return {sample for subset in gds_info['subsets'] for sample in subset['sample_id']}


GDS = dataset_download

if __name__ == '__main__':
    info = GDSInfo()
    print(list(info.keys())[:5])
    print(info['GDS10']['title'])
    print(info['GDS10']['platform_organism'])

    dataset_download('GDS1001')
