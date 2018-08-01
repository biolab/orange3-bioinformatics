""" GEO Dataset (GDS) """
import gzip
import pickle


from collections import defaultdict
from Orange.data import DiscreteVariable, ContinuousVariable, StringVariable, Domain

from orangecontrib.bioinformatics.geo.utils import (
    create_table, parse_attribute_line, to_float, spots_mean, gds_ensure_downloaded
)
from orangecontrib.bioinformatics.utils import serverfiles
from orangecontrib.bioinformatics.ncbi import taxonomy
from orangecontrib.bioinformatics.geo.config import *


class MapSpotToGene:
    __slots__ = ('spot_id', 'gene_name', 'data')

    def __init__(self, spot_id, gene_name, data):
        """ Store mapping between spot id and gene.
        """
        self.spot_id = spot_id
        self.gene_name = gene_name
        self.data = data


class GDSInfo:
    def __init__(self):
        """ Retrieve infomation about `GEO DataSets <http://www.ncbi.nlm.nih.gov/sites/GDSbrowser>`_.

        The class accesses the Orange server file that either resides on the local computer or
        is automatically retrieved from Orange server. Calls to this class do not access any NCBI's servers.

        Constructor returning the object with GEO DataSets information. The constructor
        will download GEO DataSets information file (gds_info.pickled) from Orange server,
        it will first check the local copy.

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

        path = serverfiles.localpath_download(DOMAIN, GDS_INFO_FILENAME)
        with open(path, "rb") as f:
            self.info, self.excluded = pickle.load(f,  encoding='latin1')

    def keys(self): return self.info.keys()

    def items(self): return self.info.items()

    def values(self): return self.info.values()

    def clear(self): return self.info.clear()

    def __getitem__(self, key): return self.info[key]

    def __setitem__(self, key, item): self.info[key] = item

    def __len__(self): return len(self.info)

    def __iter__(self): return iter(self.info)

    def __contains__(self, key): return key in self.info


class GDS:
    def __init__(self, gds_name, remove_unknown=None):
        """ Retrieval of a specific GEO DataSet as a :obj:`Orange.data.Table`.

        Constructor returns the object that can retrieve GEO DataSet (samples and gene expressions).
        It first checks a local cache directory if the particular data file is loaded locally,
        else it downloads it from `NCBI's GEO FTP site <ftp://ftp.ncbi.nih.gov/pub/geo/DATA/SOFT/GDS/>`_.

        :param gds_name: An NCBI's ID for the data set in the form "GDSn" where "n" is a GDS ID number.

        :param remove_unknown: Remove spots with sample profiles that include unknown values. They are removed
                               if the proportion of samples with unknown values is above the threshold set by
                               ``remove_unknown``. If None, nothing is removed.

        """

        self.gds_name = gds_name
        self.filename = serverfiles.localpath(DOMAIN, self.gds_name + '.soft.gz')
        gds_ensure_downloaded(self.gds_name)

        self.spot2gene = {}
        self.gene2spots = {}

        self.info = None
        self.gds_data = None
        self.parse_file(remove_unknown=remove_unknown)

        taxid = taxonomy.search(self.info["sample_organism"], exact=True)
        self.info["taxid"] = taxid[0] if len(taxid) == 1 else None

        self.genes = sorted(self.gene2spots.keys())
        self.spots = sorted(self.spot2gene.keys())
        self.info["gene_count"] = len(self.genes)

    def parse_file(self, remove_unknown=None):
        """ Parse GDS data file. Create self.info and self.gds_data
        """

        entity_line = '^'
        entity_attribute = '!'
        attribute_separator = '='
        data_header_description = '#'

        info_dict = {'subsets': []}
        subsets_dict = {}
        data = {}

        line_indicator = {
            '^DATABASE': {},
            '^DATASET': info_dict,
            '^SUBSET': subsets_dict
        }

        with gzip.open(self.filename, 'rt') as gds_file:
            line_type = None
            subset_id = None

            for line in gds_file:
                if line.startswith(entity_line):
                    split_line = line.split('=')
                    line_type = split_line[0].strip()

                    if line_type == '^SUBSET':
                        subset_id = split_line[1].strip()
                        line_indicator[line_type][subset_id] = {}

                elif line.startswith(entity_attribute):
                    if attribute_separator in line:
                        key, value = parse_attribute_line(line.strip())
                        if not subset_id:
                            line_indicator[line_type][key] = value
                        else:
                            line_indicator[line_type][subset_id][key] = value

                elif line.startswith(data_header_description):
                    continue

                # data table row
                else:
                    d = line.rstrip().split('\t')

                    if remove_unknown and (float(d[2:].count('null')) / len(d[2:]) > remove_unknown):
                        continue

                    try:
                        data[d[0]] = MapSpotToGene(d[0], d[1], [to_float(v) for v in d[2:]])
                        # gene to spot and spot to genes mappings
                        spot, gene = d[0:2]
                        self.spot2gene[spot] = gene
                        self.gene2spots[gene] = self.gene2spots.get(gene, []) + [spot]
                    except ValueError:
                        # header
                        info_dict['samples'] = line.rstrip().split('\t')[2:]

            for sub_id, sub_dict in subsets_dict.items():
                sub_dict['id'] = sub_id
                sub_dict['sample_id'] = sub_dict['sample_id'].split(',')
                info_dict['subsets'].append(sub_dict)

            self.info = info_dict
            self.gds_data = data

    def _sample_annotations(self, sample_type=None):
        """ Return a dictionary with sample annotation.
        """
        annotation = {}
        for a in self.info['samples']:
            annotation[a] = {}
        for info in self.info['subsets']:
            if not sample_type or info['type'] == sample_type:
                for sample in info['sample_id']:
                    annotation[sample][info['type']] = info['description']
        return annotation

    def _sample_types(self):
        """ Return a set of sample types.
        """
        return set([info["type"] for info in self.info["subsets"]])

    def _sample_to_class(self, sample_type=None):
        """ Return class values for GDS samples.
        """
        annotations = self._sample_annotations(sample_type)
        return dict([(sample, "|".join([a for t, a in sorted(ann.items())])) for sample, ann in annotations.items()])

    def _to_orange_data_table(self, report_genes=True, merge_function=spots_mean, sample_type=None, transpose=False):
        """ Convert parsed GEO format to orange, save by genes or by spots.
        """
        if transpose:  # samples in rows
            sample2class = self._sample_to_class(sample_type)
            cvalues = sorted(set(sample2class.values()))
            if None in cvalues:
                cvalues.remove(None)

            samp_ann = self._sample_annotations()

            ad = defaultdict(set)
            for d in samp_ann.values():
                for n, v in d.items():
                    ad[n].add(v)

            # auto-select sample type if there is only one
            if len(ad) == 1:
                sample_type = list(ad.keys())[0]

            classvar = DiscreteVariable(name=sample_type or "class", values=cvalues)
            spots = self.genes if report_genes else self.spots
            atts = [ContinuousVariable(name=gene) for gene in spots]

            metasvar = [DiscreteVariable(name=n, values=sorted(values))
                        for n, values in ad.items() if n != sample_type]

            X = []
            Y = []
            metas = []
            for (i, sampleid) in enumerate(self.info["samples"]):
                vals = [((merge_function([self.gds_data[spot].data[i] for spot in self.gene2spots[gene]]))
                         if report_genes else self.gds_data[gene].data[i]) for gene in spots]
                X.append(vals)
                Y.append(sample2class.get(sampleid, None))
                metas.append([samp_ann[sampleid].get(n, None) for n, _ in ad.items() if n != sample_type])

            domain = Domain(atts, classvar, metas=metasvar)
            return create_table(domain, X, Y, metas)
        else:  # genes in rows
            annotations = self._sample_annotations(sample_type)
            atts = [ContinuousVariable(name=ss) for ss in self.info["samples"]]
            for i, a in enumerate(atts):
                setattr(a, "attributes", annotations[self.info["samples"][i]])

            geneatname = "gene" if report_genes else "spot"
            metasvar = [StringVariable(geneatname)]
            nameval = self.genes if report_genes else self.spots

            X = []
            metas = []
            for g in nameval:
                if report_genes:
                    X.append(list(map(lambda *x: merge_function(x), *[self.gds_data[spot].data
                                                                      for spot in self.gene2spots[g]])))
                else:
                    X.append(self.gds_data[g].data)
            metas = [[a] for a in nameval]
            domain = Domain(atts, [], metas=metasvar)
            return create_table(domain, X, None, metas)

    def get_data(self, report_genes=True, merge_function=spots_mean, sample_type=None, transpose=False):
        """

        :param bool report_genes: Microarray spots reported in the GEO data set can either be merged according to their
                             gene ids (if True) or can be left as spots.
        :param sample_type: the type of annotation, or (if :obj:`transpose` is True) the type of class labels to be
                            included in the data set. The entire annotation of samples will be included either
                            in the class value or in the ``.attributes`` field of each data set attribute.
        :param transpose: The output table can have spots/genes in rows and samples in columns
                         (False, default) or samples in rows and  spots/genes in columns (True).

        :param merge_function:

        :return: the GEO DataSet as an :obj:`Orange.data.Table`.
        """

        return self._to_orange_data_table(merge_function=merge_function, sample_type=sample_type,
                                          transpose=transpose, report_genes=report_genes)

    def __str__(self):
        return "%s (%s), samples=%s, features=%s, genes=%s, subsets=%d" % \
              (self.info["dataset_id"],
               self.info["sample_organism"],
               self.info['sample_count'],
               self.info['feature_count'],
               self.info['gene_count'],
               len(self.info['subsets']),)
