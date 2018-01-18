import os
import gzip
import re
import io
import urllib.request
import pickle
import numpy

from collections import defaultdict
from Orange.data import DiscreteVariable, ContinuousVariable, StringVariable, Domain, Table
from orangecontrib.bioinformatics.utils import serverfiles
from orangecontrib.bioinformatics.ncbi import taxonomy
from orangecontrib.bioinformatics.geo.config import *


unknown = NAN = float("nan")
isunknown = numpy.isnan


def create_domain(at, cl, metas):
    return Domain(at, cl, metas=metas)


def create_table(domain, X, Y, metas):
    classvar = domain.class_var
    metaatts = domain.metas

    if Y:
        Y = numpy.array([[classvar.to_val(row)] for row in Y],
                        dtype=float)
    if metas:
        metas = numpy.array([[c.to_val(v) for c, v in zip(metaatts, row)]
                             for row in metas],
                            dtype=object)
    data = Table(domain, numpy.asarray(X), Y=Y, metas=metas)

    return data


def spots_mean(x):
    vs = [v for v in x if not isunknown(v)]
    if len(vs) == 0:
        return unknown
    else:
        return sum(vs) / len(vs)


def spots_median(x):
    vs = [v for v in x if not isunknown(v)]
    if len(vs) == 0:
        return unknown
    return numpy.median(vs)


def spots_min(x):
    vs = [v for v in x if not isunknown(v)]
    if len(vs) == 0:
        return unknown
    else:
        return min(vs)


def spots_max(x):
    vs = [v for v in x if not isunknown(v)]
    if len(vs) == 0:
        return unknown
    else:
        return max(vs)


p_assign = re.compile(" = (.*$)")
p_tagvalue = re.compile("![a-z]*_([a-z_]*) = (.*)$")
tagvalue = lambda x: p_tagvalue.search(x).groups()


class GDSInfo:
    """ Retrieve infomation about `GEO DataSets
    <http://www.ncbi.nlm.nih.gov/sites/GDSbrowser>`_.  The class accesses
    the Orange server file that either resides on the local computer or
    is automatically retrieved from Orange server. Calls to
    this class do not access any NCBI's servers.

    Constructor returning the object with GEO DataSets information. If
    `force_update` is True, the constructor will download GEO DataSets
    information file (gds_info.pickled) from Orange server, otherwise
    it will first check the local copy.

    An instance behaves like a dictionary: the keys are GEO DataSets
    IDs, and the dictionary values for is a dictionary providing various
    information about the particular data set.
    """

    def __init__(self, force_update=False):
        path = serverfiles.localpath(DOMAIN, GDS_INFO_FILENAME)
        if not os.path.exists(path) or force_update:
            serverfiles.download(DOMAIN, GDS_INFO_FILENAME)
        f = open(path, "rb")
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


class GeneData:
    """ Store mapping between spot id and gene. """
    def __init__(self, spot_id, gene_name, d):
        self.spot_id = spot_id
        self.gene_name = gene_name
        self.data = d


class GDS:
    """
    Retrieval of a specific GEO DataSet as a :obj:`Orange.data.Table`.

    Constructor returns the object that can retrieve
    GEO DataSet (samples and gene expressions). It first checks
    a local cache directory if the particular data file is
    loaded locally, else it downloads it from `NCBI's GEO FTP site
    <ftp://ftp.ncbi.nih.gov/pub/geo/DATA/SOFT/GDS/>`_.

    :param gdsname: An NCBI's ID for the data set in the form "GDSn"
      where "n" is a GDS ID number.

    :param force_download: Force the download.

    """

    def __init__(self, gdsname, verbose=False, force_download=False):
        self.gdsname = gdsname
        self.verbose = verbose
        self.force_download = force_download
        self.filename = serverfiles.localpath(DOMAIN, self.gdsname + ".soft.gz")
        d = os.path.dirname(self.filename)
        if not os.path.exists(d):
            os.makedirs(d)
        self._getinfo()  # to get the info
        taxid = taxonomy.search(self.info["sample_organism"], exact=True)
        self.info["taxid"] = taxid[0] if len(taxid) == 1 else None
        self._getspotmap()  # to get gene->spot and spot->gene mapping
        self.genes = sorted(self.gene2spots.keys())
        self.spots = sorted(self.spot2gene.keys())
        self.info["gene_count"] = len(self.genes)
        self.gdsdata = None
        self.data = None

    def _download(self):
        """Download GDS data file if not in local cache or forced download requested."""
        localpath = serverfiles.localpath(DOMAIN)
        if self.force_download or not os.path.exists(self.filename):
            base_url = ("ftp://{FTP_NCBI}/{FTP_DIR}"
                        .format(FTP_NCBI=FTP_NCBI, FTP_DIR=FTP_DIR))
            url = base_url + self.gdsname + ".soft.gz"
            targetfn = os.path.join(localpath, self.gdsname + ".soft.gz")
            r = urllib.request.urlopen(url)
            with open(targetfn + "2", 'wb') as f:
                 f.write(r.read())
            with gzip.open(targetfn + "2") as f:
                f.read() #verify the download
            os.rename(targetfn + "2", targetfn)

    def _getinfo(self):
        """Parse GDS data file and return a dictionary with info."""
        getstate = lambda x: x.split(" ")[0][1:]
        getid = lambda x: x.rstrip().split(" ")[2]
        self._download()
        f = gzip.open(self.filename, "rb")
        f = io.TextIOWrapper(f, encoding=SOFT_ENCODING)

        state = None
        previous_state = None

        info = {"subsets": []}
        subset = None

        # GDS information part
        for line in f:
            if line[0] == "^":
                previous_state = state; state = getstate(line)
                if state == "SUBSET":
                    if subset:
                        info["subsets"] += [subset]
                    subset = {"id" : getid(line)}
                if state == "DATASET":
                    info["dataset_id"] = getid(line)
                continue
            if state == "DATASET":
                if previous_state == "DATABASE":
                    tag, value = tagvalue(line)
                    info[tag] = value
                else:
                    if subset:
                        info["subsets"] += [subset]
                    break
            if state == "SUBSET":
                tag, value = tagvalue(line)
                if tag == "description" or tag == "type":
                    subset[tag] = value
                if tag == "sample_id":
                    subset[tag] = value.split(",")
        for t,v in info.items():
            if "count" in t:
                info[t] = int(v)

        # sample information
        state = None
        for line in f:
            if state == "header":
                info["samples"] = line.rstrip().split("\t")[2:]
                break
            if line.startswith("!dataset_table_begin"):
                state = "header"
        self.info = info

    def _getspotmap(self, include_spots=None):
        """Return gene to spot and spot to genes mapings."""
        f = gzip.open(self.filename, "rb")
        f = io.TextIOWrapper(f, encoding=SOFT_ENCODING)
        for line in f:
            if line.startswith("!dataset_table_begin"):
                break
        f.readline()  # skip header
        spot2gene = {}
        gene2spots = {}
        for line in f:
            if line.startswith("!dataset_table_end"):
                break
            spot, gene = line.rstrip().split("\t")[0:2]
            if include_spots and (spot not in include_spots):
                continue
            spot2gene[spot] = gene
            gene2spots[gene] = gene2spots.get(gene, []) + [spot]

        self.spot2gene = spot2gene
        self.gene2spots = gene2spots

    def sample_annotations(self, sample_type=None):
        """Return a dictionary with sample annotation."""
        annotation = {}
        for a in self.info["samples"]:
            annotation[a] = {}
        for info in self.info["subsets"]:
            if not sample_type or info["type"] == sample_type:
                for id in info["sample_id"]:
                    annotation[id][info["type"]] = info["description"]
        return annotation

    def sample_to_class(self, sample_type=None):
        """Return class values for GDS samples."""
        annotations = self.sample_annotations(sample_type)
        return dict([(sample, "|".join([a for t, a in sorted(ann.items())])) for sample, ann in annotations.items()])

    def sample_types(self):
        """Return a set of sample types."""
        return set([info["type"] for info in self.info["subsets"]])

    def _parse_soft(self, remove_unknown=None):
        """Parse GDS data, returns data dictionary."""
        f = gzip.open(self.filename, "rb")
        f = io.TextIOWrapper(f, encoding=SOFT_ENCODING)

        mfloat = lambda x: float(x) if x != 'null' else unknown

        data = {}
        # find the start of the data part
        for line in f:
            if line.startswith("!dataset_table_begin"):
                break
        f.readline()

        # load data
        for line in f:
            if line.startswith("!dataset_table_end"):
                break
            d = line.rstrip().split("\t")
            if remove_unknown and (float(d[2:].count('null')) / len(d[2:]) > remove_unknown):
                continue
            data[d[0]] = GeneData(d[0], d[1], [mfloat(v) for v in d[2:]])
        self.gdsdata = data

    def _to_ExampleTable(self, report_genes=True, merge_function=spots_mean, sample_type=None, transpose=False):
        """Convert parsed GEO format to orange, save by genes or by spots."""
        if transpose:  # samples in rows
            sample2class = self.sample_to_class(sample_type)
            cvalues = sorted(set(sample2class.values()))
            if None in cvalues:
                cvalues.remove(None)

            samp_ann = self.sample_annotations()

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
                vals = [((merge_function([self.gdsdata[spot].data[i] for spot in self.gene2spots[gene]]))
                         if report_genes else self.gdsdata[gene].data[i]) for gene in spots]
                X.append(vals)
                Y.append(sample2class.get(sampleid, None))
                metas.append([samp_ann[sampleid].get(n, None) for n, _ in ad.items() if n != sample_type])

            domain = create_domain(atts, classvar, metasvar)
            return create_table(domain, X, Y, metas)
        else:  # genes in rows
            annotations = self.sample_annotations(sample_type)
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
                    X.append(list(map(lambda *x: merge_function(x),
                        *[self.gdsdata[spot].data for spot in self.gene2spots[g]])))
                else:
                    X.append(self.gdsdata[g].data)
            metas = [[a] for a in nameval]
            domain = create_domain(atts, None, metasvar)
            return create_table(domain, X, None, metas)

    def getdata(self, report_genes=True, merge_function=spots_mean, sample_type=None,
                transpose=False, remove_unknown=None):
        """
        Returns the GEO DataSet as an :obj:`Orange.data.Table`.

        :param report_genes: Microarray spots reported in the GEO data set can
          either be merged according to their gene ids
          (if True) or can be left as spots.

        :param transpose: The output
          table can have spots/genes in rows and samples in columns
          (False, default) or samples in rows and  spots/genes in columns
          (True).

        :param sample_type: the type of annotation, or (if :obj:`transpose` is True)
          the type of class labels to be included in the data set.
          The entire annotation of samples will
          be included either in the class value or in
          the ``.attributes`` field of each data set
          attribute.

        :param remove_unknown: Remove spots with sample profiles that
          include unknown values. They are removed if the proportion
          of samples with unknown values is above the threshold set by
          ``remove_unknown``. If None, nothing is removed.
        """
        if self.verbose: print("Reading data ...")
        # if not self.gdsdata:
        self._parse_soft(remove_unknown=remove_unknown)
        # if remove_unknown:
        # some spots were filtered out, need to revise spot<>gene mappings
        self._getspotmap(include_spots=set(self.gdsdata.keys()))
        if self.verbose: print("Converting to example table ...")
        self.data = self._to_ExampleTable(merge_function=merge_function,
                                          sample_type=sample_type, transpose=transpose,
                                          report_genes=report_genes)
        return self.data

    def __str__(self):
        return "%s (%s), samples=%s, features=%s, genes=%s, subsets=%d" % \
              (self.info["dataset_id"],
               self.info["sample_organism"],
               self.info['sample_count'],
               self.info['feature_count'],
               self.info['gene_count'],
               len(self.info['subsets']),
               )


def _float_or_na(x):
    if isunknown(x):
        return unknown
    else:
        return float(x)


def transpose_class_to_labels(data, attcol="sample"):
    """Converts data with genes as attributes to data with genes in rows."""
    if attcol in [v.name for v in data.domain.getmetas().values()]:
        atts = [ContinuousVariable(str(d[attcol])) for d in data]
    else:
        atts = [ContinuousVariable("S%d" % i) for i in range(len(data))]
    for i, d in enumerate(data):
        atts[i].setattr("class", str(d.getclass()))
    domain = Domain(atts, None)

    newdata = []
    for a in data.domain.attributes:
        newdata.append([_float_or_na(d[a]) for d in data])

    gene = StringVariable("gene")
    id = StringVariable.new_meta_id()
    new = Table(domain, newdata)
    new.domain.addmeta(id, gene)
    for i, d in enumerate(new):
        d[gene] = data.domain.attributes[i].name

    return new


def transpose_labels_to_class(data, class_label=None, gene_label="gene"):
    """Converts data with genes in rows to data with genes as attributes."""
    # if no class_label (attribute type) given, guess it from the data
    if not class_label:
        l = []
        for a in data.domain.attributes:
            l.extend(list(a.attributes.keys()))
        l = list(set(l))
        class_label = l[0]
        if len(set(l)) > 1:
            import warnings
            warnings.warn("More than single attribute label types (%s), took %s"
                          % (", ".join(l), class_label))

    if gene_label in [v.name for v in data.domain.getmetas().values()]:
        atts = [ContinuousVariable(str(d[gene_label])) for d in data]
    else:
        atts = [ContinuousVariable("A%d" % i) for i in range(len(data))]

    classvalues = list(set([a.attributes[class_label] for a in data.domain.attributes]))

    if all([isinstance(x, (int, float, complex)) for x in classvalues]):
        classvar = ContinuousVariable(class_label)
    else:
        classvar = DiscreteVariable(class_label, values=classvalues)

    domain = Domain(atts, classvar)

    newdata = []
    for a in data.domain.attributes:
        newdata.append([_float_or_na(d[a]) for d in data] + [a.attributes[class_label]])

    sample = StringVariable("sample")
    id = StringVariable.new_meta_id()
    new = Table(domain, newdata)
    new.domain.addmeta(id, sample)
    for i, d in enumerate(new):
        d[sample] = data.domain.attributes[i].name

    return new


def transpose(data):
    """ Transposes data matrix, converts class information to attribute label and back. """
    if data.domain.class_var:
        return transpose_class_to_labels(data)
    else:
        return transpose_labels_to_class(data)
