import unittest
import tempfile
import shutil

try:
    from unittest import mock
except ImportError:
    import backports.unittest_mock
    backports.unittest_mock.install()
    from unittest import mock

try:
    from types import SimpleNamespace as namespace
except ImportError:
    class namespace(object):
        def __init__(self, **kwargs): self.__dict__.update(kwargs)
        def __repr__(self):
            contents = ",".join("{}={!r}".format(*it)
                                for it in sorted(self.__dict__.items()))
            return "namespace(" + contents + ")"

import doctest

from orangecontrib.bio import kegg
from orangecontrib.bio.kegg import api as keggapi
from orangecontrib.bio.kegg import conf as keggconf


list_organism = """\
T01001\thsa\tHomo sapiens (human)\tEukaryotes;Animals;Vertebrates;Mammals
T00005\tsce\tSaccharomyces cerevisiae (budding yeast)\tEukaryotes;Fungi;Ascomycetes;Saccharomycetes
T00245\tddi\tDictyostelium discoideum (cellular slime mold)\tEukaryotes;Protists;Amoebozoa;Dictyostelium\
"""

list_pathway_hsa = """\
path:hsa00010\tGlycolysis / Gluconeogenesis - Homo sapiens (human)
path:hsa00020\tCitrate cycle (TCA cycle) - Homo sapiens (human)
path:hsa00030\tPentose phosphate pathway - Homo sapiens (human)\
"""

info_pathway = """\
pathway          KEGG Pathway Database
path             Release 81.0+/01-18, Jan 17
                 Kanehisa Laboratories
                 479,620 entries
"""

genome_T01001 = """\
ENTRY       T01001            Complete  Genome
NAME        hsa, HUMAN, 9606
DEFINITION  Homo sapiens (human)
ANNOTATION  manual
TAXONOMY    TAX:9606
  LINEAGE   Eukaryota; Metazoa; Chordata; Craniata; Vertebrata; Euteleostomi; Mammalia; Eutheria; Euarchontoglires; Primates; Haplorrhini; Catarrhini; Hominidae; Homo
DATA_SOURCE RefSeq (Assembly:GCF_000001405.31)
ORIGINAL_DB NCBI
            OMIM
            HGNC
            HPRD
            Ensembl
STATISTICS  Number of protein genes:       20234
            Number of RNA genes:           18981
///
"""


def mock_service():
    s = namespace(
        list=namespace(
            organism=namespace(get=lambda: list_organism),
            pathway=lambda org: {
                "hsa": namespace(get=lambda: list_pathway_hsa)
            }[org],
        ),
        info=lambda db:
            {"pathway": namespace(get=lambda: info_pathway)}[db],
        get=lambda key: {
            "genome:T01001": namespace(get=lambda: genome_T01001)
        }[key],
    )
    return s


def mock_kegg_api():
    api = keggapi.KeggApi()
    api.service = mock_service()
    return api


def load_tests(loader, tests, ignore):
    def setUp(testcase):
        # testcase._tmpdir = tempfile.TemporaryDirectory(prefix="kegg-tests")
        testcase._tmpdir = tempfile.mkdtemp(prefix="kegg-tests")
        testcase._old_cache_path = keggconf.params["cache.path"]
        keggconf.params["cache.path"] = testcase._tmpdir
        testcase._mock_ctx = mock.patch(
            "orangecontrib.bio.kegg.api.web_service",
            mock_service)
        testcase._mock_ctx.__enter__()
        s = keggapi.web_service()
        assert isinstance(s, namespace)

    def tearDown(testcase):
        testcase._mock_ctx.__exit__(None, None, None)
        keggconf.params["cache.path"] = testcase._old_cache_path
        shutil.rmtree(testcase._tmpdir)

    api = mock_kegg_api()

    tests.addTests(
        doctest.DocTestSuite(
            keggapi, optionflags=doctest.ELLIPSIS,
            extraglobs={"api": api},
            setUp=setUp, tearDown=tearDown
        )
    )
    return tests
