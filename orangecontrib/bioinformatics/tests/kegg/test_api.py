import shutil
import tempfile
import unittest
from types import SimpleNamespace
from unittest import mock

from orangecontrib.bioinformatics.kegg import api as keggapi
from orangecontrib.bioinformatics.kegg import conf as keggconf

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

genome_T01001 = """ # noqa \
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
    s = SimpleNamespace(
        list=SimpleNamespace(
            organism=SimpleNamespace(get=lambda: list_organism),
            pathway=lambda org: {"hsa": SimpleNamespace(get=lambda: list_pathway_hsa)}[org],
        ),
        info=lambda db: {"pathway": SimpleNamespace(get=lambda: info_pathway)}[db],
        get=lambda key: {"genome:T01001": SimpleNamespace(get=lambda: genome_T01001)}[key],
    )
    return s


class TestKeggApi(unittest.TestCase):
    def setUp(self):
        # testcase._tmpdir = tempfile.TemporaryDirectory(prefix="kegg-tests")
        self._tmpdir = tempfile.mkdtemp(prefix="kegg-tests")
        self._old_cache_path = keggconf.params["cache.path"]
        keggconf.params["cache.path"] = self._tmpdir
        self._mock_ctx = mock.patch("orangecontrib.bioinformatics.kegg.api.web_service", mock_service)
        self._mock_ctx.__enter__()
        s = keggapi.web_service()
        assert isinstance(s, SimpleNamespace)

    def tearDown(self):
        self._mock_ctx.__exit__(None, None, None)
        keggconf.params["cache.path"] = self._old_cache_path
        shutil.rmtree(self._tmpdir)

    def test_api(self):
        def mock_kegg_api():
            api = keggapi.KeggApi()
            api.service = mock_service()
            return api

        api = mock_kegg_api()
        self.assertIsNotNone(api.list_organisms())
        self.assertIsNotNone(api.list_pathways("hsa"))
        self.assertIsNotNone(api.info("pathway"))


if __name__ == '__main__':
    unittest.main()
