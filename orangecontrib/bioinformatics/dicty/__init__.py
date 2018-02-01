import pickle


from urllib.request import urlopen
from orangecontrib.bioinformatics.geneset.utils import serverfiles
from . import phenotypes


def txt2ll(s, separ=' ', lineSepar='\n'):
    return [a.split(separ) for a in s.split(lineSepar)]


def download_url(url, repeat=2):
    def do():
        return urlopen(url)

    if repeat <= 0:
        do()
    else:
        try:
            return do()
        except:
            return download_url(url, repeat=repeat-1)


def empty_none(s):
    if s:
        return s
    else:
        return None


class DictyBase:
    domain = "dictybase"
    filename = "information_mappings.pck"
    tags = ["Dictyostelium discoideum", "gene", "essential", "dictyBase"]

    @classmethod
    def download_information(cls):
        """
        Downloads gene information and parses it.
        Returns a dictionary {ID: (name, synonyms, products)}
        """
        s = download_url(
            "http://www.dictybase.org/db/cgi-bin/dictyBase/download/download.pl?area=general&ID=gene_information.txt").read()
        out = []
        s = s.decode("utf-8")  # bytes to string
        for l in txt2ll(s, separ='\t', lineSepar='\n')[1:]:
            if len(l) == 4:
                id = l[0]
                name = l[1]
                synonyms = filter(None, l[2].split(", "))
                products = l[3]
                out.append((id, name, synonyms, products))
        return dict((a, (b, c, d)) for a, b, c, d in out)

    @classmethod
    def download_mappings(cls):
        """
        Downloads DDB-GeneID-UniProt mappings and parses them.
        Returns a list of (ddb, ddb_g, uniprot) triplets.

        2009/04/07: ddb's appear unique
        """
        s = download_url(
            "http://www.dictybase.org/db/cgi-bin/dictyBase/download/download.pl?area=general&ID=DDB-GeneID-UniProt.txt").read()
        out = []
        s = s.decode("utf-8")  # bytes to string
        for l in txt2ll(s, separ='\t', lineSepar='\n')[1:]:
            if len(l) == 4:
                ddb = empty_none(l[0])
                ddb_g = empty_none(l[1])
                name = empty_none(l[2])
                uniprot = empty_none(l[3])
                out.append((ddb, ddb_g, name, uniprot))
        return out

    @classmethod
    def pickle_data(cls):
        info = cls.download_information()
        mappings = cls.download_mappings()
        return pickle.dumps((info, mappings), -1)

    def __init__(self):
        fn = serverfiles.localpath_download(self.domain, self.filename)
        self.info, self.mappings = pickle.load(open(fn, 'rb'))

