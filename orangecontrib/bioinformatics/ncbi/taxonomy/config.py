DOMAIN = "taxonomy"
FILENAME = "ncbi-taxonomy.sqlite"
TAXDUMP_URL = "http://ftp.ncbi.nih.gov/pub/taxonomy/taxdump.tar.gz"

DEFAULT_ORGANISM = '9606'

COMMON_NAMES = (
    ("6500",   "Aplysia californica"),
    ("3702",   "Arabidopsis thaliana"),
    ("9913",   "Bos taurus"),
    ("6239",   "Caenorhabditis elegans"),
    ("5476",   "Candida albicans"),
    ("3055",   "Chlamydomonas reinhardtii"),
    ("7955",   "Danio rerio"),
    ("352472", "Dictyostelium discoideum AX4"),
    ("7227",   "Drosophila melanogaster"),
    ("562",    "Escherichia coli"),
    ("11103",  "Hepatitis C virus"),
    ("9606",   "Homo sapiens"),
    ("10090",  "Mus musculus"),
    ("2104",   "Mycoplasma pneumoniae"),
    ("4530",   "Oryza sativa"),
    ("5833",   "Plasmodium falciparum"),
    ("4754",   "Pneumocystis carinii"),
    ("10116",  "Rattus norvegicus"),
    ("4932",   "Saccharomyces cerevisiae"),
    ("4896",   "Schizosaccharomyces pombe"),
    ("31033",  "Takifugu rubripes"),
    ("8355",   "Xenopus laevis"),
    ("4577",   "Zea mays")

)

SHORT_NAMES = {
    "6500":   ["aplysia"],
    "3702":   ["arabidopsis", "thaliana", "plant"],
    "9913":   ["cattle", "cow"],
    "6239":   ["nematode", "roundworm"],
    "5476":   ["thrush", "candidiasis", "candida"],
    "3055":   ["algae"],
    "7955":   ["zebrafish"],
    "352472": ["dicty", "amoeba", "slime mold"],
    "7227":   ["fly", "fruit fly", "vinegar fly"],
    "562":    ["ecoli", "coli", "bacterium"],
    "11103":  ["virus, hepatitis"],
    "9606":   ["human"],
    "10090":  ["mouse", "mus"],
    "2104":   ["bacterium", "mycoplasma"],
    "4530":   ["asian rice", "rice", "cereal", "plant"],
    "5833":   ["plasmodium", "malaria", "parasite"],
    "4754":   ["pneumonia", "fungus"],
    "10116":  ["rat", "laboratory rat"],
    "4932":   ["yeast", "baker yeast", "brewer yeast"],
    "4896":   ["yeast", "fission yeast"],
    "31033":  ["fish", "pufferfish"],
    "8355":   ["frog", "african clawed frog"],
    "4577":   ["corn", "cereal grain", "plant"]
}


COMMON_NAMES_MAPPING = dict(COMMON_NAMES)

