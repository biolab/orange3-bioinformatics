DOMAIN = 'gene_info'
FILENAME = 'gene_info.sqlite'
TITLE = 'NCBI Gene Information'
TAGS = ['NCBI', 'Gene names', 'Gene Ids', 'Description', 'Essential']
FTP_FILENAME = 'gene_info.gz'
FTP_URL = 'ftp://ftp.ncbi.nlm.nih.gov/gene/DATA/gene_info.gz'

# dictyBase gene names
DICTY_GENE_NAMES = \
    'http://www.dictybase.org/db/cgi-bin/dictyBase/download/download.pl?area=general&ID=gene_information.txt'

# NCBI_INFO - order is important
GENE_INFO_TAGS = ['tax_id',
                  'gene_id',
                  'symbol',
                  'synonyms',
                  'db_refs',
                  'description',
                  'locus_tag',
                  'chromosome',
                  'map_location',
                  'type_of_gene',
                  'symbol_from_nomenclature_authority',
                  'full_name_from_nomenclature_authority',
                  'nomenclature_status',
                  'other_designations',
                  'modification_date']

GENE_MATCHER_TAGS = ['input_name', 'type_of_match', 'ncbi_id', '_possible_hits']

# PICKLED GENE MAPPER
MATCHER_FILENAME = 'gene_mapper_{}'
MATCHER_TITLE = 'Gene mapper'
MATCHER_TAGS = ['source id', 'ncbi id', 'symbols', 'synonyms', 'genes', 'mapping']

MATCHER_TUPLE_TAGS = ['tax_id', 'gene_id', 'symbol', 'synonyms', 'sources', 'locus_tag']
MAP_GENE_IDS = 'gene_ids'
MAP_SOURCES = 'sources'
MAP_SYMBOLS = 'symbols'
MAP_SYNONYMS = 'synonyms'
MAP_LOCUS = 'locus_tag'

# OWGeneInfo settings
GENE_INFO_HEADER_LABELS = ['NCBI ID', 'Ensembl ID', 'Symbol', 'Locus Tag', 'Chromosome', 'Description', 'Synonyms', 'Nomenclature']
NCBI_DETAIL_LINK = 'http://www.ncbi.nlm.nih.gov/sites/entrez?Db=gene&Cmd=ShowDetailView&TermToSearch={}'
ENSEMBL_DETAIL_LINK = 'http://www.ensembl.org/id/{}'
