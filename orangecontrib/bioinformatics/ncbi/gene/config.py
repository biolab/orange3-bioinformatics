DOMAIN = 'GeneInfo'
FILENAME = 'NCBI_GeneInfo.sqlite'
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

# PICKLED NAME MAPPER
TITLE_SOURCE_MAPPER = 'Gene mapper from source id to ncbi id'
TAGS_SOURCE_MAPPER = ['source id', 'map', 'ncbi id']
SOURCE_MAPPER_FILENAME = 'source_mapper.pck'

TITLE_SYMBOL_MAPPER = 'Gene mapper from symbol to ncbi id'
TAGS_SYMBOL_MAPPER = ['source id', 'map', 'ncbi id']
SYMBOL_MAPPER_FILENAME = 'symbol_mapper.pck'

TITLE_SYNONYM_MAPPER = 'Gene mapper from synonym to ncbi id'
TAGS_SYNONYM_MAPPER = ['source id', 'map', 'ncbi id']
SYNONYM_MAPPER_FILENAME = 'synonym_mapper.pck'

MAPPER_GENE_TAGS = ['tax_id', 'gene_id', 'symbol', 'synonyms', 'source']

GENE_MATCHER_REQUIRED_FILES = [SOURCE_MAPPER_FILENAME, SYMBOL_MAPPER_FILENAME, SYNONYM_MAPPER_FILENAME]

