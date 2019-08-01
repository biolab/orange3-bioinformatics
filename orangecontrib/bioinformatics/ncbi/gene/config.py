DOMAIN = 'gene'

gene_info_attributes = ('species', 'tax_id', 'gene_id', 'symbol', 'synonyms', 'db_refs', 'description', 'locus_tag',
                        'chromosome', 'map_location', 'type_of_gene', 'symbol_from_nomenclature_authority',
                        'full_name_from_nomenclature_authority', 'nomenclature_status', 'other_designations',
                        'modification_date', 'homology_group_id', 'homologs')

_select_gene_info_columns = """
    gene_info.species, gene_info.tax_id, gene_info.gene_id, gene_info.symbol, gene_info.synonyms,  gene_info.db_refs, 
    gene_info.description, gene_info.locus_tag, gene_info.chromosome,  gene_info.map_location, gene_info.type_of_gene,
    gene_info.symbol_from_nomenclature_authority, gene_info.full_name_from_nomenclature_authority,
    gene_info.nomenclature_status, gene_info.other_designations, gene_info.modification_date, 
    gene_info.homology_group_id, gene_info.homologs
    """

# Pretty strings
ENTREZ_ID = 'Entrez ID'
ENSEMBl_ID = 'Ensembl ID'
GENE_SYMBOL = 'Symbol'
GENE_SYNONYMS = 'Synonyms'

# OWGenes settings
NCBI_DETAIL_LINK = 'http://www.ncbi.nlm.nih.gov/sites/entrez?Db=gene&Cmd=ShowDetailView&TermToSearch={}'
ENSEMBL_DETAIL_LINK = 'http://www.ensembl.org/id/{}'
owgenes_header = [
    ['Input ID', 'Entrez ID', 'Name', 'Description', 'Synonyms', 'Other IDs'],
    ['input_identifier', 'gene_id', 'symbol', 'description', 'synonyms', 'db_refs']
]