DOMAIN = 'gene'

gene_info_attributes = (
    'species',
    'tax_id',
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
    'modification_date',
    'homology_group_id',
    'homologs',
)

_select_gene_info_columns = """
    gene_info.species, gene_info.tax_id, gene_info.gene_id, gene_info.symbol, gene_info.synonyms,  gene_info.db_refs,
    gene_info.description, gene_info.locus_tag, gene_info.chromosome,  gene_info.map_location, gene_info.type_of_gene,
    gene_info.symbol_from_nomenclature_authority, gene_info.full_name_from_nomenclature_authority,
    gene_info.nomenclature_status, gene_info.other_designations, gene_info.modification_date,
    gene_info.homology_group_id, gene_info.homologs
    """

query_exact = f"""
    SELECT  {_select_gene_info_columns}
    FROM gene_info
    JOIN gene_info_fts on gene_info.rowid = gene_info_fts.rowid
    WHERE gene_info_fts
        MATCH ?
            AND (lower(gene_info_fts.gene_id)=?
             OR  lower(gene_info_fts.symbol)=?
             OR  lower(gene_info_fts.locus_tag)=?
             OR  lower(gene_info_fts.symbol_from_nomenclature_authority)=?)
"""

query = f"""
    SELECT {_select_gene_info_columns}
    FROM gene_info
    JOIN gene_info_fts on gene_info.rowid = gene_info_fts.rowid
    WHERE gene_info_fts MATCH ?
"""

# Pretty strings
ENTREZ_ID = 'Entrez ID'
ENSEMBl_ID = 'Ensembl ID'
GENE_SYMBOL = 'Symbol'
GENE_SYNONYMS = 'Synonyms'
