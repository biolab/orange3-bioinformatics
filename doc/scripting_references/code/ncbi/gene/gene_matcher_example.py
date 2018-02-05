from orangecontrib.bioinformatics.ncbi.gene import GeneMatcher, GENE_INFO_TAGS

# specify input
organism = 9606
genes_symbols_to_match = ['HB1', 'BCKDHB', 'TWIST1']

# initialize gene matcher object
gene_matcher = GeneMatcher(organism)
gene_matcher.genes = genes_symbols_to_match

# run matching process
gene_matcher.run_matcher()

# inspect results
for gene in gene_matcher.genes:
    print("\ninput name: " + gene.input_name,
          "\nid from ncbi: ", gene.ncbi_id,
          "\nmatch type: ", gene.type_of_match
          )
    if gene.ncbi_id is None and gene.possible_hits:
        print('possible_hits: ', [hit.ncbi_id for hit in gene.possible_hits])
