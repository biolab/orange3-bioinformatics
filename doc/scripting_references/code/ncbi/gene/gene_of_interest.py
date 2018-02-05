gene_of_interest = gene_matcher.genes[0].possible_hits[0]
gene_of_interest.load_ncbi_info()

for tag in GENE_INFO_TAGS:
    print(tag + ':', getattr(gene_of_interest, tag))