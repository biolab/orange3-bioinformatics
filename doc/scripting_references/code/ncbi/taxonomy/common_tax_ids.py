import orangecontrib.bioinformatics.ncbi.taxonomy

for taxid in orangecontrib.bioinformatics.ncbi.taxonomy.common_taxids():
    print("%-6s %s" % (taxid, orangecontrib.bioinformatics.ncbi.taxonomy.name(taxid)))
