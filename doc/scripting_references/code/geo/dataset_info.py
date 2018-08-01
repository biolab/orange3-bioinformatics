""" Documentation script """
import textwrap


from orangecontrib.bioinformatics.geo.dataset import GDSInfo

gds_info = GDSInfo()
gds = gds_info["GDS10"]

print("ID:")
print(gds["dataset_id"])
print("Features: ")
print(gds["feature_count"])
print("Genes:")
print(gds["gene_count"])
print("Organism:")
print(gds["platform_organism"])
print("PubMed ID:")
print(gds["pubmed_id"])
print("Sample types:")

for sample_type in set([sinfo["type"] for sinfo in gds["subsets"]]):
    ss = [sinfo["description"] for sinfo in gds["subsets"] if sinfo["type"] == sample_type]
    print("  %s (%s)" % (sample_type, ", ".join(ss)))

print("")
print("Description:")
print("\n".join(textwrap.wrap(gds["description"], 70)))
