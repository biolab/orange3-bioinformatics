from orangecontrib.bioinformatics import go

ontology = go.Ontology()
annotations = go.Annotations("4932", ontology=ontology)

# keys are symbol names, values are Entrez IDs
genes_ids = {'Yta7p': '853186', 'RPN2': '854735', 'RPT2': '851557'}
res = annotations.get_enriched_terms(genes_ids.values())

print(res)
print("Enriched terms:")
for go_id, (genes, p_value, ref) in res.items():
    if p_value < 0.05:
        print(ontology[go_id].name + " with p-value: %.4f " % p_value + ", ".join(genes))

# And again for slims
annotations.ontology.set_slims_subset('goslim_yeast')

res = annotations.get_enriched_terms(genes_ids.values(), slims_only=True)
print("\n\nEnriched slim terms:")
for go_id, (genes, p_value, _) in res.items():
    if p_value < 0.2:
        print(ontology[go_id].name + " with p-value: %.4f " % p_value + ", ".join(genes))
