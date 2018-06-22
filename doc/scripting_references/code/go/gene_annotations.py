from orangecontrib.bioinformatics import go

ontology = go.Ontology()

# Load annotations for yeast.
annotations = go.Annotations("4932", ontology=ontology)

# keys are symbol names, values are Entrez IDs
genes = {'RRB1': '855161', 'OST4': '851366', 'VID27': '855509'}
res = annotations.get_enriched_terms(genes.values())


print(annotations.gene_annotations['855161'])
for a in annotations.gene_annotations['855161']:
    print(ontology[a.go_id].name + " with evidence code " + a.evidence)


# Get all genes annotated to the same terms as RRB1
ids = set([a.go_id for a in annotations.gene_annotations['855161']])
for term_id in ids:
    ants = annotations.get_annotations_by_go_id(term_id)
    genes = set([a.gene_id for a in ants])
    print(", ".join(genes) + " annotated to " + term_id + " " + ontology[term_id].name)

