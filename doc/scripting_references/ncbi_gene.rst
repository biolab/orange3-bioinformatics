.. py:currentmodule:: orangecontrib.bioinformatics.ncbi.gene

.. index:: gene name matching
.. index:: gene info
.. index:: NCBI

=====
Genes
=====

This module is a wrapper around `Gene database provided from NCBI <https://www.ncbi.nlm.nih.gov/gene>`_. It exposes
a simple interface for working with genes in Python. Additionally, it provides a way to map any (almost) kind of gene
identifier to its corresponding Entrez Id.


Usage
=====

.. code-block:: python3

	from orangecontrib.bioinformatics.ncbi.gene import GeneMatcher

	# Notice that we have symbols, synonyms and Ensembel ID here
	genes_of_interest = ['CD4', 'ENSG00000205426', "2'-PDE", 'HB-1Y']

	# Initialize GeneMatcher. Human is our organism of interest.
	gm = GeneMatcher('9606')
	# this will automatically start the process of name matching
	gm.genes = genes_of_interest

	# print results
	for gene, gene_obj in  zip(genes_of_interest, gm.genes):
		print(f"{gene:<20} {gene_obj}")


We are lucky all of the gene names have a unique match in the Gene database. That's great!

.. code-block:: python3

	CD4                  <Gene symbol=CD4, tax_id=9606, gene_id=920>
	ENSG00000205426      <Gene symbol=KRT81, tax_id=9606, gene_id=3887>
	2'-PDE               <Gene symbol=PDE12, tax_id=9606, gene_id=201626>
	HB-1Y                <Gene symbol=HMHB1, tax_id=9606, gene_id=57824>



Now that we have identified our genes, we can explore further. Genes get automatically populated
with additional information from the NCBI database.

.. code-block:: python3

	g = gm.genes[0]

	print(g.synonyms)
	['CD4mut']

	print(g.db_refs)
	{'MIM': '186940', 'HGNC': 'HGNC:1678', 'Ensembl': 'ENSG00000010610'}

	print(g.type_of_gene)
	protein-coding

	print(g.description)
	CD4 molecule


	# look at all the available Gene attributes
	print(g.__slots__)
	('species', 'tax_id', 'gene_id', 'symbol', 'synonyms', 'db_refs', 'description', 'locus_tag', 'chromosome',
	'map_location', 'type_of_gene', 'symbol_from_nomenclature_authority', 'full_name_from_nomenclature_authority',
	'nomenclature_status', 'other_designations', 'modification_date', 'homology_group_id',
	'homologs', 'input_identifier')


We can also access homologs directly from Gene interface:

.. code-block:: python3

	print(g.homologs)
	{'9913': '407098', '10090': '12504', '10116': '24932'}

	print(g.homology_group_id)
	'513'

	# Find homolog in mouse.
	print(g.homolog_gene(taxonomy_id='10090'))
	'12504'



Class References
================

.. autoclass:: Gene()
	:members:
	:special-members: __init__


.. autoclass:: GeneMatcher()
	:members:
	:special-members: __init__


.. autoclass:: GeneInfo()
	:members:
	:special-members: __init__