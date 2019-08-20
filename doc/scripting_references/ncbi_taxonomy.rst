========
Taxonomy
========

.. py:currentmodule:: orangecontrib.bioinformatics.ncbi.taxonomy


This module provides access to the `NCBI's organism taxonomy information
<http://www.ncbi.nlm.nih.gov/Taxonomy/>`_ and organism name unification across different modules.

.. autofunction:: name
.. autofunction:: other_names
.. autofunction:: search
.. autofunction:: lineage
.. autofunction:: common_taxids
.. autofunction:: common_taxid_to_name

Usage
=====

The following script takes the list of taxonomy IDs and prints out their name:

.. code-block:: python3

	from orangecontrib.bioinformatics.ncbi import taxonomy

	for taxid in taxonomy.common_taxids():
		print(f"{taxid:<10}{taxonomy.name(taxid):<30}{','.join(name for name in taxonomy.shortname(taxid))}")


The output of the script is::

	6500      Aplysia californica           aplysia
	3702      Arabidopsis thaliana          arabidopsis,thaliana,plant
	9913      Bos taurus                    cattle,cow
	6239      Caenorhabditis elegans        nematode,roundworm
	5476      Candida albicans              thrush,candidiasis,candida
	3055      Chlamydomonas reinhardtii     algae
	7955      Danio rerio                   zebrafish
	44689     Dictyostelium discoideum      dicty,amoeba,slime mold
	7227      Drosophila melanogaster       fly,fruit fly,vinegar fly
	562       Escherichia coli              ecoli,coli,bacterium
	11103     Hepatitis C virus             virus, hepatitis
	9606      Homo sapiens                  human
	10090     Mus musculus                  mouse,mus
	2104      Mycoplasma pneumoniae         bacterium,mycoplasma
	4530      Oryza sativa                  asian rice,rice,cereal,plant
	5833      Plasmodium falciparum         plasmodium,malaria,parasite
	4754      Pneumocystis carinii          pneumonia,fungus
	10116     Rattus norvegicus             rat,laboratory rat
	4932      Saccharomyces cerevisiae      yeast,baker yeast,brewer yeast
	4896      Schizosaccharomyces pombe     yeast,fission yeast
	31033     Takifugu rubripes             fish,pufferfish
	8355      Xenopus laevis                frog,african clawed frog
	4577      Zea mays                      corn,cereal grain,plant
