===================================
Organism Taxonomy (:mod:`taxonomy`)
===================================

.. py:currentmodule:: orangecontrib.bioinformatics.ncbi.taxonomy


This module provides access to the `NCBI's organism taxonomy information
<http://www.ncbi.nlm.nih.gov/Taxonomy/>`_ and organism name unification across different modules.

.. autofunction:: name
.. autofunction:: other_names
.. autofunction:: search
.. autofunction:: lineage
.. autofunction:: common_taxids
.. autofunction:: common_taxid_to_name
.. autofunction:: taxname_to_taxid

Examples
--------

The following script takes the list of taxonomy IDs and prints out their name:

.. literalinclude:: code/ncbi/taxonomy/common_tax_ids.py

The output of the script is::

    3702   Arabidopsis thaliana
    9913   Bos taurus
    6239   Caenorhabditis elegans
    5476   Candida albicans
    3055   Chlamydomonas reinhardtii
    7955   Danio rerio
    352472 Dictyostelium discoideum AX4
    7227   Drosophila melanogaster
    562    Escherichia coli
    11103  Hepatitis C virus
    9606   Homo sapiens
    10090  Mus musculus
    2104   Mycoplasma pneumoniae
    4530   Oryza sativa
    5833   Plasmodium falciparum
    4754   Pneumocystis carinii
    10116  Rattus norvegicus
    4932   Saccharomyces cerevisiae
    4896   Schizosaccharomyces pombe
    31033  Takifugu rubripes
    8355   Xenopus laevis
    4577   Zea mays
