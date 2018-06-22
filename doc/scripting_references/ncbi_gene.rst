.. py:currentmodule:: orangecontrib.bioinformatics.ncbi.gene

.. index:: gene name matching
.. index:: gene info
.. index:: NCBI

********************************************************
Gene name matching and ncbi info (:mod:`gene`)
********************************************************

Example
--------

.. literalinclude:: code/ncbi/gene/gene_matcher_example.py

Output::

   input name: HB1
   id from ncbi:  None
   match type:  None
   possible_hits:  [3887, 6331, 8184]

   input name: BCKDHB
   id from ncbi:  594
   match type:  Symbol match

   input name: TWIST1
   id from ncbi:  7291
   match type:  Symbol match


Two out of three genes had a unique match with corresponding ncbi gene id.  Symbol 'HB1' is used in multiple genes
so we store them for further analysis.

One can also display gene information from NCBI database.

.. literalinclude:: code/ncbi/gene/gene_of_interest.py

Output::

    tax_id: 9606
    gene_id: 3887
    symbol: KRT81
    synonyms: |HB1|Hb-1|KRTHB1|MLN137|ghHkb1|hHAKB2-1|
    db_refs: MIM:602153|HGNC:HGNC:6458|Ensembl:ENSG00000205426|Vega:OTTHUMG00000167574
    description: keratin 81
    locus_tag: -
    chromosome: 12
    map_location: 12q13.13
    type_of_gene: protein-coding
    symbol_from_nomenclature_authority: KRT81
    full_name_from_nomenclature_authority: keratin 81
    nomenclature_status: O
    other_designations: keratin, type II cuticular Hb1|K81|MLN 137|ghHb1|hair keratin K2.9|hard keratin, type II, 1|keratin 81, type II|keratin, hair, basic, 1|metastatic lymph node 137 gene protein|type II hair keratin Hb1|type-II keratin Kb21
    modification_date: 20171105


Class References
----------------

.. autoclass:: Gene()
   :members:
   :special-members: __init__

.. autoclass:: GeneInfo()
   :members:
   :special-members: __init__

.. autoclass:: GeneMatcher()
   :members:
   :special-members: __init__
