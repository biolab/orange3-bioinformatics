==============================================
KEGG - Kyoto Encyclopedia of Genes and Genomes
==============================================

.. py:currentmodule:: orangecontrib.bioinformatics.kegg

.. automodule:: orangecontrib.bioinformatics.kegg
   :members:
   :member-order: bysource

DBEntry (:mod:`entry`)
----------------------

The :class:`entry.DBEntry` represents a DBGET databas entry.
The individual KEGG Database interfaces below provide their own
specialization for this base class.

.. autoclass:: orangecontrib.bioinformatics.kegg.entry.DBEntry
   :members:
   :member-order: bysource
   :show-inheritance:


KEGG Databases interface (:mod:`databases`)
-------------------------------------------

.. autoclass:: orangecontrib.bioinformatics.kegg.databases.DBDataBase
   :members:
   :member-order: bysource
   :show-inheritance:


.. autoclass:: orangecontrib.bioinformatics.kegg.databases.GenomeEntry
   :members:
   :member-order: bysource
   :show-inheritance:


.. autoclass:: orangecontrib.bioinformatics.kegg.databases.Genome
   :members:
   :member-order: bysource
   :show-inheritance:


.. autoclass:: orangecontrib.bioinformatics.kegg.databases.GeneEntry
   :members:
   :exclude-members:
      alt_names
   :member-order: bysource
   :show-inheritance:

.. autoclass:: orangecontrib.bioinformatics.kegg.databases.Genes
   :members:
   :member-order: bysource
   :show-inheritance:


.. autoclass:: orangecontrib.bioinformatics.kegg.databases.CompoundEntry
   :members:
   :member-order: bysource
   :show-inheritance:

.. autoclass:: orangecontrib.bioinformatics.kegg.databases.Compound
   :members:
   :member-order: bysource
   :show-inheritance:


.. autoclass:: orangecontrib.bioinformatics.kegg.databases.ReactionEntry
   :members:
   :member-order: bysource
   :show-inheritance:

.. autoclass:: orangecontrib.bioinformatics.kegg.databases.Reaction
   :members:
   :member-order: bysource
   :show-inheritance:


.. autoclass:: orangecontrib.bioinformatics.kegg.databases.EnzymeEntry
   :members:
   :member-order: bysource
   :show-inheritance:

.. autoclass:: orangecontrib.bioinformatics.kegg.databases.Enzyme
   :members:
   :member-order: bysource
   :show-inheritance:


.. autoclass:: orangecontrib.bioinformatics.kegg.databases.PathwayEntry
   :members:
   :member-order: bysource
   :show-inheritance:


.. autoclass:: orangecontrib.bioinformatics.kegg.databases.Pathway
   :members:
   :member-order: bysource
   :show-inheritance:


KEGG Pathway (:mod:`pathway`)
-----------------------------

.. autoclass:: orangecontrib.bioinformatics.kegg.pathway.Pathway
   :members:
   :exclude-members:
      entrys
   :member-order: bysource
   :show-inheritance:


Utilities
---------

.. autoclass:: orangecontrib.bioinformatics.kegg.entry.parser.DBGETEntryParser
   :members:
