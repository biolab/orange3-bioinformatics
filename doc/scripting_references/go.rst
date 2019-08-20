.. py:currentmodule:: orangecontrib.bioinformatics.go
.. py:module:: orangecontrib.bioinformatics.go

=============
Gene Ontology
=============


Provides access to `Gene Ontology`_ and its gene annotations.

.. _Gene Ontology: http://geneontology.org/

Class References
=================

.. autoclass:: Ontology()
   :members:
      defined_slims_subsets, named_slims_subset, set_slims_subset, slims_for_term,
      extract_super_graph, extract_sub_graph, __getitem__, __len__, __iter__, __contains__


.. autoclass:: Term()

   .. attribute:: id

      The term id.

   .. attribute:: namespace

      The namespace of the term.

   .. attribute:: def_

      The term definition (Note the use of trailing underscore
      to avoid conflict with a python keyword).

   .. attribute:: is_a

      List of term ids this term is a subterm of (parent terms).

   .. attribute:: related

      List of (rel_type, term_id) tuples with rel_type specifying
      the relationship type with term_id.


.. autoclass:: Annotations()
   :members:
   :member-order: bysource
   :exclude-members:
      set_ontology, load, Load, parse_file, ParseFile, RemapGenes,
      AddAnnotation, DrawEnrichmentGraph, GetAllAnnotations, GetAllGenes,
      GetAnnotatedTerms, GetEnrichedTerms, GetGeneNamesTranslator,
      GetOntology, SetOntology, aliasMapper, allAnnotations,
      geneAnnotations, geneNames, geneNamesDict, termAnnotations


.. autoclass:: AnnotationRecord()
   :members:


Usage
=====

Load the ontology and print out some terms::

   from orangecontrib.bioinformatics import go
   ontology = go.Ontology()
   term = ontology["GO:0097194"] # execution phase of apoptosis

   # print a term
   print(term)

   # access fields by name
   print(term.id, term.name)
   # note the use of underscore due to a conflict with a python def keyword
   print(term.def_)


Searching the annotation (part of :download:`code/go/gene_annotations.py`)

.. literalinclude:: code/go/gene_annotations.py

Term enrichment (part of :download:`code/go/enrichment.py`)

.. literalinclude:: code/go/enrichment.py
