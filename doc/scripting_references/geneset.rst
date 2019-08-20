.. py:currentmodule:: orangecontrib.bioinformatics.geneset
.. py:module:: orangecontrib.bioinformatics.geneset

.. index:: gene set
.. index:: gene sets

==========
Gene sets
==========

This module can load either gene sets distributed with Orange
or custom gene sets in the `GMT file format <http://www.molmine.com/magma/fileformats.html>`_.


Loading gene sets
=================

.. autofunction:: orangecontrib.bioinformatics.geneset.list_all

.. autofunction:: orangecontrib.bioinformatics.geneset.load_gene_sets


Supporting functionality
========================

.. autoclass:: orangecontrib.bioinformatics.geneset.GeneSets
   :members:
   :show-inheritance:

.. autoclass:: orangecontrib.bioinformatics.geneset.GeneSet
   :members:


Helper functions to work with serverfiles
==========================================

.. autofunction:: orangecontrib.bioinformatics.geneset.filename

.. autofunction:: orangecontrib.bioinformatics.geneset.filename_parse

