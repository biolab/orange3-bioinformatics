.. py:currentmodule:: orangecontrib.bioinformatics.dicty.phenotypes

.. index:: mutant phenotypes
.. index:: phenotypes
.. index::
   single: D. dictyostelium; mutants

===============================
D. discoideum Mutant Phenotypes
===============================

This modules provides an interface to `Dictyostelium mutant
phenotypes <http://dictybase.org/Downloads/>`_ data from the
`dictyBase <http://dictybase.org/>`_.  The mutants are presented as
:obj:`DictyMutant` objects with their respective name, strain descriptor,
associated genes and associated phenotypes.


Usage
=====

>>> from orangecontrib.bio.dicty.phenotypes import *
>>> # Create a set of all mutant objects
>>> dicty_mutants = mutants()
>>> # List a set of all genes referenced by a single mutant
>>> print(mutant_genes(dicty_mutants[0]))
  ['acbA']
>>> # List a set of all phenotypes referenced by a single mutant
>>> print(mutant_phenotypes(dicty_mutants[0]))
  ['decreased lipid binding']

Function References
===================

.. automodule:: orangecontrib.bioinformatics.dicty.phenotypes
   :members:
   :member-order: bysource


