.. py:currentmodule:: orangecontrib.bioinformatics.geo.dataset

.. index:: NCBI
.. index:: GEO
.. index:: Gene Expression Omnibus
.. index:: microarray data sets

===============================
NCBI's Gene Expression Omnibus
===============================

This module provides an interface to `NCBI
<http://www.ncbi.nlm.nih.gov/>`_'s `Gene Expression Omnibus
<http://www.ncbi.nlm.nih.gov/geo/>`_ repository. It
supports `GEO DataSets <http://www.ncbi.nlm.nih.gov/sites/GDSbrowser>`_
query and retrieval.

In the following example :obj:`GDS.get_data`
construct a data set with genes in rows and samples in
columns. Notice that the annotation about each sample is retained
in ``.attributes``.

>>> from orangecontrib.bioinformatics.geo.dataset import GDS
>>> gds = GDS("GDS1676")
>>> data = gds.get_data()
>>> len(data)
719
>>> data[200]
[0.503, 0.690, 0.607, -2.250, 0.000, ...] {CD40}
>>> data.domain.attributes[0]
ContinuousVariable(name='GSM63816', number_of_decimals=3)
>>> data.domain.attributes[0].attributes
{'infection': 'acute', 'time': '1 d', 'dose': '20 U/ml IL-2'}

Class References
=================

.. autoclass:: GDSInfo()
   :members:
   :special-members: __init__



.. autoclass:: GDS()
   :members:


Usage
=====

The following script prints out information about a specific data
set. It does not download the data set, just uses the (local) GEO data
sets information file (:download:`dataset_info.py <code/geo/dataset_info.py>`).

.. literalinclude:: code/geo/dataset_info.py

The output of this script is::

    ID:
    GDS10
    Features:
    39114
    Genes:
    29942
    Organism:
    Mus musculus
    PubMed ID:
    11827943
    Sample types:
      tissue (spleen, thymus)
      disease state (diabetic, diabetic-resistant, nondiabetic)
      strain (NOD, Idd3, Idd5, Idd3+Idd5, Idd9, B10.H2g7, B10.H2g7 Idd3)

    Description:
    Examination of spleen and thymus of type 1 diabetes nonobese diabetic
    (NOD) mouse, four NOD-derived diabetes-resistant congenic strains and
    two nondiabetic control strains.


Samples in GEO data sets belong to sample subsets, which in turn belong
to specific types.  The above GDS10 has three sample types, of which the
subsets for the tissue type are spleen and thymus. For supervised data
mining it would be useful to find out which data sets provide enough
samples for each label. It is (semantically) convenient to perform
classification within sample subsets of the same type. The following
script goes through all data sets and finds those with enough
samples within each of the subsets for a specific type. The function
``valid`` determines which subset types (if any) satisfy our criteria
(:download:`dataset_samples.py <code/geo/dataset_samples.py>`).

.. literalinclude:: code/geo/dataset_samples.py

The requested number of samples, ``n=40``, seems to be a quite
a stringent criteria met - at the time of writing this -
by 40 data sets with 48 sample subsets. The output starts with::

    GDS1292
      tissue:raphe magnus/40, somatomotor cortex/43
    GDS1293
      tissue:raphe magnus/40, somatomotor cortex/41
    GDS1412
      protocol:no treatment/47, hormone replacement therapy/42
    GDS1490
      other:non-neural/50, neural/100
    GDS1611
      genotype/variation:wild type/48, upf1 null mutant/48
    GDS2373
      gender:male/82, female/48
    GDS2808
      protocol:training set/44, validation set/50

Let us now pick data set GDS2960 and see if we can predict the disease
state. We will use logistic regression, and within 10-fold cross
validation measure AUC, the area under ROC. AUC is the probability of
correctly distinguishing the two classes, (e.g., the disease and control).
From (:download:`predict_disease_state.py <code/geo/predict_disease_state.py>`):

.. literalinclude:: code/geo/predict_disease_state.py

The output of this script is::

    Samples: 101, Genes: 4069
    AUC = 0.996

The AUC for this data set is very high, indicating that using these gene
expression data it is almost trivial to separate the two classes.
