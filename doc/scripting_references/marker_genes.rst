============
Marker Genes
============

We provide access to two different databases (of marker genes):

- `PanglaoDB <https://panglaodb.se/>`_
- `CellMarker <http://biocc.hrbmu.edu.cn/CellMarker/>`_

Data is preprocessed in Orange readable format and it is hosted `here <http://download.biolab.si/datasets/bioinformatics/marker_genes/>`_.

To get marker genes we need access files on our http file server. We can use a `tool <https://github.com/biolab/serverfiles>`_ that handles just that.
Detailed explanation on how to use it can be found  :doc:`here. <serverfiles>`

Usage
=====

>>> from Orange.data import Table
>>> from orangecontrib.bioinformatics.utils import serverfiles
>>> markers = serverfiles.localpath_download('marker_genes', 'panglao_gene_markers.tab')
>>> data = Table(markers)
>>> data
>>>
[[] {Mouse, CTRB1, Acinar cells, PanglaoDB, https://panglaodb.se/, ...},
 [] {Mouse, KLK1, Acinar cells, PanglaoDB, https://panglaodb.se/, ...},
 [] {Mouse, RBPJL, Acinar cells, PanglaoDB, https://panglaodb.se/, ...},
 [] {Mouse, PTF1A, Acinar cells, PanglaoDB, https://panglaodb.se/, ...},
 [] {Mouse, TRY4, Acinar cells, PanglaoDB, https://panglaodb.se/, ...},
 ...
]