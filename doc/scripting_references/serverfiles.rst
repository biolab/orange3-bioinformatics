===========
ServerFiles
===========

Utility that accesses files on a HTTP server and stores them locally for reuse.


Installation
============

::

	pip install serverfiles


Usage
=====

.. code:: ipython3

	"""
	Orange datasets are hosted on http://datasets.biolab.si/

	Folder structure:

	"core" -> publicly available datasets in Orange
	"sc "  -> publicly available datasets in scOrange

	To read and download datasets we use "https://github.com/biolab/serverfiles"
	installable trough "pip install serverfiles"

	"""

	import serverfiles
	from Orange.data import Table

.. code:: ipython3

	single_cell_url = 'http://datasets.biolab.si/sc'

	# connect to the HTTP file server
	sf = serverfiles.ServerFiles(server=single_cell_url)

	# list all available files
	sf.listfiles()


.. parsed-literal::

	[('DC_expMatrix_DCnMono.tab.gz',),
	 ('DC_expMatrix_deeper.characterization.tab.gz',),
	 ('aml-1k.pickle',),
	 ('aml-8k.pickle',),
	 ('ccp_data_Tcells_normCounts.counts.all_genes.tab.gz',),
	 ('ccp_data_Tcells_normCounts.counts.cycle_genes.tab.gz',),
	 ('ccp_data_liver.counts.all_genes.tab.gz',),
	 ('ccp_data_liver.counts.cycle_genes.tab.gz',),
	 ('ccp_data_mESCbulk.counts.all_genes.tab.gz',),
	 ('ccp_data_mESCbulk.counts.cycle_genes.tab.gz',),
	 ('ccp_normCountsBuettnerEtAl.counts.all_genes.tab.gz',),
	 ('ccp_normCountsBuettnerEtAl.counts.cycle_genes.tab.gz',),
	 ('ccp_normCounts_mESCquartz.counts.all_genes.tab.gz',),
	 ('ccp_normCounts_mESCquartz.counts.cycle_genes.tab.gz',),
	 ('cdp_expression_macosko.tab.gz',),
	 ('cdp_expression_shekhar.tab.gz',),
	 ('dm_proj_neurons_li2017.pkl.gz',),
	 ('nestorawa_forcellcycle.pkl.gz',),
	 ('pbmc_kang2018_raw_control.pkl.gz',),
	 ('pbmc_kang2018_raw_stimulated.pkl.gz',),
	 ('pbmc_kang2018_sample.pkl.gz',)]


.. code:: ipython3

	# setup local folder structure
	lf = serverfiles.LocalFiles('path_to_my_data', serverfiles=sf)
	lf.listfiles()


.. parsed-literal::

	[]


.. code:: ipython3

	# Download data from serverfiles to localfiles
	lf.download('aml-1k.pickle')
	lf.listfiles()

	# to get the path of downloaded file
	# lf.localpath('aml-1k.pickle')


.. parsed-literal::

	[('aml-1k.pickle',)]


.. code:: ipython3

	# to get a path of a given file use localpath_download,
	# if the file does not exist it will be downloaded from the
	# serverfiles automatically

	aml_1k = lf.localpath_download('aml-1k.pickle')
	print(aml_1k)
	print(lf.listfiles())
	aml_8k = lf.localpath_download('aml-8k.pickle')
	print(aml_8k)
	print(lf.listfiles())


.. parsed-literal::

	path_to_my_data/aml-1k.pickle
	[('aml-1k.pickle',)]
	path_to_my_data/aml-8k.pickle
	[('aml-8k.pickle',), ('aml-1k.pickle',)]


.. code:: ipython3

	# This is useful because now we can do something like this:
	def load(experiment):
	# checks if new version available
	lf.update(experiment)
	# ensures that we have the data localy
	return Table(lf.localpath_download(experiment))

	# note that using load method we keep our localfiles and serverfiles in sync.
	# it will always get the latest version from the serverfiles.
	load('aml-1k.pickle')


.. parsed-literal::

	[[0.000, 0.000, 5.648, 5.267, 0.000, ... | healthy] {1, 6681b0788fc2b2e21975e26f588cc7a9, ACGGGAGATGACCA-1},
	 [0.983, 0.573, 0.000, 0.000, 0.000, ... | AML] {1, 7c1e27874a82e7d5c1c877fa2cba7ba7, GAACAGCTGCTTAG-1},
	 [0.000, 0.000, 0.000, 0.000, 4.069, ... | healthy] {2, b9d78fbc8b1bf9aa478e1fbd93ac883c, TCTTACGAAAAAGC-1},
	 [0.000, 0.000, 0.000, 0.000, 2.304, ... | healthy] {1, 841c0e79f017a8ece40b2532f82a9c7a, CGAGCCGACTCTAT-1},
	 [0.000, 0.000, 0.000, 0.000, 0.000, ... | healthy] {2, 08c51bbdc6bb95fc17cd3965e6d6c4fb, GCAAACTGTTGGCA-1},
	 ...
	]


.. code:: ipython3

	# if one does not want to work with LocalFiles
	# access data directly trough ServerFiles

	sf.download('aml-1k.pickle', target='./my_downloaded_file.pickle')


.. figure:: https://user-images.githubusercontent.com/15876321/53346326-16bc9700-3917-11e9-835d-00a5407310ea.png

