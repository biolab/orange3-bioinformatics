Orange3-bioinformatics
=======================

<a href='https://orange3-bioinformatics.readthedocs.io/en/latest/?badge=latest'>
    <img src='https://readthedocs.org/projects/orange3-bioinformatics/badge/?version=latest' alt='Documentation Status' />
</a>

<a href="https://codecov.io/gh/biolab/orange3-bioinformatics">
  <img src="https://codecov.io/gh/biolab/orange3-bioinformatics/branch/master/graph/badge.svg" />
</a>

<a href="https://pypi.org/project/Orange3-Bioinformatics/">
  <img alt="PyPI" src="https://img.shields.io/pypi/v/orange3-bioinformatics.svg" />
</a>

<a href="https://anaconda.org/conda-forge/orange3-bioinformatics">
  <img alt="Conda" src="https://img.shields.io/conda/v/conda-forge/orange3-bioinformatics.svg" />
</a>

<a href="https://pypi.org/project/Orange3-Bioinformatics/">
  <img alt="PyPI - License" src="https://img.shields.io/pypi/l/orange3-bioinformatics.svg" />
</a>


Orange Bioinformatics extends Orange, a data mining software
package, with common functionality for bioinformatics. The provided
functionality can be accessed as a Python library or through a visual
programming interface (Orange Canvas). The latter is also suitable for
non-programmers.

In Orange Canvas the analyst connects basic computational units, called
widgets, into data flow analytics schemas. Two units-widgets can be
connected if they share a data type. Compared to other popular tools like
Taverna, Orange widgets are high-level, integrated potentially complex
tasks, but are specific enough to be used independently. Even elaborate
analyses rarely consist of more than ten widgets; while tasks such as
clustering and enrichment analysis could be executed with up to five
widgets. While building the schema each widget is independently controlled
with settings, the settings do not conceptually burden the analyst.

Orange Bioinformatics provides access to publicly available data, like GEO data sets, GO and KEGG.
All features can be combined with powerful visualization, network exploration and
data mining techniques from the Orange data mining framework.

Installation
------------

To install the add-on with pip use

    pip install Orange3-bioinformatics

To register this add-on with Orange, but keep the code in the development directory (do not copy it to
Python's site-packages directory), run

    pip install -e .

Documentation / widget help can be built by running

    make html htmlhelp

from the doc directory.

Usage
-----

After the installation, the widgets from this add-on are registered with Orange. To run Orange from the terminal, use

    python3 -m Orange.canvas
    
or

    orange-canvas

The new widgets are in the toolbox bar under Bioinformatics section.
