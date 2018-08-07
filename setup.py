#!/usr/bin/env python
import os


from setuptools import setup, find_packages


NAME = 'Orange3-Bioinformatics'
DOCUMENTATION_NAME = 'Orange Bioinformatics'
VERSION = '3.2.0'

DESCRIPTION = 'Orange Bioinformatics add-on for Orange data mining software package.'
LONG_DESCRIPTION = open(os.path.join(os.path.dirname(__file__), 'README.md')).read()
AUTHOR = 'Bioinformatics Laboratory, FRI UL'
AUTHOR_EMAIL = 'info@biolab.si'
URL = 'https://github.com/biolab/orange3-bioinformatics'
LICENSE = 'GPL3+'

KEYWORDS = (
    # [PyPi](https://pypi.python.org) packages with keyword "orange3 add-on"
    # can be installed using the Orange Add-on Manager
    'orange3 add-on'
)

CLASSIFIERS = (
    'License :: OSI Approved :: GNU General Public License v3 or later (GPLv3+)',
    'Topic :: Scientific/Engineering :: Artificial Intelligence',
    'Topic :: Scientific/Engineering :: Visualization',
    'Topic :: Scientific/Engineering :: Bio-Informatics',
    'Topic :: Software Development :: Libraries :: Python Modules',
    'Intended Audience :: Education',
    'Intended Audience :: Science/Research',
    'Intended Audience :: Developers',
    'Development Status :: 1 - Planning',
    'Programming Language :: Python :: 3 :: Only',
    'Operating System :: OS Independent'
)

PACKAGES = find_packages()
PACKAGE_DATA = {}

requirements = ['requirements.txt']

INSTALL_REQUIRES = sorted(set(
    line.partition('#')[0].strip()
    for file in (os.path.join(os.path.dirname(__file__), file)
                 for file in requirements)
    for line in open(file)
) - {''})


ENTRY_POINTS = {
    'orange3.addon': (
        'bioinformatics = orangecontrib.bioinformatics'
    ),
    'orange.widgets': (
        'Bioinformatics = orangecontrib.bioinformatics.widgets'
    ),
    'orange.canvas.help': (
        'html-index = orangecontrib.bioinformatics.widgets:WIDGET_HELP_PATH'
    )
}

NAMESPACE_PACAKGES = ["orangecontrib"]

TEST_SUITE = "orangecontrib.bioinformatics.tests.suite"

if __name__ == '__main__':
    setup(
        name=NAME,
        version=VERSION,
        description=DESCRIPTION,
        long_description=LONG_DESCRIPTION,
        author=AUTHOR,
        author_email=AUTHOR_EMAIL,
        url=URL,
        classifiers=CLASSIFIERS,
        license=LICENSE,
        keywords=KEYWORDS,
        packages=PACKAGES,
        package_data=PACKAGE_DATA,
        install_requires=INSTALL_REQUIRES,
        entry_points=ENTRY_POINTS,
        namespace_packages=NAMESPACE_PACAKGES,
        test_suite=TEST_SUITE,
        include_package_data=True,
        zip_safe=False,
    )
