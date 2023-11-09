#!/usr/bin/env python
import setuptools

if __name__ == '__main__':
    setuptools.setup(
        use_scm_version=True,
        setup_requires=['setuptools-scm', 'setuptools>=40.0'],
        install_requires=[
            'Orange3>=3.34.0',
            'orange-widget-base>=4.19.0',
            'pyclipper>=1.2.0',
            'point-annotator~=2.0',
            'requests',
            'requests-cache>=0.8.0',
            'resdk>=20.0.0',
            'genesis-pyapi',
            'single_sample_gsea>=0.2.0',
            # Versions are determined by Orange
            'numpy',
            'scipy',
            'serverfiles',
        ],
        extras_require={
            # docutils changed html in 0.17; fixing to 0.16 until parser fixed
            # todo: remove docutils when parser fixed in widget-base and released
            'doc': ['sphinx', 'recommonmark', 'sphinx_rtd_theme', 'docutils<0.17'],
            'test': [
                'flake8',
                'flake8-comprehensions',
                'flake8-black',
                'pep8-naming',
                'isort',
                'pre-commit',
                'pytest',
                'coverage',
                'codecov',
            ],
        },
    )
