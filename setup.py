#!/usr/bin/env python
import setuptools

if __name__ == '__main__':
    setuptools.setup(
        use_scm_version=True,
        setup_requires=['setuptools-scm', 'setuptools>=40.0'],
        install_requires=[
            'Orange3>=3.28.0',
            'orange-widget-base>=4.14.1',
            'scipy>=1.5.0',
            'pyclipper>=1.2.0',
            'point-annotator~=2.0',
            'requests',
            'requests-cache>=0.8.0',
            'serverfiles',
            'resdk>=13.3.0',
            'genesis-pyapi',
            # Versions are determined by Orange
            'numpy',
        ],
        extras_require={
            'doc': ['sphinx', 'recommonmark', 'sphinx_rtd_theme'],
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
