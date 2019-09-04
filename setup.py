#!/usr/bin/env python
import setuptools

if __name__ == '__main__':
    setuptools.setup(
        use_scm_version=True,
        setup_requires=[
            'setuptools-scm',
            'setuptools>=40.0'
        ],
        install_requires=[
            'serverfiles==0.3.0',
            'requests==2.22.0',
            'requests-cache==0.5.2',
            'genesis-pyapi==1.2.1',
            'pyclipper==1.1.0.post1',
            'Orange3>=3.22.0',
            'point-annotator~=1.0',
            # Versions are determined by Orange
            'scipy',
            'numpy',
        ],
        extras_require={
            'doc': ['sphinx', 'recommonmark'],
            'test': [
                'flake8~=3.7.8',
                'flake8-comprehensions~=2.2.0',
                'flake8-black~=0.1.0',
                'pep8-naming~=0.8.2',
                'isort~=4.3.21',
                'pre-commit~=1.18.3',
                'pytest~=5.1.0',
                'pytest-cov~=2.7.1',
                'coverage~=4.5.4',
                'codecov~=2.0.15'
            ],
        },
    )
