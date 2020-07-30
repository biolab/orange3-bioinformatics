#!/usr/bin/env python
import setuptools

if __name__ == '__main__':
    setuptools.setup(
        use_scm_version=True,
        setup_requires=['setuptools-scm', 'setuptools>=40.0'],
        install_requires=[
            'serverfiles==0.3.0',
            'requests==2.22.0',
            'requests-cache==0.5.2',
            'genesis-pyapi==1.2.1',
            'pyclipper==1.1.0.post1',
            'Orange3>=3.22.0',
            'point-annotator~=2.0',
            'resdk~=12.1.1',
            'scipy>=1.5.0',
            # Versions are determined by Orange
            'numpy',
        ],
        extras_require={
            'doc': ['sphinx', 'recommonmark'],
            'test': [
                'flake8~=3.8.3',
                'flake8-comprehensions~=3.2.3',
                'flake8-black~=0.2.0',
                'pep8-naming~=0.11.1',
                'isort~=4.3.21',
                'pre-commit~=2.5.1',
                'pytest~=5.4.3',
                'coverage~=5.1',
                'codecov~=2.1.7',
            ],
        },
    )
