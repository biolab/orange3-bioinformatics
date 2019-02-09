#!/usr/bin/env bash

# https://conda.io/docs/user-guide/tasks/use-conda-with-travis-ci.html#the-travis-yml-file

# install miniconda
wget http://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh -O miniconda.sh
bash miniconda.sh -b -p $HOME/miniconda
export PATH="$HOME/miniconda/bin:$PATH"
hash -r
# update miniconda
conda config --set always_yes yes --set changeps1 no
conda config --add channels conda-forge
conda update -q conda
conda info -a