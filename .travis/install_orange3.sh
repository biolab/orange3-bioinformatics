#!/usr/bin/env bash

cd $TRAVIS_BUILD_DIR/orange

# remove if exist
rm -f -R orange3
# clone orange from git
git clone https://github.com/biolab/orange3.git
cd orange3

# install orange
pip install -e .

# go back to add-on dir
cd $TRAVIS_BUILD_DIR