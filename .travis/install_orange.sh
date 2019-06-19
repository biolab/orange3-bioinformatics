#!/usr/bin/env bash

if [ ${ORANGE_INSTALL} == "source" ];
then
    echo "Orange install for source"
    pip install --no-deps https://github.com/biolab/orange3/archive/c851e6cee25ea76a4c56178db4aa7d9a02772bf7.zip
else
    echo "Orange install for conda"
    conda install orange3

fi
