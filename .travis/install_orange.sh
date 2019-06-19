#!/usr/bin/env bash

if [ ${ORANGE_INSTALL} == "source" ];
then
    echo "Orange install for source"
    pip install --no-deps https://github.com/biolab/orange3/archive/4e240420d2349cf59f463d6cff8129ac2b5a970b.zip
else
    echo "Orange install for conda"
    conda install orange3

fi
