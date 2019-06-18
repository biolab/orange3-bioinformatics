#!/usr/bin/env bash

if [ ${ORANGE_INSTALL} == "source" ];
then
    echo "Orange install for source"
    pip install --no-deps https://github.com/biolab/orange3/archive/3.21.0.zip
else
    echo "Orange install for conda"
    conda install orange3

fi
