#!/usr/bin/env bash

if [ ${ORANGE_INSTALL} == "source" ];
then
    echo "Orange install for source"
    pip install --no-deps https://github.com/biolab/orange3/archive/95b8d5f845861b33aff4469dc123236699a4354b.zip
else
    echo "Orange install for conda"
    conda install orange3

fi
