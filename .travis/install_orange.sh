#!/usr/bin/env bash

if [ ${ORANGE_INSTALL} == "source" ];
then
    echo "Orange install for source"
    pip install --no-deps https://github.com/biolab/orange3/archive/33d9a9233e8d82cfa4a624be876377adca6284cb.zip
else
    echo "Orange install for conda"
    conda install orange3

fi
