#!/usr/bin/env bash

if [ ${ORANGE_INSTALL} == "source" ];
then
    echo "Orange install for source"
    pip install --no-deps https://github.com/biolab/orange3/archive/c57f851170f01c724a7f71df4d3b35364fb451e6.zip
else
    echo "Orange install for conda"
    conda install orange3

fi
