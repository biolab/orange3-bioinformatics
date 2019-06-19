#!/usr/bin/env bash

if [ ${ORANGE_INSTALL} == "source" ];
then
    echo "Orange install for source"
    pip install --no-deps https://github.com/biolab/orange3/archive/bc61f182f155de198f43de018d576afb88265a54.zip
else
    echo "Orange install for conda"
    conda install orange3

fi
