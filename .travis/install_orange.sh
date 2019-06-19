#!/usr/bin/env bash

if [ ${ORANGE_INSTALL} == "source" ];
then
    echo "Orange install for source"
    pip install --no-deps https://github.com/biolab/orange3/archive/55954d6ea2ee38b8bf4e835051b7ead29f3839b6.zip
else
    echo "Orange install for conda"
    conda install orange3

fi
