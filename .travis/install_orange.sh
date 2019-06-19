#!/usr/bin/env bash

if [ ${ORANGE_INSTALL} == "source" ];
then
    echo "Orange install for source"
    pip install --no-deps https://github.com/biolab/orange3/archive/0b4ebacf78db6d7766a1128c1ea06f9d30a17c42.zip
else
    echo "Orange install for conda"
    conda install orange3

fi
