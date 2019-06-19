#!/usr/bin/env bash

if [ ${ORANGE_INSTALL} == "source" ];
then
    echo "Orange install for source"
    pip install --no-deps https://github.com/PrimozGodec/orange3/archive/test-tsne-import.zip
else
    echo "Orange install for conda"
    conda install orange3

fi
