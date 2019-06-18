#!/usr/bin/env bash

if [ ${ORANGE_INSTALL} == "source" ];
then
    echo "Orange install for source"
    pip install --no-deps https://github.com/biolab/orange3/archive/9f5ca5fae6a50047ca06163cde4011cd8ec93b7f.zip
else
    echo "Orange install for conda"
    conda install orange3

fi
