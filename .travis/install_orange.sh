#!/usr/bin/env bash

if [ ${ORANGE_INSTALL} == "source" ];
then
    echo "Orange install for source"
    pip install --no-deps https://github.com/biolab/orange3/archive/d672ab78dee0449c79a942c2d8ff2dbe38956394.zip
else
    echo "Orange install for conda"
    conda install orange3

fi
