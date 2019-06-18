#!/usr/bin/env bash

if [ ${ORANGE_INSTALL} == "source" ];
then
    echo "Orange install for source"
    pip install --no-deps https://github.com/biolab/orange3/archive/67fafb09609e8cf49ebc9388f4d158e06ed77639.zip
else
    echo "Orange install for conda"
    conda install orange3

fi
