#!/usr/bin/env bash

if [ $ORANGE_INSTALL == "source" ];
then
    echo "Orange install for source"
    pip install https://github.com/biolab/orange3/archive/master.zip
else
    echo "Orange install for conda"
    conda install orange3

fi