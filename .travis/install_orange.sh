#!/usr/bin/env bash

if [ ${ORANGE_INSTALL} == "source" ];
then
    echo "Orange install for source"
    pip install --no-deps https://github.com/biolab/orange3/archive/514ae98b1409a38a9176d53d3fae12122b749668.zip
else
    echo "Orange install for conda"
    conda install orange3

fi
