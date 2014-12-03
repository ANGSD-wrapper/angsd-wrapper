#!/usr/bin/env bash
#   Script to create the directory structure for angsd-wrapper
source common.conf
#   test if our directory is writable (it should be!)
if [ -w `dirname ${PROJECT_DIR}` ]
    then mkdir -p ${RESULTS_DIR} ${TEMP_DIR} ${DATA_DIR}
    else echo "Directory is not writable!"; exit 1
fi
