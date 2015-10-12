#!/bin/bash

set -e
set -o pipefail

#   Load the variables from supplied config file
source $1

#   Are we using Common_Config? If so, source it
if [[ -f "${COMMON}" ]]
then
    source "${COMMON}"
fi

#   Where is angsd-wrapper located?
SOURCE=$2

#   Create outdirectory
OUT=${SCRATCH}/${PROJECT}/plots
mkdir -p ${OUT}

#   Where are our graphing scripts?
GRAPH_DIR=${SOURCE}/Rplots

#   Check for local R installation
if ! `command -v Rscript > /dev/null 2> /dev/null`
then
    echo "R is neither installed nor in your PATH!"
    exit 1
fi

#   Generate plots for each variable supplied

#   Site Frequency Spectrum
if ! [[ -z "${SFS_OUT}" ]]
then
    Rscript "${GRAPH_DIR}"/SFS_plot.R "${SFS_OUT}"
fi
