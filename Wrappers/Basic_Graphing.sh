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
    echo "Plotting SFS..."
    SFS_NAME=`basename "${SFS_OUT}" _DerivedSFS`
    Rscript "${GRAPH_DIR}"/SFS_plot.R "${SFS_OUT}" "${OUT}" "${SFS_NAME}"
    echo "Thetas plot can be found at ${OUT}/${SFS_NAME}_SFS.pdf"
fi

#   Thetas
if ! [[ -z "${THETAS_OUT}" ]]
then
    echo "Plotting Thetas..."
    THETAS_NAME=`basename "${THETAS_OUT}" .pest.pg`
    Rscript "${GRAPH_DIR}"/Thetas_plot.R "${THETAS_OUT}" "${OUT}" "${THETAS_NAME}"
    echo "Thetas plot can be found at ${OUT}/${THETAS_NAME}_Thetas.pdf"
fi
