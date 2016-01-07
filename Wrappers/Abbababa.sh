#! /usr/bin/env bash

set -e
set -o pipefail

#   Load variables from supplied config file
source "$1"

#   Are we using Common_Config? If so, source it
if [[ -f "${COMMON}" ]]
then
    source "${COMMON}"
fi

#   Where is angsd-wrapper located?
SOURCE=$2

#   Where is ANGSD located?
ANGSD_DIR=${SOURCE}/dependencies/angsd

#   Create outdirectory
OUT=${SCRATCH}/${PROJECT}/Abbababa
mkdir -p "${OUT}"

#   Check for local R installation
if $(command -v Rscript > /dev/null 2> /dev/null)
then
    echo "R is installed"
else
    echo "R is not installed or not in your PATH"
    exit 1
fi

#   What's our outgroup?
if [[ -z "${OUTGROUP}" ]]
then
    OUTGROUP="${ANC_SEQ}"
fi

#   Now we actually run the command
#   Do we have a regions file?
if [[ -f "${REGIONS}" ]]
then
    if ! `command -v python > /dev/null 2> /dev/null`; then echo "Please install Python and place in your PATH" >&2; exit 1; fi
    echo "Sorting ${REGIONS} to match the order in ${OUTGROUP}" >&2
    FAI=$(ls $( dirname "${OUTGROUP}" )| grep -E "$( basename ${OUTGROUP} )\.fai|$( basename ${OUTGROUP} .fasta )\.fai")
    python "${SOURCE}"/Wrappers/sortRegions.py --fai $(dirname "${OUTGROUP}")/"${FAI}" --regions "${REGIONS}" --project "${PROJECT}"
    REGIONS=$(dirname "${REGIONS}")/"${PROJECT}"_SortedRegions.txt
    echo "Running Abbababa" >&2
    "${ANGSD_DIR}"/angsd \
        -doAbbababa "${DO_ABBABABA}" \
        -rmTrans "${REMOVE_TRANS}" \
        -blockSize "${BLOCKSIZE}" \
        -doCounts "${DO_COUNTS}" \
        -anc "${OUTGROUP}" \
        -bam "${SAMPLE_LIST}" \
        -uniqueOnly "${UNIQUE_ONLY}" \
        -minMapQ "${MIN_MAPQ}" \
        -minQ "${MIN_BASEQUAL}" \
        -minInd "${MIN_IND}" \
        -nThreads "${N_CORES}" \
        -checkBamHeaders "${CHECK_BAM_HEADERS}" \
        -rf "${REGIONS}" \
        -out "${OUT}"/"${PROJECT}".D
#   Are we missing a definiton for regions?
elif [[ -z "${REGIONS}" ]]
then
    echo "Running Abbababa" >&2
    "${ANGSD_DIR}"/angsd \
        -doAbbababa "${DO_ABBABABA}" \
        -rmTrans "${REMOVE_TRANS}" \
        -blockSize "${BLOCKSIZE}" \
        -doCounts "${DO_COUNTS}" \
        -anc "${OUTGROUP}" \
        -bam "${SAMPLE_LIST}" \
        -uniqueOnly "${UNIQUE_ONLY}" \
        -minMapQ "${MIN_MAPQ}" \
        -minQ "${MIN_BASEQUAL}" \
        -minInd "${MIN_IND}" \
        -nThreads "${N_CORES}" \
        -checkBamHeaders "${CHECK_BAM_HEADERS}" \
        -out "${OUT}"/"${PROJECT}".D
#   Assuming a single reigon was defined in config file
else
    echo "Running Abbababa" >&2
    "${ANGSD_DIR}"/angsd \
        -doAbbababa "${DO_ABBABABA}" \
        -rmTrans "${REMOVE_TRANS}" \
        -blockSize "${BLOCKSIZE}" \
        -doCounts "${DO_COUNTS}" \
        -anc "${OUTGROUP}" \
        -bam "${SAMPLE_LIST}" \
        -uniqueOnly "${UNIQUE_ONLY}" \
        -minMapQ "${MIN_MAPQ}" \
        -minQ "${MIN_BASEQUAL}" \
        -minInd "${MIN_IND}" \
        -nThreads "${N_CORES}" \
        -checkBamHeaders "${CHECK_BAM_HEADERS}" \
        -r "${REGIONS}" \
        -out "${OUT}"/"${PROJECT}".D
fi

#   jackKnife.R is provided with angsd.
echo "Using jackKnife.R to finish Abbababa" >&2
Rscript "${ANGSD_DIR}"/R/jackKnife.R \
    file="${OUT}"/"${PROJECT}".D.abbababa \
    indNames="${SAMPLE_LIST}" \
    outfile="${OUT}"/"${PROJECT}".abbababa
