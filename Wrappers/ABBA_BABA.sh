#! /usr/bin/env bash

set -e
set -o pipefail

#   Load variables from supplied config file
source $1

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
OUT=${SCRATCH}/${PROJECT}/ABBABABA
mkdir -p ${OUT}

#   Check for local R installation
if `command -v Rscript > /dev/null 2> /dev/null`
then
    echo "R is insalled"
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
        -P "${N_CORES}" \
        -checkBamHeaders "${CHECK_BAM_HEADERS}" \
        -out "${OUT}"/"${PROJECT}".D
#   Assuming a single reigon was defined in config file
else
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
        -P "${N_CORES}" \
        -checkBamHeaders "${CHECK_BAM_HEADERS}" \
        -r "${REGIONS}" \
        -out "${OUT}"/"${PROJECT}".D
fi

#   jackKnife.R is provided with angsd.
Rscript "${ANGSD_DIR}"/R/jackKnife.R \
    file="${OUT}"/"${PROJECT}".D.abbababa \
    indNames="${SAMPLE_LIST}" \
    outfile="${OUT}"/"${PROJECT}".abbababa
