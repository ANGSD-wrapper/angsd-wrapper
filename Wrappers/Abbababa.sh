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
if command -v Rscript > /dev/null 2> /dev/null
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
    echo "Sorting ${REGIONS} to match the order in ${REF_SEQ}" >&2
    REF_EXT=$(echo "${REF_SEQ}" | rev | cut -f 1 -d '.' | rev)
    FAI=$(find "$(dirname "${REF_SEQ}")" -name "$(basename "$REF_SEQ" ".${REF_EXT}")*.fai")
    Rscript "${SOURCE}"/Wrappers/sortRegions.R "${REGIONS}" "${FAI}"
    REGIONS="$(find "$(dirname "$REGIONS")" -name '*_sorted.txt')"
        echo "Running Abbababa" >&2
    WRAPPER_ARGS=$(echo -doAbbababa "${DO_ABBABABA}" \
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
        -rf "${REGIONS}" \
        -out "${OUT}"/"${PROJECT}".D \
        -useLast "${USE_LAST}" \
        -enhance "${ENHANCE}")
#   Are we missing a definiton for regions?
elif [[ -z "${REGIONS}" ]]
then
    echo "Running Abbababa" >&2
    WRAPPER_ARGS=$(echo -doAbbababa "${DO_ABBABABA}" \
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
        -out "${OUT}"/"${PROJECT}".D \
        -useLast "${USE_LAST}" \
        -enhance "${ENHANCE}")
#   Assuming a single reigon was defined in config file
else
    echo "Running Abbababa" >&2
    WRAPPER_ARGS=$(echo -doAbbababa "${DO_ABBABABA}" \
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
        -r "${REGIONS}" \
        -out "${OUT}"/"${PROJECT}".D \
        -useLast "${USE_LAST}" \
        -enhance "${ENHANCE}")
fi
# Check for advanced arguments, and overwrite any overlapping definitions
FINAL_ARGS=$(source "${SOURCE}"/Wrappers/Arg_Zipper.sh "${WRAPPER_ARGS}" "${ADVANCED_ARGS}")
echo "Final arguments: ${FINAL_ARGS}" 1<&2
"${ANGSD_DIR}"/angsd ${FINAL_ARGS}

#   jackKnife.R is provided with angsd.
echo "Using jackKnife.R to finish Abbababa" >&2
Rscript "${ANGSD_DIR}"/R/jackKnife.R \
    file="${OUT}"/"${PROJECT}".D.abbababa \
    indNames="${SAMPLE_LIST}" \
    outfile="${OUT}"/"${PROJECT}".abbababa

#   Move ${PROJECT}.abbababa.txt to ${PROJECT}_Abbababa.graph.me
mv "${OUT}"/"${PROJECT}".abbababa.txt "${OUT}"/"${PROJECT}"_Abbababa.graph.me
