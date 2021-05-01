#!/usr/bin/env bash

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

#   Where is ANGSD?
ANGSD_DIR="${SOURCE}"/dependencies/angsd
NGS_POPGEN_DIR="${SOURCE}"/dependencies/ngsPopGen

#   Make the outdirectory
OUT=${SCRATCH}/${PROJECT}/PCA
mkdir -p "${OUT}"

#   Find number of individuals
N_IND=$(wc -l < "${SAMPLE_LIST}")

#   Do we have a regions file?
if [[ -f "${REGIONS}" ]]
then
    WRAPPER_ARGS=$(echo -bam "${SAMPLE_LIST}" \
        -GL "${GT_LIKELIHOOD}" \
        -out "${OUT}"/"${PROJECT}"_PCA \
        -doMajorMinor "${DO_MAJORMINOR}" \
        -doMaf "${DO_MAF}" \
        -doGeno "${DO_GENO}"\
        -doPost "${DO_POST}" \
        -nInd "${N_IND}" \
        -P "${N_CORES}" \
        -rf "${REGIONS}")
#   Are we missing a definiton for regions?
elif [[ -z "${REGIONS}" ]]
then
    WRAPPER_ARGS=$(echo -bam "${SAMPLE_LIST}" \
        -GL "${GT_LIKELIHOOD}" \
        -out "${OUT}"/"${PROJECT}"_PCA \
        -doMajorMinor "${DO_MAJORMINOR}" \
        -doMaf "${DO_MAF}" \
        -doGeno "${DO_GENO}"\
        -doPost "${DO_POST}" \
        -nInd "${N_IND}" \
        -P "${N_CORES}")
#   Assuming a single region was defined in config file
else
    WRAPPER_ARGS=$(echo -bam "${SAMPLE_LIST}" \
        -GL "${GT_LIKELIHOOD}" \
        -out "${OUT}"/"${PROJECT}"_PCA \
        -doMajorMinor "${DO_MAJORMINOR}" \
        -doMaf "${DO_MAF}" \
        -doGeno "${DO_GENO}" \
        -doPost "${DO_POST}" \
        -nInd "${N_IND}" \
        -P "${N_CORES}" \
        -r "${REGIONS}")
fi

# Check for advanced arguments, and overwrite any overlapping definitions
FINAL_ARGS=( $(source ${SOURCE}/Wrappers/Arg_Zipper.sh "${WRAPPER_ARGS}" "${ADVANCED_ARGS}") )
# echo "Final arguments: ${FINAL_ARGS}" 1<&2
"${ANGSD_DIR}"/angsd "${FINAL_ARGS[@]}"
#
gunzip "${OUT}/${PROJECT}_PCA.geno.gz" -f

"${NGS_POPGEN_DIR}"/ngsCovar -probfile "${OUT}/${PROJECT}"_PCA.geno \
    -outfile "${OUT}/${PROJECT}"_PCA.graph.me \
    -nind "${N_IND}" \
    -nsites "${N_SITES}" \
    -norm "${NORM}" \
    -call "${CALL}"
