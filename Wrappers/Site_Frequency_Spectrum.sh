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
ANGSD_DIR=${SOURCE}/dependencies/angsd

#   Variables created from transforming other variables
#       The number of individuals in the taxon we are analyzing
N_IND=$(wc -l < "${SAMPLE_LIST}")
#       How many inbreeding coefficients are supplied?
N_F=$(wc -l < "${SAMPLE_INBREEDING}")
#       For ANGSD, the actual sample size is twice the number of individuals, since each individual has two chromosomes.
#       The individual inbreeding coefficents take care of the mismatch between these two numbers

#   Perform a check to see if number of individuals matches number of inbreeding coefficients
if [ "${N_IND}" -ne "${N_F}" ]
then
    echo "Mismatch between number of samples in ${SAMPLE_LIST} and ${SAMPLE_INBREEDING}"
    exit 1
fi

#   Check to see if ancestral state is supplied: If not, polarize samples using
#   the reference sequence and generate folded saf.
if [ ! -f "${ANC_SEQ}" ]
then
    echo "Ancestral state data not found, using reference sequence to polarize alignment data. BAQ will likewise not be calculated."
    if [ ! -f "${REF_SEQ}" ]
    then
        echo "No reference sequence supplied, unable to perform calculations."
        exit 2
    else
        ANC_SEQ=$REF_SEQ
        REF_SEQ=
        BAQ=0  # ANGSD can't calculate BAQ without a supplied reference
        FOLD=1
    fi
else
    FOLD=0
fi

#   Create outdirectory
OUT="${SCRATCH}"/"${PROJECT}"/SFS
mkdir -p "${OUT}"

#   Now we actually run the command, this creates a binary file that contains the prior SFS
if [[ -f "${OUT}"/"${PROJECT}"_SFSOut.mafs.gz ]] && [ "$OVERRIDE" = "false" ]
then
    echo "WRAPPER:maf already exists and OVERRIDE=false, skipping angsd -bam..."
else
    #   Do we have a regions file?
    if [[ -f "${REGIONS}" ]]
    then
	WRAPPER_ARGS=$(echo -bam "${SAMPLE_LIST}" \
            -out "${OUT}"/"${PROJECT}"_SFSOut \
            -indF "${SAMPLE_INBREEDING}" \
            -doSaf "${DO_SAF}" \
            -uniqueOnly "${UNIQUE_ONLY}" \
            -anc "${ANC_SEQ}" \
            -minMapQ "${MIN_MAPQ}" \
            -minQ "${MIN_BASEQUAL}" \
            -nInd "${N_IND}" \
            -minInd "${MIN_IND}"\
            -baq "${BAQ}" \
            -ref "${REF_SEQ}" \
            -GL "${GT_LIKELIHOOD}" \
            -P "${N_CORES}" \
            -doMajorMinor "${DO_MAJORMINOR}" \
            -doMaf "${DO_MAF}" \
            -doGeno "${DO_GENO}" \
            -rf "${REGIONS}" \
            -doPost "${DO_POST}")
    #   Are we missing a definiton for regions?
    elif [[ -z "${REGIONS}" ]]
    then
	WRAPPER_ARGS=$(echo -bam "${SAMPLE_LIST}" \
            -out "${OUT}"/"${PROJECT}"_SFSOut \
            -indF "${SAMPLE_INBREEDING}" \
            -doSaf "${DO_SAF}" \
            -uniqueOnly "${UNIQUE_ONLY}" \
            -anc "${ANC_SEQ}" \
            -minMapQ "${MIN_MAPQ}" \
            -minQ "${MIN_BASEQUAL}" \
            -nInd "${N_IND}" \
            -minInd "${MIN_IND}"\
            -baq "${BAQ}" \
            -ref "${REF_SEQ}" \
            -GL "${GT_LIKELIHOOD}" \
            -P "${N_CORES}" \
            -doMajorMinor "${DO_MAJORMINOR}" \
            -doMaf "${DO_MAF}" \
            -doGeno "${DO_GENO}" \
            -doPost "${DO_POST}")
    #   Assuming a single region was defined in config file
    else
	WRAPPER_ARGS=$(echo -bam "${SAMPLE_LIST}" \
            -out "${OUT}"/"${PROJECT}"_SFSOut \
            -indF "${SAMPLE_INBREEDING}" \
            -doSaf "${DO_SAF}" \
            -uniqueOnly "${UNIQUE_ONLY}" \
            -anc "${ANC_SEQ}" \
            -folded "${FOLD}" \
            -minMapQ "${MIN_MAPQ}" \
            -minQ "${MIN_BASEQUAL}" \
            -nInd "${N_IND}" \
            -minInd "${MIN_IND}" \
            -baq "${BAQ}" \
            -ref "${REF_SEQ}" \
            -GL "${GT_LIKELIHOOD}" \
            -P "${N_CORES}" \
            -doMajorMinor "${DO_MAJORMINOR}" \
            -doMaf "${DO_MAF}" \
            -doGeno "${DO_GENO}" \
            -doPost "${DO_POST}" \
            -r "${REGIONS}")
    fi
fi
# Check for advanced arguments, and overwrite any overlapping definitions
# Creates a single argument array from both inputs
FINAL_ARGS=( "$(bash "${SOURCE}/Wrappers/Arg_Zipper.sh" "${WRAPPER_ARGS}" "${ADVANCED_ARGS}")" )
# DEBUGGING
# echo "Wrapper arguments: ${WRAPPER_ARGS}" 1<&2
# echo -e "Final arguments:" ${FINAL_ARGS} 1<&2

"${ANGSD_DIR}"/angsd ${FINAL_ARGS[@]}

"${ANGSD_DIR}"/misc/realSFS \
    "${OUT}/${PROJECT}"_SFSOut.saf.idx \
    -P "${N_CORES}" \
    -fold "${FOLD}" \
    > "${OUT}/${PROJECT}"_DerivedSFS.graph.me
