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
echo "${SAMPLE_LIST}" 1<&2
N_IND=$(wc -l < "${SAMPLE_LIST}" | tr -d '[[:space:]]')

#       How many inbreeding coefficients are supplied?
N_F=$(wc -l < "${SAMPLE_INBREEDING}")
#       For ANGSD, the actual sample size is twice the number of individuals, since each individual has two chromosomes.
#       The individual inbreeding coefficents take care of the mismatch between these two numbers
N_CHROM=$(expr 2 \* "${N_IND}")

#   Perform a check to see if number of individuals matches number of inbreeding coefficients
if [ "${N_IND}" -ne "${N_F}" ]
then
    echo "WRAPPER: Mismatch between number of samples in ${SAMPLE_LIST} and ${SAMPLE_INBREEDING}"
    exit 1
fi

#   Create outdirectory
OUT=${SCRATCH}/${PROJECT}/Thetas
mkdir -p "${OUT}"

if [[  -f "${OUT}"_Diversity.mafs.gz ]] && [ "${OVERRIDE}" = "false" ]; then
    echo "WRAPPER: maf already exists and OVERRIDE=false, skipping angsd -bam...";
else
    #   Now we actually run the command, this creates a binary file that contains the prior SFS
    #   Do we have a regions file?
    if [[ -f "${REGIONS}" ]]
    then
	WRAPPER_ARGS=$(echo -bam "${SAMPLE_LIST}" \
            -out "${OUT}"/"${PROJECT}"_Diversity \
            -indF "${SAMPLE_INBREEDING}" \
            -doSaf "${DO_SAF}" \
            -doThetas 1 \
            -uniqueOnly "${UNIQUE_ONLY}" \
            -anc "${ANC_SEQ}" \
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
            -pest "${PEST}" \
            -rf "${REGIONS}")
    #   Are we missing a definiton for regions?
    elif [[ -z "${REGIONS}" ]]
    then
	WRAPPER_ARGS=$(echo -bam "${SAMPLE_LIST}" \
            -out "${OUT}"/"${PROJECT}"_Diversity \
            -indF "${SAMPLE_INBREEDING}" \
            -doSaf "${DO_SAF}" \
            -doThetas 1 \
            -uniqueOnly "${UNIQUE_ONLY}" \
            -anc "${ANC_SEQ}" \
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
            -pest "${PEST}")
    #   Assuming a single region was defined in config file
    else
	WRAPPER_ARGS=$(echo -bam "${SAMPLE_LIST}" \
        -out "${OUT}"/"${PROJECT}"_Diversity \
        -indF "${SAMPLE_INBREEDING}" \
        -doSaf "${DO_SAF}" \
        -doThetas 1 \
        -uniqueOnly "${UNIQUE_ONLY}" \
        -anc "${ANC_SEQ}" \
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
        -pest "${PEST}" \
        -r "${REGIONS}")
    fi
fi
# Check for advanced arguments, and overwrite any overlapping definitions
FINAL_ARGS=$(source ${SOURCE}/Wrappers/Arg_Zipper.sh "${WRAPPER_ARGS}" "${ADVANCED_ARGS}")
# echo "Final arguments: ${FINAL_ARGS}" 1<&2
"${ANGSD_DIR}"/angsd ${FINAL_ARGS}


# Either whole chromosome stats, or create many windowed stats
if [ "${SLIDING_WINDOW}" = "false" ]
then
    "${ANGSD_DIR}"/misc/thetaStat do_stat \
        "${OUT}"/"${PROJECT}"_Diversity.thetas.idx \
        -nChr "${N_CHROM}"
else
    "${ANGSD_DIR}"/misc/thetaStat do_stat \
        "${OUT}"/"${PROJECT}"_Diversity.thetas.idx \
        -nChr "${N_CHROM}" \
        -win "${WIN}" \
        -step "${STEP}"
fi

# Decoding the thetas.idx file, storing output into a text file
${SOURCE}/dependencies/angsd/misc/thetaStat print \
 	 ${OUT}/${PROJECT}_Diversity.thetas.idx > \
	 ${OUT}/${PROJECT}_Diversity.thetas.txt

# Filter pestPG file for invariant sites
echo "WRAPPER: Creating files for Shiny graphing..." >&2
Rscript ${SOURCE}/Wrappers/thetas_filtering.R \
    ${SOURCE} \
    ${OUT}/${PROJECT}_Diversity.thetas.idx.pestPG \
    ${OUT}/"${PROJECT}"_Thetas.graph.me
