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

#   Variables created from transforming other variables
#       The number of individuals in the taxon we are analyzing
N_IND=$(wc -l < "${SAMPLE_LIST}")
#       How many inbreeding coefficients are supplied?
N_F=$(wc -l < "${SAMPLE_INBREEDING}")

#   Perform a check to see if number of individuals matches number of inbreeding coefficients
if [ "${N_IND}" -ne "${N_F}" ]
then
    echo "Mismatch between number of samples in ${SAMPLE_LIST} and ${SAMPLE_INBREEDING}"
    exit 1
fi

#   Create outdirectory
OUT="${SCRATCH}"/"${PROJECT}"/GenotypeLikelihoods
mkdir -p "${OUT}"

#   Where is ANGSD?
ANGSD_DIR="${SOURCE}"/dependencies/angsd

#   Now we actually run the command
#   Do we have a regions file?
if [[ -f "${REGIONS}" ]]
then
    "${ANGSD_DIR}"/angsd \
        -bam "${SAMPLE_LIST}" \
        -out "${OUT}"/"${PROJECT}"_snps \
        -indF "${SAMPLE_INBREEDING}" \
        -doMajorMinor "${DO_MAJORMINOR}" \
        -uniqueOnly "${UNIQUE_ONLY}" \
        -minMapQ "${MIN_MAPQ}" \
        -minQ "${MIN_BASEQUAL}" \
        -GL "${GT_LIKELIHOOD}" \
        -doGlf "${DO_GLF}" \	
        -rf "${REGIONS}" \
        -doGeno "${DO_GENO}" \
        -doPost "${DO_POST}"\
        -postCutoff "${POST_CUTOFF}" \
        -doMaf "${DO_MAF}" \
        -SNP_pval "${SNP_PVAL}" \
        -nInd "${N_IND}" \
        -minInd "${MIN_IND}" \
        -P "${N_CORES}"
#   Are we missing a definiton for regions?
elif [[ -z "${REGIONS}" ]]
then
    "${ANGSD_DIR}"/angsd \
        -bam "${SAMPLE_LIST}" \
        -out "${OUT}"/"${PROJECT}"_snps \
        -indF "${SAMPLE_INBREEDING}" \
        -doMajorMinor "${DO_MAJORMINOR}" \
        -uniqueOnly "${UNIQUE_ONLY}" \
        -minMapQ "${MIN_MAPQ}" \
        -minQ "${MIN_BASEQUAL}" \
        -GL "${GT_LIKELIHOOD}" \
	-doGlf "${DO_GLF}" \
        -doGeno "${DO_GENO}" \
        -doPost "${DO_POST}"\
        -postCutoff "${POST_CUTOFF}" \
        -doMaf "${DO_MAF}" \
        -SNP_pval "${SNP_PVAL}" \
        -nInd "${N_IND}" \
        -minInd "${MIN_IND}" \
        -P "${N_CORES}"
#   Assuming a single reigon was defined in config file
else
    "${ANGSD_DIR}"/angsd \
        -bam "${SAMPLE_LIST}" \
        -out "${OUT}"/"${PROJECT}"_snps \
        -indF "${SAMPLE_INBREEDING}" \
        -doMajorMinor "${DO_MAJORMINOR}" \
        -uniqueOnly "${UNIQUE_ONLY}" \
        -minMapQ "${MIN_MAPQ}" \
        -minQ "${MIN_BASEQUAL}" \
        -GL "${GT_LIKELIHOOD}" \
	-doGlf "${DO_GLF}" \
        -r "${REGIONS}" \
        -doGeno "${DO_GENO}" \
        -doPost "${DO_POST}"\
        -postCutoff "${POST_CUTOFF}" \
        -doMaf "${DO_MAF}" \
        -SNP_pval "${SNP_PVAL}" \
        -nInd "${N_IND}" \
        -minInd "${MIN_IND}" \
        -P "${N_CORES}"
fi
