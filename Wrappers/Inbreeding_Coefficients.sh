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

#   Where is ngsF?
NGSF_DIR="${SOURCE}"/dependencies/ngsF

N_IND=$(wc -l < "${SAMPLE_LIST}")

#   Make the outdirectory
OUT=${SCRATCH}/${PROJECT}/Inbreeding_Coefficients
mkdir -p "${OUT}"

#   Can't use beagle format
if [[ "${DO_GLF}" == 2 ]]; then echo "Can't use beagle format, please specify a different value for DO_GLF!"; exit 1; fi

if [[ -f "${OUT}"/"${PROJECT}"_.mafs.gz ]] && [[ "${OVERRIDE}" = "false" ]]
then
    echo "mafs already exists and OVERRIDE=false, skipping angsd -bam..."
else
#   Do we have a regions file?
    if [[ -f "${REGIONS}" ]]
    then
        "${ANGSD_DIR}"/angsd \
            -bam "${SAMPLE_LIST}" \
            -rf "${REGIONS}" \
            -doGLF "${DO_GLF}" \
            -GL "${GT_LIKELIHOOD}" \
            -out "${OUT}"/"${PROJECT}" \
            -ref "${REF_SEQ}" \
            -anc "${ANC_SEQ}" \
            -doMaf "${DO_MAF}" \
            -SNP_pval "${SNP_PVAL}" \
            -doMajorMinor "${DO_MAJORMINOR}" \
            -minMapQ "${MIN_MAPQ}" \
            -minQ "${MIN_BASEQUAL}" \
            -nThreads "${N_CORES}"
    #   Are we missing a definiton for regions?
    elif [[ -z "${REGIONS}" ]]
    then
        "${ANGSD_DIR}"/angsd \
            -bam "${SAMPLE_LIST}" \
            -doGLF "${DO_GLF}" \
            -GL "${GT_LIKELIHOOD}" \
            -out "${OUT}"/"${PROJECT}" \
            -ref "${REF_SEQ}" \
            -anc "${ANC_SEQ}" \
            -doMaf "${DO_MAF}" \
            -SNP_pval "${SNP_PVAL}" \
            -doMajorMinor "${DO_MAJORMINOR}" \
            -minMapQ "${MIN_MAPQ}" \
            -minQ "${MIN_BASEQUAL}" \
            -nThreads "${N_CORES}"
    #   Assuming a single reigon was defined in config file
    else
        "${ANGSD_DIR}"/angsd \
            -bam "${SAMPLE_LIST}" \
            -r "${REGIONS}" \
            -doGLF "${DO_GLF}" \
            -GL "${GT_LIKELIHOOD}" \
            -out "${OUT}"/"${PROJECT}" \
            -ref "${REF_SEQ}" \
            -anc "${ANC_SEQ}" \
            -doMaf "${DO_MAF}" \
            -SNP_pval "${SNP_PVAL}" \
            -doMajorMinor "${DO_MAJORMINOR}" \
            -minMapQ "${MIN_MAPQ}" \
            -minQ "${MIN_BASEQUAL}" \
            -nThreads "${N_CORES}"
    fi
fi

N_SITES="`expr $(zcat "${OUT}"/${PROJECT}.mafs.gz | wc -l) - 1`"


zcat "${OUT}"/"${PROJECT}".glf.gz | "${NGSF_DIR}"/ngsF \
    -glf - \
    -out "${OUT}"/"${PROJECT}".approx_indF \
    -n_ind "${N_IND}" \
    -n_sites "${N_SITES}" \
    -min_epsilon "${MIN_EPSILON}" \
    -approx_EM \
    -seed "${SEED}" \
    -init_values r \
    -n_threads "${N_CORES}"

zcat "${OUT}"/"${PROJECT}".glf.gz | "${NGSF_DIR}"/ngsF \
    -glf - \
    -out "${OUT}"/"${PROJECT}".indF \
    -n_ind "${N_IND}" \
    -n_sites "${N_SITES}" \
    -min_epsilon "${MIN_EPSILON}" \
    -init_values "${OUT}"/"${PROJECT}".approx_indF.pars \
    -n_threads "${N_CORES}"
