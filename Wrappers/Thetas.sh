#!/usr/bin/env bash

set -e
set -u
set -o pipefail

#   Load variables from supplied config file
source $1

#   Variables created from transforming other variables
#       The number of individuals in the taxon we are analyzing
N_IND = `wc -l < "${SAMPLE_LIST}"`
#       How many inbreeding coefficients are supplied?
N_F = `wc -l < "${SAMPLE_INBREEDING}"`
#       For ANGSD, the actual sample size is twice the number of individuals, since each individual has two chromosomes.
#       The individual inbreeding coefficents take care of the mismatch between these two numbers
N_CHROM = `expr 2 \* "${N_IND}"`

#   Perform a check to see if number of individuals matches number of inbreeding coefficients
if [ "${N_IND}" -ne "${N_F}" ]
then
    echo "Mismatch between number of samples in ${SAMPLE_LIST} and ${SAMPLE_INBREEDING}"
    exit 1
fi

#   Create outdirectory
mkdir -p ${OUTDIR}

if [[  -f "${OUTDIR}/${PROJECT}_Diversity.mafs.gz" ]] && [ "${OVERRIDE}" = "false" ]; then
    echo "maf already exists and OVERRIDE=false, skipping angsd -bam...";
else
    #   Now we actually run the command, this creates a binary file that contains the prior SFS
    if [[ "${REGIONS}" == */* ]]
    then
        "${ANGSD_DIR}"/angsd \
            -bam "${SAMPLE_LIST}"\
            -out "${OUTDIR}"/"${PROJECT}"_Diversity\
            -indF "${SAMPLE_INBREEDING}"\
            -doSaf "${DO_SAF}"\
            -doThetas "${DO_THETAS}"\
            -uniqueOnly "${UNIQUE_ONLY}"\
            -anc "${ANC_SEQ}"\
            -minMapQ "${MIN_MAPQ}"\
            -minQ "${MIN_BASEQUAL}"\
            -nInd "${N_IND}"\
            -minInd "${MIN_IND}"\
            -baq "${BAQ}"\
            -ref "${REF_SEQ}"\
            -GL "${GT_LIKELIHOOD}"\
            -P "${N_CORES}"\
            -doMajorMinor "${DO_MAJORMINOR}"\
            -doMaf "${DO_MAF}"\
            -pest "${PEST}"\
            -rf "${REGIONS}" \
            -doPost "${DO POST}"
    elif [[ -z "${REGIONS}" ]]
    then
        "${ANGSD_DIR}"/angsd \
            -bam "${SAMPLE_LIST}"\
            -out "${OUTDIR}"/"${PROJECT}"_Diversity\
            -indF "${SAMPLE_INBREEDING}"\
            -doSaf "${DO_SAF}"\
            -doThetas "${DO_THETAS}"\
            -uniqueOnly "${UNIQUE_ONLY}"\
            -anc "${ANC_SEQ}"\
            -minMapQ "${MIN_MAPQ}"\
            -minQ "${MIN_BASEQUAL}"\
            -nInd "${N_IND}"\
            -minInd "${MIN_IND}"\
            -baq "${BAQ}"\
            -ref "${REF_SEQ}"\
            -GL "${GT_LIKELIHOOD}"\
            -P "${N_CORES}"\
            -doMajorMinor "${DO_MAJORMINOR}"\
            -doMaf "${DO_MAF}"\
            -pest "${PEST}"\
            -doPost "${DO POST}"
    else
        "${ANGSD_DIR}"/angsd \
        -bam "${SAMPLE_LIST}"\
        -out "${OUTDIR}"/"${PROJECT}"_Diversity\
        -indF "${SAMPLE_INBREEDING}"\
        -doSaf "${DO_SAF}"\
        -doThetas "${DO_THETAS}"\
        -uniqueOnly "${UNIQUE_ONLY}"\
        -anc "${ANC_SEQ}"\
        -minMapQ "${MIN_MAPQ}"\
        -minQ "${MIN_BASEQUAL}"\
        -nInd "${N_IND}"\
        -minInd "${MIN_IND}"\
        -baq "${BAQ}"\
        -ref "${REF_SEQ}"\
        -GL "${GT_LIKELIHOOD}"\
        -P "${N_CORES}"\
        -doMajorMinor "${DO_MAJORMINOR}"\
        -doMaf "${DO_MAF}"\
        -pest "${PEST}" \
        -r "${REGIONS}"
    fi
fi


"${ANGSD_DIR}"/misc/thetaStat make_bed \
    "${OUTDIR}"/"${PROJECT}"_Diversity.thetas.gz

if [ "${SLIDING_WINDOW}" = "false" ]
then
    "${ANGSD_DIR}"/misc/thetaStat do_stat \
        "${OUTDIR}"/"${PROJECT}"_Diversity.thetas.gz \
        -nChr "${N_CHROM}"
else
    "${ANGSD_DIR}"/misc/thetaStat do_stat \
        "${OUTDIR}"/"${PROJECT}"_Diversity.thetas.gz \
        -nChr "${N_CHROM}" \
        -win "${WIN}" \
        -step "${STEP}"
fi
