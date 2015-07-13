#!/usr/bin/env bash

set -e
set -u
set -o pipefail

#   Load variables from supplied config file
source $1
N_IND=`wc -l < "${SAMPLE_LIST}"`

mkdir -p ${OUTDIR}/${PROJECT}

#   Do we have a regions file?
if [[ -f "${REGIONS}" ]]
then
    "${ANGSD_DIR}"/angsd \
        -bam "${SAMPLE_LIST}" \
        -rf "${REGIONS}" \
        -doGLF "${DO_GLF}" \
        -GL "${GT_LIKELIHOOD}" \
        -out "${OUTDIR}"/"${PROJECT}" \
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
        -out "${OUTDIR}"/"${PROJECT}" \
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
        -out "${OUTDIR}"/"${PROJECT}" \
        -ref "${REF_SEQ}" \
        -anc "${ANC_SEQ}" \
        -doMaf "${DO_MAF}" \
        -SNP_pval "${SNP_PVAL}" \
        -doMajorMinor "${DO_MAJORMINOR}" \
        -minMapQ "${MIN_MAPQ}" \
        -minQ "${MIN_BASEQUAL}" \
        -nThreads "${N_CORES}"
fi

N_SITES=$(`zcat "${OUTDIR}"/"${PROJECT}".mafs.gz | wc -l`-1)


zcat "${OUTDIR}"/"${PROJECT}".glf.gz | "${NGSF_DIR}"/ngsF \
    -n_ind "${N_IND}" \
    -n_sites "${N_SITES}" \
    -min_epsilon "${MIN_EPSILON}" \
    -glf - \
    -out "${OUTDIR}"/"${PROJECT}".approx_indF \
    -approx_EM \
    -seed "${SEED}" \
    -init_values r \
    -n_threads "${N_CORES}"

zcat "${OUTDIR}"/"${PROJECT}".glf.gz | "${NGSF_DIR}"/ngsF \
    -n_ind "${N_IND}" \
    -n_sites "${N_SITES}" \
    -min_epsilon "${MIN_EPSILON}" \
    -glf - \
    -out "${OUTDIR}"/"${PROJECT}".indF \
    -init_values "${OUTDIR}"/"${PROJECT}".approx_indF.pars \
    -n_threads "${N_CORES}"

