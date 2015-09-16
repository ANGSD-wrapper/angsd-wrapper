#!/usr/bin/env bash
set -e
set -u

#   Source the common configuration file
source scripts/common.conf
# load utils functions
source ${SCRIPTS_DIR}/utils.sh

DO_MAJORMINOR=1
UNIQUE_ONLY=0
MIN_MAPQ=30
MIN_BASEQUAL=20
GT_LIKELIHOOD=2
DO_GENO=7
DO_POST=1
POST_CUTOFF=0.95
DO_MAF=2
SNP_PVAL=1e-6
MIN_IND=1
N_CORES=32

load_config $1

#check if regions exist
if [[ ${REGIONS} == */*]]; then
    if file_exists "${REGIONS}" && file_not_empty "${REGIONS}"; then
        >&2 echo "WRAPPER: regions file exists and not empty, starting analysis..."
elif variable_exists "${REGIONS}"; then
    >&2 echo "WRAPPER: regions variable set, starting analysis..."
else
    >&2 echo "WRAPPER: regions not set, file does not exist, or file is empty. Exiting..." >&2; exit 1
fi

N_IND=`wc -l < ${TAXON_LIST}`

if [[ ${REGIONS} == */* ]]; then
  >&2 echo "WRAPPER: $TAXON GL starting..."
  ${ANGSD_DIR}/angsd\
    -bam ${TAXON_LIST}\
    -out ./results/${TAXON}_snps\
    -doMajorMinor ${DO_MAJORMINOR}\
    -uniqueOnly ${UNIQUE_ONLY}\
    -minMapQ ${MIN_MAPQ}\
    -minQ ${MIN_BASEQUAL}\
    -GL ${GT_LIKELIHOOD}\
    -rf ${REGIONS}\
    -doGeno ${DO_GENO}\
    -doPost ${DO_POST}\
    -postCutoff ${POST_CUTOFF}\
    -doMaf ${DO_MAF}\
    -SNP_pval ${SNP_PVAL}\
    -nInd ${N_IND}\
    -minInd ${MIN_IND}\
    -P ${N_CORES}
else
  >&2 echo "WRAPPER: $TAXON GL starting"
  ${ANGSD_DIR}/angsd\
    -bam ${TAXON_LIST}\
    -out ./results/${TAXON}_snps\
    -doMajorMinor ${DO_MAJORMINOR}\
    -uniqueOnly ${UNIQUE_ONLY}\
    -minMapQ ${MIN_MAPQ}\
    -minQ ${MIN_BASEQUAL}\
    -GL ${GT_LIKELIHOOD}\
    -r ${REGIONS}\
    -doGeno ${DO_GENO}\
    -doPost ${DO_POST}\
    -postCutoff ${POST_CUTOFF}\
    -doMaf ${DO_MAF}\
    -SNP_pval ${SNP_PVAL}\
    -nInd ${N_IND}\
    -minInd ${MIN_IND}\
    -P ${N_CORES}
