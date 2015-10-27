#!/usr/bin/env bash

set -e
set -u

#   Source the common configuration file
source scripts/common.conf
# load utils functions
source ${SCRIPTS_DIR}/utils.sh

DO_MAF=2
DO_MAJORMINOR=1
DO_GENO=32
DO_POST=1
N_CORES=8
CALL=0
GT_LIKELIHOOD=1
N_SITES=100000

load_config $1

N_IND=`wc -l < ${TAXON_LIST}`

#check if regions exist
if [[ ${REGIONS} == */* ]]; then
    if file_exists "${REGIONS}" && file_not_empty "${REGIONS}"; then
        >&2 echo "WRAPPER: regions file exists and not empty, starting analysis..."
    elif variable_exists "${REGIONS}"; then
	>&2 echo "WRAPPER: regions variable set, starting analysis..."
    else
	>&2 echo "WRAPPER: regions not set, file does not exist, or file is empty. Exiting..." >&2; exit 1
    fi
fi

if [[ ${REGIONS} == */* ]]; then
	${ANGSD_DIR}/angsd \
        -bam ${TAXON_LIST}\
        -GL ${GT_LIKELIHOOD}\
        -out ${RESULTS_DIR}/${TAXON}_PCA\
        -doMajorMinor ${DO_MAJORMINOR}\
        -doMaf ${DO_MAF}\
        -doGeno ${DO_GENO}\
        -doPost ${DO_POST}\
        -nInd ${N_IND}\
        -P ${N_CORES}\
        -rf ${REGIONS}
else
	${ANGSD_DIR}/angsd \
        -bam ${TAXON_LIST}\
        -GL ${GT_LIKELIHOOD}\
        -out ${RESULTS_DIR}/${TAXON}_PCA\
        -doMajorMinor ${DO_MAJORMINOR}\
        -doMaf ${DO_MAF}\
        -doGeno ${DO_GENO}\
        -doPost ${DO_POST}\
        -nInd ${N_IND}\
        -P ${N_CORES}\
        -r ${REGIONS}
fi

gunzip ${RESULTS_DIR}/${TAXON}_PCA.geno.gz

${NGS_POPGEN_DIR}/ngsCovar \
	-probfile ${RESULTS_DIR}/${TAXON}_PCA.geno\
	-outfile ${RESULTS_DIR}/${TAXON}_PCA.covar\
	-nind ${N_IND}\
	-nsites ${N_SITES}\
	-call ${CALL}
