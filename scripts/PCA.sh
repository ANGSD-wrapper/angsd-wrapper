#!/usr/bin/env bash

set -e
set -u

#   Source the common configuration file
source scripts/Common_Variables.conf
# load utils functions
source ${SCRIPTS_DIR}/utils.sh

DO_MAF=2
DO_MAJORMINOR=1
DO_GENO=32
DO_POST=1
N_IND=`wc -l < ${TAXON_LIST}`
N_CORES=8
CALL=0
N_SITES=100000

load_config $1

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
        -rf ${REGIONS}
fi

gunzip ${RESULTS_DIR}/${TAXON}_PCA.geno.gz

ngsPopGen/ngsCovar
 -probfile results/all.test.geno
 -outfile pop.covar
 -nind 21
 -nsites 100000
 -call 0

${NGS_POPGEN_DIR}/ngsCovar \
	-probfile ${RESULTS_DIR}/${TAXON}_PCA.geno\
	-outfile ${RESULTS_DIR}/${TAXON}_PCA.covar\
	-nind ${N_IND}\
	-nsites ${N_SITES}\
	-call ${CALL}
