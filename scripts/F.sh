#!/usr/bin/env bash

set -e
set -u

#   Source the common configuration file
source scripts/common.conf
# load utils functions
source ${SCRIPTS_DIR}/utils.sh

#Defaults
MIN_BASEQUAL=20
GT_LIKELIHOOD=1
DO_GLF=3
MIN_MAPQ=30
N_CORES=32
DO_MAJORMINOR=1
DO_MAF=1
SNP_PVAL=1e-6
OVERRIDE=false
SEED=12345
MIN_EPSILON=1e-9

# load variables from supplied config file
load_config $1
N_IND=`wc -l < ${TAXON_LIST}`

#check if regions exist
if [[ ${REGIONS} == */*]]; then
    if file_exists "${REGIONS}" && file_not_empty "${REGIONS}"; then
        >&2 echo "WRAPPER: regions file exists and not empty, starting analysis..."
elif variable_exists "${REGIONS}"; then
    >&2 echo "WRAPPER: regions variable set, starting analysis..."
else
    >&2 echo "WRAPPER: regions not set, file does not exist, or file is empty. Exiting..." >&2; exit 1
fi

if [[ ${REGIONS} == */* ]]; then
	${ANGSD_DIR}/angsd \
		-bam ${TAXON_LIST}\
		-rf ${REGIONS}\
		-doGLF ${DO_GLF}\
		-GL ${GT_LIKELIHOOD}\
		-out ${RESULTS_DIR}/${TAXON}\
		-ref ${REF_SEQ}\
		-anc ${ANC_SEQ}\
		-doMaf ${DO_MAF}\
		-SNP_pval ${SNP_PVAL}\
		-doMajorMinor ${DO_MAJORMINOR}\
		-minMapQ ${MIN_MAPQ}\
		-minQ ${MIN_BASEQUAL}\
		-nThreads ${N_CORES}
else
       	${ANGSD_DIR}/angsd \
		-bam ${TAXON_LIST}\
		-r ${REGIONS}\
		-doGLF ${DO_GLF}\
		-GL ${GT_LIKELIHOOD}\
		-out ${RESULTS_DIR}/${TAXON}\
		-ref ${REF_SEQ}\
		-anc ${ANC_SEQ}\
		-doMaf ${DO_MAF}\
		-SNP_pval ${SNP_PVAL}\
		-doMajorMinor ${DO_MAJORMINOR}\
		-minMapQ ${MIN_MAPQ}\
		-minQ ${MIN_BASEQUAL}\
		-nThreads ${N_CORES}
fi

N_SITES=$((`zcat ${RESULTS_DIR}/${TAXON}.mafs.gz | wc -l`-1))


zcat ${RESULTS_DIR}/${TAXON}.glf.gz | ./ngsF/ngsF \
	-n_ind ${N_IND}\
	-n_sites ${N_SITES}\
	-min_epsilon ${MIN_EPSILON}\
	-glf -\
	-out ${RESULTS_DIR}/${TAXON}.approx_indF\
	-approx_EM\
	-seed ${SEED}\
	-init_values r\
	-n_threads ${N_CORES}

zcat ${RESULTS_DIR}/${TAXON}.glf.gz | ./ngsF/ngsF \
	-n_ind ${N_IND}\
	-n_sites ${N_SITES}\
	-min_epsilon ${MIN_EPSILON}\
	-glf -\
	-out ${RESULTS_DIR}/${TAXON}.indF\
	-init_values ${RESULTS_DIR}/${TAXON}.approx_indF.pars\
	-n_threads ${N_CORES}

