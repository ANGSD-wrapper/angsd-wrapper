#!/usr/bin/env bash

set -e
set -u

#   source common config file
source scripts/common.conf
# load utils functions
source ${SCRIPTS_DIR}/utils.sh

#NEED TO DO THIS: create a priorfile from pops (2D sfs): ngsPopGen/ngs2dSFS -postfiles results/og_allopatric_SFSOut.saf results/og_SFSOut.saf -outfile results/og_2dSFS.txt -relative 1 -nind 4 4


NSITES=100000
BLOCK_SIZE=20000
IS_LOG=1
RELATIVE=1

load_config $1

${NGS_POPGEN_DIR}/ngs2dSFS\
    -postfiles ${POP1_SFS} ${POP2_SFS}\
    -outfile ${RESULTS_DIR}/${TAXON}_2dSFS.txt\
    -relative ${RELATIVE}\
    -nind ${NIND1} ${NIND2}\
    -nsites ${NSITES}


${NGS_POPGEN_DIR}/ngsFST\
    -postfiles ${POP1_SFS} ${POP2_SFS}\
    -priorfile ${RESULTS_DIR}/${TAXON}_2dSFS.txt\
    -nind ${NIND1} ${NIND2}\
    -nsites ${NSITES}\
    -block_size ${BLOCK_SIZE}\
    -outfile ${RESULTS_DIR}/${TAXON}_pops.fst\
    -islog ${IS_LOG}
