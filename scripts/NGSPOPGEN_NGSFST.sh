#!/usr/bin/env bash

set -e
set -u

#   source common config file
source common.conf
# load utils functions
source ${SCRIPTS_DIR}/utils.sh

NSITES=100000
NIND1=20
NIND2=20
BLOCK_SIZE=20000
IS_LOG=1

load_config $1

${NGS_POPGEN_DIR}/ngsFST\
    -postfiles ${POP1_SFS} ${POP2_SFS}\
    -priorfiles ${SPECTRUM_1} ${SPECTRUM_2}\
    -nind ${NIND1} ${NIND2}\
    -nsites ${NSITES}\
    -block_size ${BLOCK_SIZE}\
    -outfile ${RESULTS_DIR}/pops.fst\ 
    -islog ${IS_LOG}
