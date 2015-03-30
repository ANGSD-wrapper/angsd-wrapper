#!/usr/bin/env bash

set -e
set -u

#   source common config file
source scripts/common.conf
# load utils functions
source ${SCRIPTS_DIR}/utils.sh

#NEED TO DO THIS: create a priorfile from pops (2D sfs): ngsPopGen/ngs2dSFS -postfiles results/og_allopatric_SFSOut.saf results/og_SFSOut.saf -outfile results/og_2dSFS.txt -relative 1 -nind 4 4
TAXON1_LIST=${DATA_DIR}/${TAXON1}_samples.txt
TAXON2_LIST=${DATA_DIR}/${TAXON2}_samples.txt

NCORES=32
GT_LIKELIHOOD=1
DO_SAF=1

NSITES=100000
BLOCK_SIZE=20000
IS_LOG=1
RELATIVE=1


load_config $1


${ANGSD}/angsd 
    -b ${TAXON1_LIST} 
    -anc ${ANC_SEQ} 
    -out ${TAXON1} 
    -P {$N_CORES} 
    -r ${REGIONS} 
    -GL ${GT_LIKELIHOOD} 
    -doSaf ${DO_SAF} 
    -sites ${RESULTS_DIR}/${TAXON}_intersect.txt

${ANGSD}/angsd 
    -b ${TAXON2_LIST} 
    -anc ${ANC_SEQ} 
    -out ${TAXON1} 
    -P {$N_CORES} 
    -r ${REGIONS} 
    -GL ${GT_LIKELIHOOD} 
    -doSaf ${DO_SAF} 
    -sites ${RESULTS_DIR}/${TAXON}_intersect.txt

N_SITES=`wc -l ${RESULTS_DIR}/intersect.txt | cut -f 1 -d " "`
NIND1=`wc -l ${TAXON1_LIST} | cut -f 1 -d " "`
NIND2=`wc -l ${TAXON2_LIST} | cut -f 1 -d " "`

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
