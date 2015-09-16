#!/usr/bin/env bash

set -e
set -u

#   source common config file
source scripts/common.conf
# load utils functions
source ${SCRIPTS_DIR}/utils.sh

#Need to have run 2DSFS or have files in correct name format

DO_SAF=2
UNIQUE_ONLY=0
MIN_BASEQUAL=20
BAQ=1
MIN_IND1=1
MIN_IND2=1
GT_LIKELIHOOD=2
MIN_MAPQ=30
N_CORES=16
DO_MAJORMINOR=1
DO_MAF=1
BLOCK_SIZE=20000

load_config $1

TAXON1_LIST=${DATA_DIR}/${TAXON1}_samples.txt
TAXON2_LIST=${DATA_DIR}/${TAXON2}_samples.txt
POP1_SFS=${RESULTS_DIR}/${TAXON1}_Intergenic_Conditioned.saf
POP2_SFS=${RESULTS_DIR}/${TAXON2}_Intergenic_Conditioned.saf
INTERSECT=${RESULTS_DIR}/intersect.${TAXON1}.${TAXON2}_intergenic.txt
TWODSFS=${RESULTS_DIR}/2DSFS_Intergenic.${TAXON1}.${TAXON2}.sfs

N_IND_1=`wc -l < ${TAXON1_LIST}`
N_IND_2=`wc -l < ${TAXON2_LIST}`

#check for sfs, intersect file, and 2dsfs from 2DSFS.sh
#exit with error if any don't exist
if file_exists "${POP1_SFS}"; then
    >&2 echo "WRAPPER: saf for ${TAXON1} exists, continuing to check for ${TAXON2} saf..."
else >&2 echo "WRAPPER: saf for ${TAXON1} does not exist, exiting..." >&2; exit 1
fi

if file_exists "${POP2_SFS}"; then
    >&2 echo "WRAPPER: saf for ${TAXON2} exists, continuing to check for intersect..."
else >&2 echo "WRAPPER: saf for ${TAXON2} does not exist, exiting..." >&2; exit 1
fi

if file_exists "${INTERSECT}"; then
    >&2 echo "WRAPPER: intersect exists, continuing to check for 2dsfs..."
else >&2 echo "WRAPPER: intersect does not exist, exiting..." >&2; exit 1
fi

if file_exists "${TWODSFS}"; then
    >&2 echo "WRAPPER: 2dsfs exists, continuing to check for intersect..."
else >&2 echo "WRAPPER: 2dsfs does not exist, exiting..." >&2; exit 1
fi

# get number of sites and individuals
N_SITES=`wc -l ${INTERSECT} | cut -f 1 -d " "`

# convert ANGSD 2DSFS for ngsPopGen use
Rscript ${SCRIPTS_DIR}/SFS_to_FST.R ${TWODSFS}\
    > ${RESULTS_DIR}/2DSFS_Intergenic.${TAXON1}.${TAXON2}.converted.sfs

# get FST
# this is what needs the new ngsPopGen
# make sure to change this path when wrapper is updated
${NGS_TOOLS_DIR}/ngsPopGen/ngsFST\
    -postfiles ${POP1_SFS} ${POP2_SFS}\
    -priorfile ${RESULTS_DIR}/2DSFS_Intergenic.${TAXON1}.${TAXON2}.converted.sfs\
    -nind ${N_IND_1} ${N_IND_2}\
    -nsites ${N_SITES}\
    -outfile ${RESULTS_DIR}/${TAXON1}.${TAXON2}.fst
