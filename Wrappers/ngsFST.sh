#!/usr/bin/env bash

set -e
set -u
set -o pipefail

#   Load variables from supplied config file
source $1

#   Variables created from transforming other variables
#       The number of individuals in the groups we are analyzing
N_IND_1=`wc -l < "${G1_SAMPLE_LIST}"`
N_IND_2=`wc -l < "${G2_SAMPLE_LIST}"`

#   Check for SFS, intersect file, and 2DSFS from 2D_Site_Frequency_Spectrum
#   Exit with error if any of them don't exist
if [[ -f "${G1_SFS}" ]]
then
    >&2 echo "WRAPPER: saf for ${GROUP_1} exists, continuing to check for ${GROUP_2} saf..."
else
    >&2 echo "WRAPPER: saf for ${GROUP_1} does not exist, exiting..." >&2
    exit 1
fi

if [[ -f "${G2_SFS}" ]]
then
    >&2 echo "WRAPPER: saf for ${GROUP_2} exists, continuing to check for intersect..."
else
    >&2 echo "WRAPPER: saf for ${GROUP_2} does not exist, exiting..." >&2
    exit 1
fi

if [[ -f "${INTERSECT}" ]]
then
    >&2 echo "WRAPPER: intersect exists, continuing to check for 2dsfs..."
else
    >&2 echo "WRAPPER: intersect does not exist, exiting..." >&2
    exit 1
fi

if [[ -f "${2DSFS}" ]]
then
    >&2 echo "WRAPPER: 2dsfs exists, continuing to check for intersect..."
else
    >&2 echo "WRAPPER: 2dsfs does not exist, exiting..." >&2
    exit 1
fi

#   Get number of sites and individuals
N_SITES=`wc -l < "${INTERSECT}"`

#   Convert 2DSFS file for ngsPopGen use
Rscript ${ANGSD_WRAPPER}/Wrappers/SFS_to_FST.R ${2DSFS} \
    > ${SCRATCH}/${PROJECT}/2DSFS_Intergenic.${GROUP_1}.${GROUP_2}.converted.sfs

#   Get the FST
#       This is what needs the new ngsPopGen
#       Make sure to change this path when wrapper is updated
${NGS_POPGEN}/ngsPopGen/ngsFST \
    -postfiles ${G1_SFS} ${POP2_SFS} \
    -priorfile ${RESULTS_DIR}/2DSFS_Intergenic.${GROUP_1}.${GROUP_2}.converted.sfs \
    -nind ${N_IND_1} ${N_IND_2} \
    -nsites ${N_SITES} \
    -outfile ${RESULTS_DIR}/${GROUP_1}.${GROUP_2}.fst
