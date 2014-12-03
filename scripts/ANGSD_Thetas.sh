#!/usr/bin/env bash

set -e
set -u

# load utils functions
source ./scripts/utils.sh

#Defaults
DO_SAF=2
UNIQUE_ONLY=0
MIN_BASEQUAL=20
BAQ=1
MIN_IND=1
GT_LIKELIHOOD=2
MIN_MAPQ=30
N_CORES=32
DO_MAJORMINOR=1
DO_MAF=1
DO_THETAS=1
REGIONS="1:"
OVERRIDE=false
SLIDING_WINDOW=false
WIN=50000
STEP=10000

# load variables from supplied config file
load_config $1

# Variables created from transforming other variables
#   The number of individuals in the taxon we are analyzing
#   We use an embedded command to do this
#   ( wc -l < FILE will return just the line count of FILE,
#   rather than the line count and the filename. More efficient than piping
#   to a separate 'cut' process!)
TAXON_LIST=data/${TAXON}_samples.txt
TAXON_INBREEDING=data/${TAXON}_F.txt
N_IND=`wc -l < ${TAXON_LIST}`
#   For ANGSD, the actual sample size is twice the number of individuals, since
#   each individual has two chromosomes. The individual inbreeding coefficents
#   take care of the mismatch between these two numbers
N_CHROM=`expr 2 \* ${N_IND}`

# if directories don't exist, create them 
if directory_exists "./results"; then 
    echo "directories exist, skipping init.sh"; 
else 
    echo "creating directories..."; 
    bash ./scripts/init.sh; 
fi


if file_exists "./results/${TAXON}_Diversity.mafs.gz" && [ "$OVERRIDE" = "false" ]; then 
    echo "maf already exists and OVERRIDE=false, skipping angsd -bam...";
else
    #   Now we actually run the command, this creates a binary file that contains the prior SFS
    if [[ ${REGIONS} == */* ]]; then #if regions contains a / (is a path) then it's probably a region file. Therefore, use -rf
	${ANGSD_DIR}/angsd \
        -bam ${TAXON_LIST}\
        -out results/${TAXON}_Diversity\
        -indF ${TAXON_INBREEDING}\
        -doSaf ${DO_SAF}\
        -doThetas ${DO_THETAS}\
        -uniqueOnly ${UNIQUE_ONLY}\
        -anc ${ANC_SEQ}\
        -minMapQ ${MIN_MAPQ}\
        -minQ ${MIN_BASEQUAL}\
        -nInd ${N_IND}\
        -minInd ${MIN_IND}\
        -baq ${BAQ}\
        -ref ${REF_SEQ}\
        -GL ${GT_LIKELIHOOD}\
        -P ${N_CORES}\
        -doMajorMinor $DO_MAJORMINOR\
        -doMaf $DO_MAF\
        -pest ${PEST}\
        -rf ${REGIONS}
    else
        ${ANGSD_DIR}/angsd \
        -bam ${TAXON_LIST}\
        -out results/${TAXON}_Diversity\
        -indF ${TAXON_INBREEDING}\
        -doSaf ${DO_SAF}\
        -doThetas ${DO_THETAS}\
        -uniqueOnly ${UNIQUE_ONLY}\
        -anc ${ANC_SEQ}\
        -minMapQ ${MIN_MAPQ}\
        -minQ ${MIN_BASEQUAL}\
        -nInd ${N_IND}\
        -minInd ${MIN_IND}\
        -baq ${BAQ}\
        -ref ${REF_SEQ}\
        -GL ${GT_LIKELIHOOD}\
        -P ${N_CORES}\
        -doMajorMinor $DO_MAJORMINOR\
        -doMaf $DO_MAF\
        -pest ${PEST}\
        -r ${REGIONS}
    fi
fi


${ANGSD_DIR}/misc/thetaStat make_bed\
    results/${TAXON}_Diversity.thetas.gz
#    results/${TAXON}_Tajimas

if [ "$SLIDING_WINDOW" = "false" ]; then
    ${ANGSD_DIR}/misc/thetaStat do_stat\
        results/${TAXON}_Diversity.thetas.gz\
        -nChr ${N_CHROM}
else
    ${ANGSD_DIR}/misc/thetaStat do_stat\
        results/${TAXON}_Diversity.thetas.gz\
        -nChr ${N_CHROM}\
        -win ${WIN}\
        -step ${STEP}
fi
