#!/bin/sh

set -e
set -u

# load utils functions
source ./scripts/utils.sh

#Defaults
DO_SAF=2
UNIQUE_ONLY=1
MIN_BASEQUAL=20
BAQ=1
MIN_IND=1
GT_LIKELIHOOD=2
MIN_MAPQ=30
N_CORES=32
DO_MAJORMINOR=1
DO_MAF=1
REGIONS="1:"
OVERRIDE=false

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

#   Now we actually run the command, this creates a binary file that contains the prior SFS
if file_exists "./results/${TAXON}_SFSOut.mafs.gz" && [ "$OVERRIDE" = "false" ]; then 
    echo "maf already exists and OVERRIDE=false, skipping angsd -bam...";
else
    if [[ ${REGIONS} == */* ]]; then
	${ANGSD_DIR}/angsd \
        -bam ${TAXON_LIST}\
        -out results/${TAXON}_SFSOut\
        -indF ${TAXON_INBREEDING}\
        -doSaf ${DO_SAF}\
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
        -rf ${REGIONS}
    else
        ${ANGSD_DIR}/angsd \
        -bam ${TAXON_LIST}\
        -out results/${TAXON}_SFSOut\
        -indF ${TAXON_INBREEDING}\
        -doSaf ${DO_SAF}\
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
        -r ${REGIONS}
	fi
fi

${ANGSD_DIR}/misc/realSFS\
    results/${TAXON}_SFSOut.saf\
    ${N_CHROM}\
    -P ${N_CORES}\
    > results/${TAXON}_DerivedSFS
