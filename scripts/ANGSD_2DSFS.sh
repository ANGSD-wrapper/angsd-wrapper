#!/bin/sh

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
REGIONS="1:"
OVERRIDE=false

# load variables from supplied config file
load_config $1

# calculated variables
TAXON_LIST1=data/${TAXON1}_samples.txt
TAXON_INBREEDING1=data/${TAXON1}_F.txt
N_IND1=`wc -l < ${TAXON_LIST1}`

N_CHROM1=`expr 2 \* ${N_IND1}`

TAXON_LIST2=data/${TAXON2}_samples.txt
TAXON_INBREEDING2=data/${TAXON2}_F.txt
N_IND2=`wc -l < ${TAXON_LIST2}`

N_CHROM2=`expr 2 \* ${N_IND2}`

# For 1st taxon
if file_exists "./results/${TAXON1}_Intergenic.mafs.gz" && [ "$OVERRIDE" = "false" ]; then 
    echo "maf already exists and OVERRIDE=false, skipping angsd -bam...";
else
    ${ANGSD_DIR}/angsd\
        -bam {$TAXON_LIST1}\
        -out results/${TAXON1}_Intergenic\
        -doMajorMinor ${DO_MAJORMINOR}\
        -doMaf ${DO_MAF}\
        -indF ${TAXON_INBREEDING1}\
        -doSaf ${DO_SAF}\
        -uniqueOnly ${UNIQUE_ONLY}\
        -anc ${ANC_SEQ}\
        -minMapQ ${MIN_MAPQ}\
        -minQ ${MIN_BASEQUAL}\
        -nInd ${N_IND1}\
        -minInd ${MIN_IND}\
        -baq ${BAQ}\
        -ref ${REF_SEQ}\
        -GL ${GT_LIKELIHOOD}\
        -P ${N_CORES}\
        -rf ${REGIONS}
fi

${ANGSD_DIR}/misc/realSFS\
    results/${TAXON1}_Intergenic.saf\
    ${N_CHROM1}\
    -P ${N_CORES}\
    > results/${TAXON1}_Intergenic.sfs 

# For 2nd taxon:
if file_exists "./results/${TAXON2}_Intergenic.mafs.gz" && [ "$OVERRIDE" = "fal\
se" ]; then
    echo "maf already exists and OVERRIDE=false, skipping angsd -bam...";
else
    ${ANGSD_DIR}/angsd\
        -bam {$TAXON_LIST2}\
        -out results/${TAXON2}_Intergenic\
        -doMajorMinor ${DO_MAJORMINOR}\
        -doMaf ${DO_MAF}\
        -indF ${TAXON_INBREEDING2}\
        -doSaf ${DO_SAF}\
        -uniqueOnly ${UNIQUE_ONLY}\
        -anc ${ANC_SEQ}\
        -minMapQ ${MIN_MAPQ}\
        -minQ ${MIN_BASEQUAL}\
        -nInd ${N_IND2}\
        -minInd ${MIN_IND}\
        -baq ${BAQ}\
        -ref ${REF_SEQ}\
        -GL ${GT_LIKELIHOOD}\
        -P ${N_CORES}\
        -rf ${REGIONS}
fi

${ANGSD_DIR}/misc/realSFS\
    results/${TAXON2}_Intergenic.saf\
    ${N_CHROM2}\
    -P ${N_CORES}\
    > results/${TAXON2}_Intergenic.sfs

# extract compressed files
gunzip -k results/${TAXON1}_Intergenic.saf.pos.gz
gunzip -k results/${TAXON2}_Intergenic.saf.pos.gz

# find positions that occur in both taxa with uniq
cat results/${TAXON1}_Intergenic.saf.pos results/${TAXON2}_Intergenic.saf.pos | sort | uniq -d > results/intersect.${TAXON1}.${TAXON2}_intergenic.txt

# calculate allele frequencies only on sites in both populations
    ${ANGSD_DIR}/angsd\
        -bam {$TAXON_LIST1}\
        -out results/${TAXON1}_Intergenic_Conditioned\
        -doMajorMinor ${DO_MAJORMINOR}\
        -doMaf ${DO_MAF}\
        -indF ${TAXON_INBREEDING1}\
        -doSaf ${DO_SAF}\
        -uniqueOnly ${UNIQUE_ONLY}\
        -anc ${ANC_SEQ}\
        -minMapQ ${MIN_MAPQ}\
        -minQ ${MIN_BASEQUAL}\
        -nInd ${N_IND1}\
        -minInd ${MIN_IND}\
        -baq ${BAQ}\
        -ref ${REF_SEQ}\
        -GL ${GT_LIKELIHOOD}\
        -P ${N_CORES}\
        -rf ${REGIONS}\
        -sites results/intersect.${TAXON1}.${TAXON2}_intergenic.txt

    ${ANGSD_DIR}/angsd\
        -bam {$TAXON_LIST2}\
        -out results/${TAXON2}_Intergenic_Conditioned\
        -doMajorMinor ${DO_MAJORMINOR}\
        -doMaf ${DO_MAF}\
        -indF ${TAXON_INBREEDING2}\
        -doSaf ${DO_SAF}\
        -uniqueOnly ${UNIQUE_ONLY}\
        -anc ${ANC_SEQ}\
        -minMapQ ${MIN_MAPQ}\
        -minQ ${MIN_BASEQUAL}\
        -nInd ${N_IND2}\
        -minInd ${MIN_IND}\
        -baq ${BAQ}\
        -ref ${REF_SEQ}\
        -GL ${GT_LIKELIHOOD}\
        -P ${N_CORES}\
        -rf ${REGIONS}\
        -sites results/intersect.${TAXON1}.${TAXON2}_intergenic.txt

# estimate joint SFS using realSFS
${ANGSD_DIR}/misc/realSFS 2dsfs\
    results/${TAXON1}_Intergenic_Conditioned.saf results/${TAXON2}_Intergenic_Conditioned.saf\
    ${N_IND1}\
    ${N_IND2}\
    -P ${N_CORES}\
    > results/2DSFS_Intergenic.${TAXON1}.${TAXON2}.sfs
