#! /usr/bin/env bash

set -e
set -u

source scripts/Common_Variables.conf

source ${SCRIPTS_DIR}/utils.sh

DO_FASTA=2
DO_COUNTS=1
DO_ABBABABA=1
UNIQUE_ONLY=0
MIN_BASEQUAL=20
BAQ=1
MIN_IND=1
MIN_MAPQ=30
N_CORES=32
CHECK_BAM_HEADERS=0
BLOCKSIZE=1000


load_config $1

#check if regions exist
if [[ ${REGIONS} == */*]]; then
    if file_exists "${REGIONS}" && file_not_empty "${REGIONS}"; then
        >&2 echo "WRAPPER: regions file exists and not empty, starting analysis..."
elif variable_exists "${REGIONS}"; then
    >&2 echo "WRAPPER: regions variable set, starting analysis..."
else
    >&2 echo "WRAPPER: regions not set, file does not exist, or file is empty. Exiting..." >&2; exit 1
fi

# extract consensus sequence of Tripsacum to be treated as outgroup
#angsd -doFasta 2 -doCounts 1 -i Tripsacum.bam -out Tripsacum
if [[ ${DO_CONSENSUS} == 1 ]]; then
    ${ANGSD_DIR}/angsd \
	-doFasta ${DO_FASTA}\
	-doCounts ${DO_COUNTS}\
	-i ${OUTGROUP}\
	-out ${RESULTS_DIR}/${TAXON}
fi

#-checkBamHeaders 0 -rf scaffoldNamesWithData.txt
#Those two flags are important, as in Tripsacum some small scaffolds have no reads mapped on. In that case, the consensus sequence file of Tripsacum (Tripsacum.fa) will have no sequence for those scaffolds.
#While those scaffold exist in the three ingroup bam files, angsd will scream "chromosome length inconsistency error".
#-checkBamHeaders 0 just ask it to skip checking the bam headers. The analyzed chromosomes are given by "-rf scaffoldNamesWithData.txt", which is generated by "cut -f1 Tripsacum.fa.fai" (the first column of samtools faidxed Tripsacum.fa file)
#-bam P3_luxurians.txt simply gives the directory of the three ingroup bam files
if [[ ${REGIONS} == */* ]]; then
	${ANGSD_DIR}/angsd \
		-doAbbababa ${DO_ABBABABA}\
		-blockSize ${BLOCKSIZE}\
		-doCounts ${DO_COUNTS}\
		-anc ${ANC_SEQ}\
		-bam ${TAXON_LIST}\
		-uniqueOnly ${UNIQUE_ONLY}\
		-minMapQ ${MIN_MAPQ}\
		-minQ ${MIN_BASEQUAL}\
		-minInd ${MIN_IND}\
		-P ${N_CORES}\
		-checkBamHeaders ${CHECK_BAM_HEADERS}\
		-rf ${REGIONS}\
		-out ${RESULTS_DIR}/${TAXON}.D
else
		${ANGSD_DIR}/angsd \
		-doAbbababa ${DO_ABBABABA}\
		-blockSize ${BLOCKSIZE}\
		-doCounts ${DO_COUNTS}\
		-anc ${ANC_SEQ}\
		-bam ${TAXON_LIST}\
		-uniqueOnly ${UNIQUE_ONLY}\
		-minMapQ ${MIN_MAPQ}\
		-minQ ${MIN_BASEQUAL}\
		-minInd ${MIN_IND}\
		-P ${N_CORES}\
		-checkBamHeaders ${CHECK_BAM_HEADERS}\
		-r ${REGIONS}\
		-out ${RESULTS_DIR}/${TAXON}.D
fi

#jackKnife.R is provided with angsd.
Rscript ${ANGSD_DIR}/R/jackKnife.R \
	file=${RESULTS_DIR}/${TAXON}.D.abbababa\
	indNames=${TAXON_LIST}\
	outfile=${RESULTS_DIR}/${TAXON}.abbababa
