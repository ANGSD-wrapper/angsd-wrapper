#!/usr/bin/env bash
set -e
set -u

#   Source the common configuration file
source scripts/common.conf
# load utils functions
source ${SCRIPTS_DIR}/utils.sh

load_config $1

N_IND=`wc -l < ${TAXON_LIST}`


${ANGSD_DIR}/angsd\
  -bam ${TAXON_LIST}\
  -out ./results/${TAXON}_snps\
  -doMajorMinor 1\
  -uniqueOnly 0\
  -minMapQ 30\
  -minQ 20\
  -GL 2\
  -r 12:\
  -doGeno 7\
  -doPost 1\
  -postCutoff 0.95\
  -doMaf 2\
  -SNP_pval 1e-6\
  -nInd ${N_IND}\
  -minInd 4\
  -P 16
#  -ref ${REF_SEQ}\
#  -anc ${ANC_SEQ}
#  -baq 1
