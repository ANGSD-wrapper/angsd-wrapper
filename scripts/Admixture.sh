#!/usr/bin/env bash

set -e
set -u

#   source the common config file
source scripts/Common_Variables.conf
# load utils functions
source ${SCRIPTS_DIR}/utils.sh

K=5
N_CORES=32
MIN_MAF=0.05

load_config $1

LIKELIHOOD=${RESULTS_DIR}/${TAXON}_SFSOut.beagle.gz

for ((k=2; k<=K; k++))
do
    echo "Calculating for K="$k
	${NGS_POPGEN_DIR}/NGSadmix\
		-likes ${LIKELIHOOD}\
		-K ${k}\
		-P ${N_CORES}\
		-o ${TAXON}.${k}\
		-minMaf ${MIN_MAF}
done
