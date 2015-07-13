#!/usr/bin/env bash

set -e
set -u
set -o pipefail

K=5
N_CORES=32
MIN_MAF=0.05

source $1

LIKELIHOOD=${RESULTS_DIR}/${TAXON}_SFSOut.beagle.gz

for ((k=2; k<=K; k++))
do
    echo "Calculating for K="$k
	${NGS_DIR}/NGSadmix\
		-likes ${LIKELIHOOD}\
		-K ${k}\
		-P ${N_CORES}\
		-o ${TAXON}.${k}\
		-minMaf ${MIN_MAF}
done
