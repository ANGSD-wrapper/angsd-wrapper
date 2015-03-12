#!/usr/bin/env bash

set -e
set -u

#   source the common config file
source scripts/common.conf
# load utils functions
source ${SCRIPTS_DIR}/utils.sh

K=5
N_CORES=32
MIN_MAF=0.05

load_config $1

for k in {2..$K}
	${ngsPopGen}/NGSadmix \
		-likes ${LIKELIHOOD}\
		-K ${k}\
		-P ${N_CORES}\
		-o ${TAXON}.$k\
		-minMaf ${MIN_MAF}
done