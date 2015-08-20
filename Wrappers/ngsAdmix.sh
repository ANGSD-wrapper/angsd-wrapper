#!/usr/bin/env bash

set -e
set -u
set -o pipefail

#   Load variables from supplied config file
source $1

#   Are we using Common_Config? If so, source it
if [[ -f "${COMMON}" ]]
then
    source "${COMMON}"
fi

#   Where is angsd-wrapper located?
SOURCE=$2

#   Where is ngsAdmix
ADMIX_DIR="${SOURCE}"/dependencies/ngsAdmix

for ((k=2; k<=K; k++))
do
    echo "Calculating for K=$k"
    "${ADMIX_DIR}"/NGSadmix \
        -likes "${LIKELIHOOD}" \
        -K "${k}" \
        -P "${N_CORES}" \
        -o "${PROJECT}"."${k}" \
        -minMaf "${MIN_MAF}"
done
