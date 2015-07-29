#!/usr/bin/env bash

set -e
set -u
set -o pipefail

#   Load variables from supplied config file
source $1

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
