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

 "${SOURCE}"/dependencies/angsd/angsd \
    -i "${ANC_BAM}"\
    -doFasta 1
