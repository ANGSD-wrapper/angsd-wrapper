#!/usr/bin/env bash

set -e
set -u
set -o pipefail

#   Load variables from supplied config file
source $1

 ${ANGSD_DIR}/angsd \
    -i ${ANC_BAM}\
    -doFasta 1