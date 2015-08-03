#!/usr/bin/env bash

set -e
set -u

#   Source the common configuration file
source scripts/Common_Variables.conf
# load utils functions
source ${SCRIPTS_DIR}/utils.sh

 ${ANGSD_DIR}/angsd \
	-i ${ANC_BAM}\
	-doFasta 1
