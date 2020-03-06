#!/usr/bin/env bash

set -e
set -o pipefail

#   Load variables from supplied config file
source "$1"

#   Where is angsd-wrapper located?
SOURCE=$2

#   Where is ANGSD located?
ANGSD_DIR=${SOURCE}/dependencies/angsd

#   Create the fasta file
#       Check to see if we're using -doCounts
if [[ -z "${DO_COUNTS}" ]]
then
    WRAPPER_ARGS=$(echo -doFasta "${DO_FASTA}" \
        -i "${ANC_BAM}"\
        -out "${OUT}"
else
    WRAPPER_ARGS=$(echo -doFasta "${DO_FASTA}" \
        -doCounts "${DO_COUNTS}" \
        -i "${ANC_BAM}"\
        -out "${OUT}"
fi
# Check for advanced arguments, and overwrite any overlapping definitions
FINAL_ARGS=$(source ${SOURCE}/Wrappers/Arg_Zipper.sh "${WRAPPER_ARGS}" "${ADVANCED_ARGS}")
# echo "Final arguments: ${FINAL_ARGS}" 1<&2
"${ANGSD_DIR}"/angsd ${FINAL_ARGS}

#   If we have SAMTools, might as well index
if `command -v samtools > /dev/null 2> /dev/null`
then
    echo "Indexing fasta file..."
    #   Do we need to unzip the fasta file?
    find "${OUT_DIR}" -name "${OUT_NAME}.fa*" | grep .gz > /dev/null 2> /dev/null
    if [[ "$?" -eq 0 ]]
    then
        gzip -d "${OUT}".fa.gz
    fi
    samtools faidx "${OUT}".fa
else
    echo "Can't find SAMTools, we aren't indexing for you"
    echo "If you'd like to index yourself,"
    echo "Then install SAMTools and run the following command:"
    echo "samtools faidx ${OUT}/${OUT_NAME}.fa"
fi
