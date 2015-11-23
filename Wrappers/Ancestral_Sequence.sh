#!/usr/bin/env bash

set -e
set -o pipefail

#   Load variables from supplied config file
source $1

#   Where is angsd-wrapper located?
SOURCE=$2

#   Where is ANGSD located?
ANGSD_DIR=${SOURCE}/dependencies/angsd

#   Create the fasta file
#       Check to see if we're using -doCounts
if [[ -z "${DO_COUNTS}" ]]
then
    "${ANGSD_DIR}"/angsd \
        -doFasta "${DO_FASTA}" \
        -i "${ANC_BAM}"\
        -out "${OUT}"
else
    "${ANGSD_DIR}"/angsd \
        -doFasta "${DO_FASTA}" \
        -doCounts "${DO_COUNTS}" \
        -i "${ANC_BAM}"\
        -out "${OUT}"
fi

#   If we have SAMTools, might as well index
if `command -v samtools > /dev/null 2> /dev/null`
then
    echo "Indexing fasta file..."
    #   Do we need to unzip the fasta file?
    find "${OUT}" -name "${OUTNAME}.fa*" | grep .gz > /dev/null 2> /dev/null
    if [[ "$?" -eq 0 ]]
    then
        gzip -d "${OUT}".fa.gz
    fi
    samtools faidx "${OUT}".fa
else
    echo "Can't find SAMTools, we aren't indexing for you"
    echo "If you'd like to index yourself,"
    echo "Then install SAMTools and run the following command:"
    echo "samtools faidx ${OUT}/${OUTNAME}.fa"
