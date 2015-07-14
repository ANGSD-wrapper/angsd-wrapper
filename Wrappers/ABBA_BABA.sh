#! /usr/bin/env bash

set -e
set -u
set -o pipefail

source $1

#   Extract consensus sequence to be treated as outgroup
if [[ "${DO_CONSENSUS}" == 1 ]]
then
    "${ANGSD_DIR}"/angsd \
        -doFasta "${DO_FASTA}" \
        -doCounts "${DO_COUNTS}" \
        -i "${OUTGROUP}" \
        -out "${SCRATCH}"/"${PROJECT}"
fi

#   Check for local R installation
if `command -v Rscript > /dev/null 2> /dev/null`
then
    echo "R is insalled"
else
    echo "R is not installed or not in your PATH"
    exit 1
fi

#   Create outdirectory
mkdir -p "${SCRATCH}"/"${PROJECT}"

#   Now we actually run the command
#   Do we have a regions file?
if [[ -f "${REGIONS}" ]]
then
    "${ANGSD_DIR}"/angsd \
        -doAbbababa "${DO_ABBABABA}" \
        -blockSize "${BLOCKSIZE}" \
        -doCounts "${DO_COUNTS}" \
        -anc "${ANC_SEQ}" \
        -bam "${SAMPLE_LIST}" \
        -uniqueOnly "${UNIQUE_ONLY}" \
        -minMapQ "${MIN_MAPQ}" \
        -minQ "${MIN_BASEQUAL}" \
        -minInd "${MIN_IND}" \
        -P "${N_CORES}" \
        -checkBamHeaders "${CHECK_BAM_HEADERS}" \
        -rf "${REGIONS}" \
        -out "${SCRATCH}"/"${PROJECT}"/"${PROJECT}".D
#   Are we missing a definiton for regions?
elif [[ -z "${REGIONS}" ]]
then
    "${ANGSD_DIR}"/angsd \
        -doAbbababa "${DO_ABBABABA}" \
        -blockSize "${BLOCKSIZE}" \
        -doCounts "${DO_COUNTS}" \
        -anc "${ANC_SEQ}" \
        -bam "${SAMPLE_LIST}" \
        -uniqueOnly "${UNIQUE_ONLY}" \
        -minMapQ "${MIN_MAPQ}" \
        -minQ "${MIN_BASEQUAL}" \
        -minInd "${MIN_IND}" \
        -P "${N_CORES}" \
        -checkBamHeaders "${CHECK_BAM_HEADERS}" \
        -out "${SCRATCH}"/"${PROJECT}"/"${PROJECT}".D
#   Assuming a single reigon was defined in config file
else
    "${ANGSD_DIR}"/angsd \
        -doAbbababa "${DO_ABBABABA}" \
        -blockSize "${BLOCKSIZE}" \
        -doCounts "${DO_COUNTS}" \
        -anc "${ANC_SEQ}" \
        -bam "${SAMPLE_LIST}" \
        -uniqueOnly "${UNIQUE_ONLY}" \
        -minMapQ "${MIN_MAPQ}" \
        -minQ "${MIN_BASEQUAL}" \
        -minInd "${MIN_IND}" \
        -P "${N_CORES}" \
        -checkBamHeaders "${CHECK_BAM_HEADERS}" \
        -r "${REGIONS}" \
        -out "${SCRATCH}"/"${PROJECT}"/"${PROJECT}".D
fi

#   jackKnife.R is provided with angsd.
Rscript "${ANGSD_DIR}"/R/jackKnife.R \
    file="${SCRATCH}"/"${PROJECT}"/"${PROJECT}".D.abbababa \
    indNames="${SAMPLE_LIST}" \
    outfile="${SCRATCH}"/"${PROJECT}"/"${PROJECT}".abbababa
