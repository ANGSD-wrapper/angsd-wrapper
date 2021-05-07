#! /usr/bin/env bash

set -e
set -o pipefail

#   Load variables from supplied config file
source "$1"

#   Are we using Common_Config? If so, source it
if [[ -f "${COMMON}" ]]
then
    source "${COMMON}"
fi

# If a regions file has been defined, determine the independent chromosomal regions
if [[ ! -z REGIONS ]]
then
    # Create an array of all unique region names
    # TODO: Doesn't account for spaces, strange characters
    regs=$(cut -f 1 -d ':' "$REGIONS" | sort -u)
fi

# Destination for copied config files
TMP_DIR="${TMP_DIR}_copied"
mkdir -p "${TMP_DIR}"

# Create and edit copies of each config file and its region
mkdir -p "${TMP_DIR}/$(basename $1)"
mkdir -p "${TMP_DIR}/Regions"
while IFS= read -r -d ' ' reg
do
    # Split regions file into multiple files in tmp_dir
    grep "${reg} "${REGIONS}" > "${TMP_DIR}/Regions/${reg}"

    # Create many duplicate config files
    cp "$1" "${TMP_DIR}/$(basename $1)/$(basename $1)_${reg}"

    # Add line to end of each file to overwrite existing region definition
    # Avoids multiple overwrites via sed or awk
    echo "REGIONS=${TMP_DIR}/Regions/${reg}" >> "${TMP_DIR}/${basename $1}/$(basename $1)_${reg}"
done < <(regs)


# If a common config file was found, then create a directory for
# copying common config files and separate them by region
if [[ -f "${COMMON}" ]]
then
    mkdir -p "${TMP_DIR}/Common"
    while IFS= read -r -d ' ' reg
    do
	cp "${COMMON}" "${TMP_DIR}/Common/$(basename "${COMMON}")_${reg}"
    done < <(regs)
fi
