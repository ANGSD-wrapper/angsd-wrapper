#!/usr/bin/env bash

set -e
set -o pipefail

#   Load variables from supplied config file
source "$1"

#   Are we using Common_Config? If so, source it
if [[ -f "${COMMON}" ]]
then
    source "${COMMON}"
fi

#   Where is angsd-wrapper located?
SOURCE="$2"

#   Where is ngsAdmix
ADMIX_DIR="${SOURCE}"/dependencies/ngsAdmix

#   Make the outdirectory
OUT="${SCRATCH}/${PROJECT}/Admixture"
mkdir -p "${OUT}"

#   Run ngsAdmix for all values of K from 2 to ${K}
for k in $(seq 2 "${K}")
do
    echo "Calculating for K=$k"
    WRAPPER_ARGS=$(echo -likes "${LIKELIHOOD}" \
        -K "${k}" \
        -P "${N_CORES}" \
        -o "${OUT}"/"${PROJECT}"."${k}" \
        -tol "${TOLERANCE}" \
        -minMaf "${MIN_MAF}")

    # Check for advanced arguments, and overwrite any overlapping definitions
    # Creates a single argument array from both inputs
    FINAL_ARGS=( "$(source "${SOURCE}/Wrappers/Arg_Zipper.sh" "${WRAPPER_ARGS}" "${ADVANCED_ARGS}")" )

    "${ADMIX_DIR}"/NGSadmix "${FINAL_ARGS[@]}"

done

#   Rename all ".qopt" files to ".qopt.graph.me"
while IFS= read -r -d '' qoptFile
do
    mv "${qoptFile}" "${qoptFile}.graph.me"
done < <(find "${OUT}" -name "*.qopt" -print0)

# for kvalue in $(eval echo {"$k"..2})
# do
#     cat "${OUT}"/"${PROJECT}"."$kvalue".qopt >> "${OUT}"/"${PROJECT}".allK.qopt
# done
