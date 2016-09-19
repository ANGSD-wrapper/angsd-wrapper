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
    "${ADMIX_DIR}"/NGSadmix \
        -likes "${LIKELIHOOD}" \
        -K "${k}" \
        -P "${N_CORES}" \
        -o "${OUT}"/"${PROJECT}"."${k}" \
        -tol "${TOLERANCE}" \
        -minMaf "${MIN_MAF}"
done

#   Rename all ".qopt" files to ".qopt.graph.me"
for qoptFile in $(find ${OUT} -name "*.qopt")
do
    mv "${qoptFile}" "${qoptFile}.graph.me"
done

# for kvalue in $(eval echo {"$k"..2})
# do
#     cat "${OUT}"/"${PROJECT}"."$kvalue".qopt >> "${OUT}"/"${PROJECT}".allK.qopt
# done
