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
SOURCE=$2

#   Where is ANGSD?
ANGSD_DIR=${SOURCE}/dependencies/angsd

#   Variables created from transforming other variables
#       The number of individuals in the taxon we are analyzing
echo "${SAMPLE_LIST}" 1<&2
N_IND=$(wc -l < "${SAMPLE_LIST}")

#       How many inbreeding coefficients are supplied?
N_F=$(wc -l < "${SAMPLE_INBREEDING}")
#       For ANGSD, the actual sample size is twice the number of individuals, since each individual has two chromosomes.
#       The individual inbreeding coefficents take care of the mismatch between these two numbers
N_CHROM=$(( 2 * "${N_IND}" ))

#   Perform a check to see if number of individuals matches number of inbreeding coefficients
if [ "${N_IND}" -ne "${N_F}" ]
then
    echo "WRAPPER: Mismatch between number of samples in ${SAMPLE_LIST} and ${SAMPLE_INBREEDING}"
    exit 1
fi

#   Create outdirectory
OUT=${SCRATCH}/${PROJECT}/Thetas
mkdir -p "${OUT}"

# If the user-defined filepaths point to existing files, skip SFS calculations
if [[ -f "${SFS}" ]] && [[ -f "${SAF}" ]] ;
then
    echo "WRAPPER: Using user-selected SFS files."

# If files exist and overriding isn't allowed, skip all calculations
elif [[ -f "${OUT}"/../SFS/"${PROJECT}"_DerivedSFS.graph.me ]] && [[ -f "${OUT}"/../SFS/"${PROJECT}"_SFSOut.saf.idx ]] && [ "${OVERRIDE}" = "false" ] ;
then
    echo "WRAPPER: Required files already exist and override set to false, please assign SAF and SFS filepath values in the theta config file $1 or change the OVERRIDE value."

    exit 1

# If no valid paths are supplied and no conflicts exist, then calculate SFS
else
    #   Now we actually run the command, this creates a binary file that contains the prior SFS
    echo "WRAPPER: Generating new SFS files to inform Theta estimates, with overwrite permissions."

    # Generate new SFS calculations
    bash "${SOURCE}"/Wrappers/Site_Frequency_Spectrum.sh "$1" "${SOURCE}"

    SFS="${OUT}/../SFS/${PROJECT}"_DerivedSFS.graph.me
    SAF="${OUT}/../SFS/${PROJECT}_SFSOut.saf.idx"
fi

"${ANGSD_DIR}"/misc/realSFS saf2theta \
    "${SAF}" \
    -sfs "${SFS}" \
    -outname "${OUT}/${PROJECT}"_Diversity

if [ "${SLIDING_WINDOW}" = "false" ]
then
    "${ANGSD_DIR}"/misc/thetaStat do_stat \
        "${OUT}"/"${PROJECT}"_Diversity.thetas.idx
else
    "${ANGSD_DIR}"/misc/thetaStat do_stat \
        "${OUT}"/"${PROJECT}"_Diversity.thetas.idx \
        -win "${WIN}" \
        -step "${STEP}"
#         -outnames "${PROJECT}_Diversity.thetasWindow.gz"
fi

# Filter pestPG file for invariant sites
echo "WRAPPER: Creating files for Shiny graphing..." >&2
Rscript "${SOURCE}"/Wrappers/thetas_filtering.R \
    "${SOURCE}" \
    "${OUT}"/"${PROJECT}"_Diversity.thetas.idx.pestPG \
    "${OUT}"/"${PROJECT}"_Thetas.graph.me
