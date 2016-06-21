#!/bin/bash

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
ANGSD_DIR="${SOURCE}"/dependencies/angsd

#   Where is ngsFST?
NGS_POPGEN="${SOURCE}"/dependencies/ngsPopGen

#   Create the out directory
OUT="${SCRATCH}"/"${PROJECT}"/Fst
mkdir -p "${OUT}"

#   Variables created from transforming other variables
#       The number of individuals in the groups we are analyzing
N_IND1=$(wc -l < "${G1_SAMPLE_LIST}")
N_IND2=$(wc -l < "${G2_SAMPLE_LIST}")
#       How many inbreeding coefficients are supplied?
N_F1=$(wc -l < "${G1_INBREEDING}")
N_F2=$(wc -l < "${G2_INBREEDING}")

#   Perform a check to see if number of individuals matches number of inbreeding coefficients
if [ "${N_IND1}" -ne "${N_F1}" ] || [ "${N_IND2}" -ne "${N_F2}" ]
then
    echo "Mismatch between number of samples and inbreeding coefficients"
    exit 1
fi

#   Sorting not working right now...
#   Sort our regions file, if provided
if [[ -f ${REGIONS} ]]
then
    if ! $(command -v Rscript > /dev/null 2> /dev/null); then echo "Please install R and Rscript and place in your PATH" >&2; exit 1; fi
    echo "Sorting ${REGIONS} to match the order in ${REF_SEQ}" >&2
    REF_EXT=$(echo "${REF_SEQ}" | rev | cut -f 1 -d '.' | rev)
    FAI=$(find $(dirname "${REF_SEQ}") -name "$(basename "$REF_SEQ" ".${REF_EXT}")*.fai")
    echo $FAI
    exit 8
    Rscript "${SOURCE}"/Wrappers/sortRegions.R "${REGIONS}" "${FAI}"
    REGIONS="$(find $(dirname ${REGIONS}) -name "*_sorted.txt")"
fi

echo "${REGIONS}"
exit 8

#   Now we actually run the command, this creates a binary file that contains the prior SFS
#       For 1st group
if [[ -f "${OUT}"/"${GROUP_1}_Intergenic.saf" ]] && [ "$OVERRIDE" = "false" ]
then
    echo "WRAPPER: ${GROUP_1} saf already exists and OVERRIDE=false, skipping angsd -bam..." >&2
else
    if [[ -f "${REGIONS}" ]]
    then
        echo "WRAPPER: $GROUP_1 sfs starting..." >&2
        "${ANGSD_DIR}"/angsd \
            -bam "${G1_SAMPLE_LIST}" \
            -out "${OUT}"/"${GROUP_1}"_Intergenic \
            -doMajorMinor "${DO_MAJORMINOR}" \
            -doMaf "${DO_MAF}" \
            -indF "${G1_INBREEDING}" \
            -doSaf "${DO_SAF}" \
            -uniqueOnly "${UNIQUE_ONLY}" \
            -anc "${ANC_SEQ}" \
            -minMapQ "${MIN_MAPQ}" \
            -minQ "${MIN_BASEQUAL}" \
            -nInd "${N_IND1}" \
            -minInd "${MIN_IND1}" \
            -baq "${BAQ}" \
            -ref "${REF_SEQ}" \
            -GL "${GT_LIKELIHOOD}" \
            -P "${N_CORES}" \
            -rf "${REGIONS}" \
            -doPost "${DO_POST}"
    elif [[ -z "${REGIONS}" ]]
    then
        echo "WRAPPER: $GROUP_1 sfs starting" >&2
        "${ANGSD_DIR}"/angsd \
            -bam "${G1_SAMPLE_LIST}" \
            -out "${OUT}"/"${GROUP_1}"_Intergenic \
            -doMajorMinor "${DO_MAJORMINOR}" \
            -doMaf "${DO_MAF}" \
            -indF "${G1_INBREEDING}" \
            -doSaf "${DO_SAF}" \
            -uniqueOnly "${UNIQUE_ONLY}" \
            -anc "${ANC_SEQ}" \
            -minMapQ "${MIN_MAPQ}" \
            -minQ "${MIN_BASEQUAL}" \
            -nInd "${N_IND1}" \
            -minInd "${MIN_IND1}" \
            -baq "${BAQ}" \
            -ref "${REF_SEQ}" \
            -GL "${GT_LIKELIHOOD}" \
            -P "${N_CORES}" \
            -doPost "${DO_POST}"
    else
        echo "WRAPPER: $GROUP_1 sfs starting" >&2
        "${ANGSD_DIR}"/angsd \
            -bam "${G1_SAMPLE_LIST}" \
            -out "${OUT}"/"${GROUP_1}"_Intergenic \
            -doMajorMinor "${DO_MAJORMINOR}" \
            -doMaf "${DO_MAF}" \
            -indF "${G1_INBREEDING}" \
            -doSaf "${DO_SAF}" \
            -uniqueOnly "${UNIQUE_ONLY}" \
            -anc "${ANC_SEQ}" \
            -minMapQ "${MIN_MAPQ}" \
            -minQ "${MIN_BASEQUAL}" \
            -nInd "${N_IND1}" \
            -minInd "${MIN_IND1}" \
            -baq "${BAQ}" \
            -ref "${REF_SEQ}" \
            -GL "${GT_LIKELIHOOD}" \
            -P "${N_CORES}" \
            -r "${REGIONS}"
    fi
fi

#   For 2nd group:
if [[ -f "${OUT}"/"${GROUP_2}_Intergenic.saf" ]] && [ "$OVERRIDE" = "false" ]
then
    echo "WRAPPER: ${GROUP_2} saf already exists and OVERRIDE=false, skipping angsd -bam..." >&2
else
    #   Do we have a regions file?
    if [[ -f "${REGIONS}" ]]
    then
        echo "WRAPPER: $GROUP_2 sfs starting..." >&2
        "${ANGSD_DIR}"/angsd \
            -bam "${G2_SAMPLE_LIST}" \
            -out "${OUT}"/"${GROUP_2}"_Intergenic \
            -doMajorMinor "${DO_MAJORMINOR}" \
            -doMaf "${DO_MAF}" \
            -indF "${G2_INBREEDING}" \
            -doSaf "${DO_SAF}" \
            -uniqueOnly "${UNIQUE_ONLY}" \
            -anc "${ANC_SEQ}" \
            -minMapQ "${MIN_MAPQ}" \
            -minQ "${MIN_BASEQUAL}" \
            -nInd "${N_IND2}" \
            -minInd "${MIN_IND2}" \
            -baq "${BAQ}" \
            -ref "${REF_SEQ}" \
            -GL "${GT_LIKELIHOOD}" \
            -P "${N_CORES}" \
            -rf "${REGIONS}" \
            -doPost "${DO_POST}"
    #   Are we missing a definiton for regions?
    elif [[ -z "${REGIONS}" ]]
    then
        echo "WRAPPER: $GROUP_2 sfs starting..." >&2
        "${ANGSD_DIR}"/angsd \
            -bam "${G2_SAMPLE_LIST}" \
            -out "${OUT}"/"${GROUP_2}"_Intergenic \
            -doMajorMinor "${DO_MAJORMINOR}" \
            -doMaf "${DO_MAF}" \
            -indF "${G2_INBREEDING}" \
            -doSaf "${DO_SAF}" \
            -uniqueOnly "${UNIQUE_ONLY}" \
            -anc "${ANC_SEQ}" \
            -minMapQ "${MIN_MAPQ}" \
            -minQ "${MIN_BASEQUAL}" \
            -nInd "${N_IND2}" \
            -minInd "${MIN_IND2}" \
            -baq "${BAQ}" \
            -ref "${REF_SEQ}" \
            -GL "${GT_LIKELIHOOD}" \
            -P "${N_CORES}" \
            -doPost "${DO_POST}"
    #   Assuming a single reigon was defined in config file
    else
        echo "WRAPPER: $GROUP_2 sfs starting..." >&2
        "${ANGSD_DIR}"/angsd \
            -bam "${G2_SAMPLE_LIST}" \
            -out "${OUT}"/"${GROUP_2}"_Intergenic \
            -doMajorMinor "${DO_MAJORMINOR}" \
            -doMaf "${DO_MAF}" \
            -indF "${G2_INBREEDING}" \
            -doSaf "${DO_SAF}" \
            -uniqueOnly "${UNIQUE_ONLY}" \
            -anc "${ANC_SEQ}" \
            -minMapQ "${MIN_MAPQ}" \
            -minQ "${MIN_BASEQUAL}" \
            -nInd "${N_IND2}" \
            -minInd "${MIN_IND2}" \
            -baq "${BAQ}" \
            -ref "${REF_SEQ}" \
            -GL "${GT_LIKELIHOOD}" \
            -P "${N_CORES}" \
            -r "${REGIONS}"
    fi
fi

#   Estimate joint SFS using realSFS
echo "WRAPPER: realSFS 2dsfs..." >&2
"${ANGSD_DIR}"/misc/realSFS "${OUT}"/"${GROUP_1}"_Intergenic.saf.idx "${OUT}"/"${GROUP_2}"_Intergenic.saf.idx \
    -P "${N_CORES}" \
    > "${OUT}"/2DSFS_Intergenic."${GROUP_1}"."${GROUP_2}".sfs

# prepare fst
echo "WRAPPER: realSFS fst prep..." >&2
"${ANGSD_DIR}"/misc/realSFS \
    fst index "${OUT}"/"${GROUP_1}"_Intergenic.saf.idx "${OUT}"/"${GROUP_2}"_Intergenic.saf.idx \
    -sfs "${OUT}"/2DSFS_Intergenic."${GROUP_1}"."${GROUP_2}".sfs \
    -fstout "${OUT}"/"${GROUP_1}"."${GROUP_2}"

#check if user wants global Fst estimate
if [ "${GLOBAL}"=="true" ]; then
    echo "WRAPPER: estimating global fst..." >&2
    "${ANGSD_DIR}"/misc/realSFS \
        fst stats "${OUT}"/"${GROUP_1}"."${GROUP_2}".fst.idx \
        > "${OUT}"/"${GROUP_1}"."${GROUP_2}".fst.global
else
    echo "WRAPPER: global fst not requested..." >&2
fi

#check if user wants windowed analysis
if [[ -n "${WIN}" ]] && [[ -n "${STEP}" ]]; then
    echo "WRAPPER: estimating windowed fst..." >&2
    "${ANGSD_DIR}"/misc/realSFS \
        fst stats2 "${OUT}"/"${GROUP_1}"."${GROUP_2}".fst.idx \
        -win "${WIN}" \
        -step "${STEP}" \
        > "${OUT}"/"${GROUP_1}"."${GROUP_2}".fst.slidingwindow
else
    echo "WRAPPER: windowed fst not requested..." >&2
fi

#########################################################################
# Deprecated
#########################################################################
#   Estimate the Fst using ngsFST
#   First, convert our 2D SFS output to the old output for ngsPopGen
#echo "WRAPPER: converting 2D SFS for Fst Estimations..." >&2
#"${ANGSD_DIR}"/misc/realSFS print \
#    "${OUT}"/"${GROUP_1}"_Intergenic.saf.idx \
#    "${OUT}"/"${GROUP_2}"_Intergenic.saf.idx \
#    -P "${N_CORES}" \
#    -oldout 1

#   Move the shared.pos.gz file to our out directory
#mv ${OUT/2DSFS}/shared.pos.gz ${OUT}

#   Unzip shared.pos.gz and get the number of shared sites
#gzip -df ${OUT}/shared.pos.gz
#N_SITES=`wc -l < "${OUT}"/shared.pos`

#   Generate a prior spectrum using ngs2dSFS from ngsPopGen
#echo "WRAPPER: generating spectrum..." >&2
#${NGS_POPGEN}/ngs2dSFS \
#    -postfiles ${OUT}/${GROUP_1}_Intergenic.saf ${OUT}/${GROUP_2}_Intergenic.saf \
#    -outfile ${OUT}/${GROUP_1}.${GROUP_2}.spectrum.txt \
#    -nind ${N_IND1} ${N_IND2} \
#    -relative ${RELATIVE} \
#    -maxlike ${MAX_LIKE} \
#    -block_size ${BLOCK_SIZE} \
#    -nsites ${N_SITES}

#   Calculate Fst using ngsFST
#echo "WRAPPER: estimating Fst..." >&2
#${NGS_POPGEN}/ngsFST \
#    -postfiles ${OUT}/${GROUP_1}_Intergenic.saf ${OUT}/${GROUP_2}_Intergenic.saf \
#    -priorfile ${OUT}/${GROUP_1}.${GROUP_2}.spectrum.txt \
#    -nind ${N_IND1} ${N_IND2} \
#    -block_size ${BLOCK_SIZE} \
#    -nsites ${N_SITES} \
#    -outfile ${OUT}/${GROUP_1}.${GROUP_2}.fst

#   Unzip the mafs files
#gzip -df ${OUT}/${GROUP_1}_Intergenic.mafs.gz
#gzip -df ${OUT}/${GROUP_2}_Intergenic.mafs.gz

#   Merge shared.pos file with Fst output file
#echo "WRAPPER: creating files for Shiny graphing..." >&2
#Rscript ${SOURCE}/Wrappers/fst_bp.R \
#    ${SOURCE} \
#    ${OUT}/shared.pos \
#    ${OUT}/${GROUP_1}.${GROUP_2}.fst \
#    ${OUT}/${GROUP_1}_Intergenic.mafs \
#    ${OUT}/${GROUP_2}_Intergenic.mafs \
#    ${OUT}/"${PROJECT}".Fst.graph.me
# why do we not zip the mafs again? sure this takes time, but the unzipped files are very large
