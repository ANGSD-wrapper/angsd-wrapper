#!/bin/bash

set -e
set -o pipefail

#   Load variables from supplied config file
source $1

#   Are we using Common_Config? If so, source it
if [[ -f "${COMMON}" ]]
then
    source "${COMMON}"
fi

#   Where is angsd-wrapper located?
SOURCE=$2

#   Variables created from transforming other variables
#       The number of individuals in the groups we are analyzing
N_IND1=`wc -l < "${G1_SAMPLE_LIST}"`
N_IND2=`wc -l < "${G2_SAMPLE_LIST}"`
#       How many inbreeding coefficients are supplied?
N_F1=`wc -l < "${G1_INBREEDING}"`
N_F2=`wc -l < "${G2_INBREEDING}"`
#       For ANGSD, the actual sample size is twice the number of individuals, since each individual has two chromosomes.
#       The individual inbreeding coefficents take care of the mismatch between these two numbers
N_CHROM1=`expr 2 \* "${N_IND1}"`
N_CHROM2=`expr 2 \* "${N_IND2}"`

#   Perform a check to see if number of individuals matches number of inbreeding coefficients
if [ "${N_IND1}" -ne "${N_F1}" ] || [ "${N_IND2}" -ne "${N_F2}" ]
then
    echo "Mismatch between number of samples and inbreeding coefficients"
    exit 1
fi

#   Where is ANGSD?
ANGSD_DIR="${SOURCE}"/dependencies/angsd

#   Where is ngsFST?
NGS_POPGEN="${SOURCE}"/dependencies/ngsPopGen

#   Create the out directory
OUT=${SCRATCH}/"${PROJECT}"/2DSFS
mkdir -p ${OUT}

#   Now we actually run the command, this creates a binary file that contains the prior SFS
#       For 1st group
if [[ -f "${OUT}"/"${GROUP_1}_Intergenic.saf" ]] && [ "$OVERRIDE" = "false" ]
then
    echo "WRAPPER: saf already exists and OVERRIDE=false, skipping angsd -bam..." >&2
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
    echo "WRAPPER: saf already exists and OVERRIDE=false, skipping angsd -bam..." >&2
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

#   Find intersecting regions
echo "WRAPPER: making intersect file..." >&2
gunzip -c "${OUT}"/"${GROUP_1}"_Intergenic.saf.pos "${OUT}"/"${GROUP_2}"_Intergenic.saf.pos | sort | uniq -d | sort -k1,1 > "${OUT}"/Intersect."${GROUP_1}"."${GROUP_2}"_intergenic.txt

#   Calculate allele frequencies only on sites in both populations
#   Do we have a regions file?
if [[ -f "${REGIONS}" ]]
then
    echo "WRAPPER: $GROUP_1 sfs round 2..." >&2
    "${ANGSD_DIR}"/angsd \
        -bam "${G1_SAMPLE_LIST}" \
        -out "${OUT}"/"${GROUP_1}"_Intergenic_Conditioned \
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
        -doPost "${DO_POST}" \
        -sites "${OUT}"/Intersect."${GROUP_1}"."${GROUP_2}"_intergenic.txt
#   Are we missing a definiton for regions?
elif [[ -z "${REGIONS}" ]]
then
    echo "WRAPPER: $GROUP_1 sfs round 2..." >&2
    "${ANGSD_DIR}"/angsd \
        -bam "${G1_SAMPLE_LIST}" \
        -out "${OUT}"/"${GROUP_1}"_Intergenic_Conditioned \
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
        -doPost "${DO_POST}" \
        -sites "${OUT}"/Intersect."${GROUP_1}"."${GROUP_2}"_intergenic.txt
#   Assuming a single reigon was defined in config file
else
    echo "WRAPPER: $GROUP_1 sfs round 2..." >&2
    "${ANGSD_DIR}"/angsd \
        -bam "${G1_SAMPLE_LIST}" \
        -out "${OUT}"/"${GROUP_1}"_Intergenic_Conditioned \
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
        -r "${REGIONS}" \
        -sites "${OUT}"/Intersect."${GROUP_1}"."${GROUP_2}"_intergenic.txt
fi

#   Do we have a regions file?
if [[ -f "${REGIONS}" ]]
then
    echo "WRAPPER: $GROUP_2 sfs round 2..." >&2
    "${ANGSD_DIR}"/angsd \
        -bam "${G2_SAMPLE_LIST}" \
        -out "${OUT}"/"${GROUP_2}"_Intergenic_Conditioned \
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
        -doPost "${DO_POST}" \
        -sites "${OUT}"/Intersect."${GROUP_1}"."${GROUP_2}"_intergenic.txt
#   Are we missing a definiton for regions?
elif [[ -z "${REGIONS}" ]]
then
    echo "WRAPPER: $GROUP_2 sfs round 2..." >&2
    "${ANGSD_DIR}"/angsd \
        -bam "${G2_SAMPLE_LIST}" \
        -out "${OUT}"/"${GROUP_2}"_Intergenic_Conditioned \
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
        -doPost "${DO_POST}" \
        -sites "${OUT}"/Intersect."${GROUP_1}"."${GROUP_2}"_intergenic.txt
#   Assuming a single reigon was defined in config file
else
    echo "WRAPPER: $GROUP_2 sfs round 2..." >&2
    "${ANGSD_DIR}"/angsd \
        -bam "${G2_SAMPLE_LIST}" \
        -out "${OUT}"/"${GROUP_2}"_Intergenic_Conditioned \
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
        -r "${REGIONS}" \
        -sites "${OUT}"/Intersect."${GROUP_1}"."${GROUP_2}"_intergenic.txt
fi

#   Estimate joint SFS using realSFS
echo "WRAPPER: realSFS 2dsfs..." >&2
"${ANGSD_DIR}"/misc/realSFS 2dsfs \
    "${OUT}"/"${GROUP_1}"_Intergenic_Conditioned.saf.idx \
    "${OUT}"/"${GROUP_2}"_Intergenic_Conditioned.saf.idx \
    -P "${N_CORES}" \
    > "${OUT}"/2DSFS_Intergenic."${GROUP_1}"."${GROUP_2}".sfs

#   Estimate the Fst using ngsFST
#   Create a directory for ngsFST results
mkdir -P ${OUT}/ngsFST

#   Get number of sites and individuals
N_SITES=`wc -l < "${OUT}"/Intersect."${GROUP_1}"."${GROUP_2}"_intergenic.txt`

#   Convert 2DSFS file for ngsPopGen use
Rscript ${SOURCE}/Wrappers/SFS_to_FST.R ${OUT}/2DSFS_Intergenic."${GROUP_1}"."${GROUP_2}".sfs \
    > ${OUT}/ngsFST/2DSFS_Intergenic.${GROUP_1}.${GROUP_2}.converted.sfs

#   Get the FST
#       This is what needs the new ngsPopGen
#       Make sure to change this path when wrapper is updated
${NGS_POPGEN}/ngsFST \
    -postfiles ${OUT}/"${GROUP_1}"_Intergenic_Conditioned.sfs ${OUT}/"${GROUP_2}"_Intergenic_Conditioned.sfs \
    -priorfile ${OUT}/ngsFST/2DSFS_Intergenic.${GROUP_1}.${GROUP_2}.converted.sfs \
    -nind ${N_IND_1} ${N_IND_2} \
    -nsites ${N_SITES} \
    -outfile ${OUT}/ngsFST/${GROUP_1}.${GROUP_2}.fst
