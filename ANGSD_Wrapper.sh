#!/bin/bash

set -e
set -u
set -o pipefail

#   A script that runs loads a config file and runs a wrapper script for ANGSD

function usage() {
    echo -e "\
Usage:  ./ANGSD_Wrapper.sh <wrapper> <config> \n\
where:  <wrapper> is one of SFS, 2DSFS, ABBA_BABA, ANC_SEQ, Genotypes, Thetas, ngsF
and:    <config> is the corresponding configuration file \n\
\n\
The following is a list of calls and what they do: \n\
    SFS         Site Frequency Spectrum \n\
    2DSFS      2D Site Frequency Spectrum \n\
    ngsFST      FST calculations, THIS NEEDS TO BE \n\
                    RUN AFTER 2D SFS
    ABBA_BABBA  ABBA_BABBA Test \n\
    ANC_SEQ     Extract Ancestral Sequence from BAM file \n\
    Genotypes   Genotype Calling \n\
    Thetas      Estimate thetas and perform neutrality test statistics \n\
    ngsF        Use ngsF to calculate per-individual inbreeding coefficients \n\
    ngsAdmix    Run ngsAdmix
\n\
For an interactive help, type 'help me' without the quotes \n\
" >&2
    exit 1
}

#   Check to see if Wrappers directory is one level below this
if ! [[ -d 'Wrappers' ]]
then
    echo "Cannot find 'Wrappers' directory, exiting" >&2
    exit 1
fi

#   Check to see if there are the correct number of arguments
if [ "$#" -lt 2 ]
then
    usage
fi

WRAPPER="$1"
CONFIG="$2"

case "${WRAPPER}" in
    "SFS" )
        #   Run the site frequency spectrum wrapper script
        bash Wrappers/Site_Frequency_Spectrum.sh "${CONFIG}"
        ;;
    "2DSFS" )
        #   Run the 2D site frequency spectrum wrapper script
        bash Wrappers/2D_Site_Frequency_Spectrum.sh "${CONFIG}"
        ;;
    "ngsFST" )
        #   Run the ngsFST wrapper script
        bash Wrappers/ngsFST.sh "${CONFIG}"
        ;;
    "ABBA_BABA" )
        #   Run the ABBA-BABA wrapper script
        bash Wrappers/ABBA_BABA.sh "${CONFIG}"
        ;;
    "ANC_SEQ" )
        #   Run the ancestral sequence wrapper script
        bash Wrappers/Ancestral_Sequence.sh "${CONFIG}"
        ;;
    "Genotypes" )
        #   Run the genotypes wrapper script
        bash Wrappers/Genotypes.sh "${CONFIG}"
        ;;
    "Thetas" )
        #   Run the thetas wrapper script
        bash Wrappers/Thetas.sh "${CONFIG}"
        ;;
    "ngsF" )
        #   Run the ngsF wrapper script
        bash Wrappers/ngsF.sh "${CONFIG}"
        ;;
    "ngsAdmix" )
        #   Run the ngsAdmix wrapper script
        bash Wrappers/ngsAdmix.sh "${CONFIG}"
        ;;
    "help" )
        bash Wrappers/Help.sh
        ;;
    * )
        usage
        ;;
esac
