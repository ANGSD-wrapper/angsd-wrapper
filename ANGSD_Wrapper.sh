#!/bin/bash

set -e
set -u
set -o pipefail

#   A script that runs loads a config file and runs a wrapper script for ANGSD

function usage() {
    echo -e "\
Usage:  ./ANGSD_Wrapper.sh <wrapper> <config> \n\
where:  <wrapper> is one of SFS, 2DSFS, ABBA_BABA
and:    <config> is the corresponding configuration file \n\
" >&2
    exit 1
}

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
    "ABBA_BABA" )
        #   Run the ABBA-BABA wrapper script
        bash Wrappers/ABBA_BABA.sh "${CONFIG}"
        ;;
    * )
        usage
        ;;
esac