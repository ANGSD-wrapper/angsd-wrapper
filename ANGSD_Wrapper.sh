#!/bin/bash

set -e
set -u
set -o pipefail

#   A script that runs loads a config file and runs a wrapper script for ANGSD

function usage() {
    echo -e "\
Usage:  ./ANGSD_Wrapper.sh <wrapper> <config> \n\
where:  <wrapper> is one of SFS
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
    * )
        usage
        ;;
esac