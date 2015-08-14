#!/bin/bash

set -e
set -u
set -o pipefail


if `command -v git > /dev/null 2> /dev/null` && `command -v wget > /dev/null 2> /dev/null`
then
    echo "Git and Wget found"
else
    echo "Cannot find either Git or Wget, please install and place in your PATH"
    exit 1
fi


function i_angsd() {
    ROOT="$1"
    HTSLIB_DIR="$2"
    cd "${ROOT}"
    git clone https://github.com/ANGSD/angsd.git
    cd angsd
    make HTSDIR="${HTSLIB_DIR}"
    echo "Path to ANGSD is `pwd`"
    echo "export PATH=$PATH:`pwd`" >> ~/.bash_profile
    source ~/.bash_profile
    cd "${ROOT}"
}

export -f i_angsd

function i_htslib() {
    ROOT="$1"
    cd "${ROOT}"
    git clone https://github.com/samtools/htslib.git
    cd htslib
    autoconf
    ./configure --prefix=`pwd`
    make
    make install
    echo "export PATH=$PATH:`pwd`" >> ~/.bash_profile
    source ~/.bash_profile
    cd "${ROOT}"
}

export -f i_htslib

function i_ngsAdmix() {
    ROOT="$1"
    cd "${ROOT}"
    mkdir ngsAdmix
    cd ngsAdmix
    wget popgen.dk/software/NGSadmix/ngsadmix32.cpp
    g++ ngsadmix32.cpp -O3 -lpthread -lz -o NGSadmix
    echo "Path to ngsAdmix is `pwd`"
    echo "export PATH=$PATH:`pwd`" >> ~/.bash_profile
    source ~/.bash_profile
    cd "${ROOT}"
}

export -f i_ngsAdmix

function i_ngsPopGen() {
    ROOT="$1"
    cd "${ROOT}"
    git clone https://github.com/mfumagalli/ngsPopGen.git
    cd ngsPopGen
    make
    echo "Path to ngsPopGen is `pwd`"
    echo "export PATH=$PATH:`pwd`" >> ~/.bash_profile
    source ~/.bash_profile
    cd "${ROOT}"
}

export -f i_ngsPopGen

function i_ngsF() {
    ROOT="$1"
    cd "${ROOT}"
    git clone https://github.com/fgvieira/ngsF.git
    cd ngsF
    make
    make test
    echo "Path to ngsF is `pwd`"
    echo "export PATH=$PATH:`pwd`" >> ~/.bash_profile
    source ~/.bash_profile
    cd "${ROOT}"
}

export -f i_ngsF

function Usage() {
    echo -e "\
Usage:  ./installer.sh program root_directory [htslib_directory] \n\
where:  program is one of the following \n\
            angsd       installs angsd \n\
            ngsAdmix    installs ngsAdmix \n\
            ngsPopGen   installs ngsPopGen \n\
            ngsF        installs ngsF \n\
            all         installs all of the above \n\
\n\
        root_direcotry is the directory in which to install everythin \n\
            defaults to `pwd` \n\
\n\
        [htslib_directory] is the directory where htslib is installed \n\
            if this is not specified, htslib will be installed automatically \n\
"
    exit 1
}

export -f Usage

if [ "$#" -lt 2 ]
then
    Usage
fi

PROG=$1
ROOTDIR=${2:-`pwd`}

case "$PROG" in
    "angsd" )
        if [ -z "$3" ]
        then
            i_htslib "${ROOTDIR}"
            HTSDIR="${ROOTDIR}"/htslib
        else
            HTSDIR="$3"
        fi
        i_angsd "${ROOTDIR}" "${HTSDIR}"
        ;;
    "ngsAdmix" )
        i_ngsAdmix "${ROOTDIR}"
        ;;
    "ngsPopGen" )
        i_ngsPopGen "${ROOTDIR}"
        ;;
    "ngsF" )
        i_ngsF "${ROOTDIR}"
        ;;
    "all" )
        if [ -z "$3" ]
        then
            i_htslib "${ROOTDIR}"
            HTSDIR="${ROOTDIR}"/htslib
        else
            HTSDIR="$3"
        fi
        i_angsd "${ROOTDIR}" "${HTSDIR}"
        i_ngsAdmix "${ROOTDIR}"
        i_ngsPopGen "${ROOTDIR}"
        i_ngsF "${ROOTDIR}"
        ;;
    * )
        Usage
        ;;
esac



