#!/bin/bash

set -e
set -u
set -o pipefail

#   Check to see if Git and Wget are installed
if `command -v git > /dev/null 2> /dev/null` && `command -v wget > /dev/null 2> /dev/null`
then
    echo "Git and Wget found"
else
    echo "Cannot find either Git or Wget, please install and place in your PATH"
    exit 1
fi

#   Make the 'dependencies' directory
mkdir dependencies
cd dependencies
ROOT=`pwd`

#   Install HTSLIB
cd "${ROOT}"
git clone https://github.com/samtools/htslib.git
cd htslib
git reset --hard c50bfe05663020c2d8f73ef2c862ee74078d578e
make
make prefix=`pwd` install
echo export PATH='$PATH':`pwd` >> ~/.bash_profile
source ~/.bash_profile
HTSLIB_DIR=`pwd`
cd "${ROOT}"

#   Install ANGSD
cd "${ROOT}"
git clone https://github.com/ANGSD/angsd.git
cd angsd
git reset --hard 6d9f209fd33898f9a5520232e726f852c3dee1c5
make HTSDIR="${HTSLIB_DIR}"
echo "Path to ANGSD is `pwd`"
echo export PATH='$PATH':`pwd` >> ~/.bash_profile
source ~/.bash_profile
cd "${ROOT}"

#   Install ngsAdmix
cd "${ROOT}"
mkdir ngsAdmix
cd ngsAdmix
wget http://popgen.dk/software/NGSadmix/ngsadmix32.cpp
g++ ngsadmix32.cpp -O3 -lpthread -lz -o NGSadmix
echo "Path to ngsAdmix is `pwd`"
echo export PATH='$PATH':`pwd` >> ~/.bash_profile
source ~/.bash_profile
cd "${ROOT}"

#   Install ngsPopGen
cd "${ROOT}"
git clone https://github.com/mfumagalli/ngsPopGen.git
cd ngsPopGen
git reset --hard abeabb73b547e067d32d620d6b58a54aad7c0070
make
echo "Path to ngsPopGen is `pwd`"
echo export PATH='$PATH':`pwd` >> ~/.bash_profile
source ~/.bash_profile
cd "${ROOT}"

#   Install ngsF
cd "${ROOT}"
git clone https://github.com/fgvieira/ngsF.git
cd ngsF
git reset --hard c39b6ad35c8512d29f09dc4ffd7b6c30afcebd16
make
echo "Path to ngsF is `pwd`"
echo export PATH='$PATH':`pwd` >> ~/.bash_profile
source ~/.bash_profile
cd "${ROOT}"
