#!/bin/bash

set -e
set -u
set -o pipefail

#   A script to upgrade angsd-wrapper dependencies

#   Where is angsd-wrapper located?
SOURCE=$1

#   Begin the upgrade process
cd ${SOURCE}/dependencies
ROOT=`pwd`

#   Upgrade HTSLIB
cd ${ROOT}
cd htslib
git checkout -f
git pull
git reset --hard c50bfe05663020c2d8f73ef2c862ee74078d578e
make
make prefix=`pwd` install
HTSLIB_DIR=`pwd`
cd ${ROOT}

#   Upgrade ANGSD
cd ${ROOT}
cd angsd
git checkout -f
git pull
git reset --hard 6d9f209fd33898f9a5520232e726f852c3dee1c5
make HTSDIR="${HTSLIB_DIR}"
cd ${ROOT}

#   Upgrade ngsAdmix
cd ${ROOT}
cd ngsAdmix
rm -rf *
wget http://popgen.dk/software/NGSadmix/ngsadmix32.cpp
g++ ngsadmix32.cpp -O3 -lpthread -lz -o NGSadmix
cd ${ROOT}

#   Upgrade ngsPopGen
cd ${ROOT}
cd ngsPopGen
git checkout -f
git pull
git reset --hard abeabb73b547e067d32d620d6b58a54aad7c0070
make
cd "${ROOT}"

#   Upgrade ngsF
cd "${ROOT}"
git checkout -f
git pull
git reset --hard c39b6ad35c8512d29f09dc4ffd7b6c30afcebd16
make
cd "${ROOT}"
