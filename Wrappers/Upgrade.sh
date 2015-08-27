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
git reset --hard 306664a776435a1c86d7263f16deb43b30db55fd
make
make prefix=`pwd` install
HTSLIB_DIR=`pwd`
cd ${ROOT}

#   Upgrade ANGSD
cd ${ROOT}
cd angsd
git checkout -f
git pull
git reset --hard 8b89ba421e97c125f474b2217b710f178c27a51e
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
