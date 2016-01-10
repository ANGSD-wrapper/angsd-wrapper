#!/bin/bash

set -e
set -u

#   Check to see if Git, Wget and SAMTools are installed
if `command -v git > /dev/null 2> /dev/null` && `command -v wget > /dev/null 2> /dev/null` && `command -v samtools > /dev/null 2> /dev/null`
then
    echo "Git, Wget, and SAMTools found"
else
    echo "Cannot find either Git, Wget, or SAMTools, please install and place in your PATH"
    exit 1
fi

#   Let angsd-wrapper be run from anywhere
echo alias "angsd-wrapper='`pwd -P`/angsd-wrapper'" >> ~/.bash_profile

#   Make the 'dependencies' directory
mkdir dependencies
cd dependencies
ROOT=$(pwd)

#   Install HTSLIB
cd "${ROOT}"
git clone https://github.com/samtools/htslib.git
cd htslib
git reset --hard 306664a776435a1c86d7263f16deb43b30db55fd
make
make prefix=`pwd` install
HTSLIB_DIR=`pwd`
cd "${ROOT}"

#   Install ANGSD
cd "${ROOT}"
git clone https://github.com/ANGSD/angsd.git
cd angsd
git reset --hard 8b89ba421e97c125f474b2217b710f178c27a51e
make HTSDIR="${HTSLIB_DIR}"
cd "${ROOT}"

#   Install ngsAdmix
cd "${ROOT}"
mkdir ngsAdmix
cd ngsAdmix
wget http://popgen.dk/software/NGSadmix/ngsadmix32.cpp
g++ ngsadmix32.cpp -O3 -lpthread -lz -o NGSadmix
cd "${ROOT}"

#   Install ngsPopGen
cd "${ROOT}"
git clone https://github.com/mfumagalli/ngsPopGen.git
cd ngsPopGen
git reset --hard abeabb73b547e067d32d620d6b58a54aad7c0070
make
cd "${ROOT}"

#   Install ngsF
cd "${ROOT}"
git clone https://github.com/fgvieira/ngsF.git
cd ngsF
git reset --hard c39b6ad35c8512d29f09dc4ffd7b6c30afcebd16
make
cd "${ROOT}"

#   Download and set up the test data
cd "${ROOT}"
cd ..
wget 'http://de.iplantcollaborative.org/dl/d/8D7AF180-508E-4639-941F-946FA14C0120/Example_Data.tar.bz2'
tar -xvjf Example_Data.tar.bz2
#       Change into the example data directory
cd Example_Data
#       Create lists of sample names
echo "Creating sample lists..." >&2
find $(pwd)/Maize/ -name "*.bam" | sort > Maize/Maize_Samples.txt
find $(pwd)/Mexicana -name "*.bam" | sort > Mexicana/Mexicana_Samples.txt
find $(pwd)/Teosinte -name "*.bam" | sort > Teosinte/Teostinte_Samples.txt
#       Index the reference and ancestral sequences
echo "Indexing reference and ancestral sequences..." >&2
cd Sequences
find $(pwd) -name "*.fa" -exec samtools faidx {} \;
#   Index the BAM files
echo "Indexing the BAM files..." >&2
cd ../Maize # Maize samples
for i in `cat Maize_Samples.txt`; do samtools index "$i"; done
cd ../Mexicana # Mexicana samples
for i in `cat Mexicana_Samples.txt`; do samtools index "$i"; done
cd ../Teosinte # Teosinte samples
for i in `cat Teosinte_Samples.txt`; do samtools index "$i"; done

#   Display final setup message
cd ..
echo
echo "Test data can be found at $(pwd)"
echo "Please run 'source ~/.bash_profile' to complete installation"
