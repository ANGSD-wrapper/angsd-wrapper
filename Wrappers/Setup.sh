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
wget http://de.iplantcollaborative.org/dl/d/3A541C91-A66A-4651-949D-4E65028C4A2F/iplant.zip
unzip iplant.zip
#       Get rid of the zip file
rm iplant.zip
#       Change into the iplant directory
cd iplant
#       Create a list of sample names
echo "Creating list of sample names..."
rm test_samples.txt
find "$(pwd)" -name "*.bam" | sort > SampleNames.txt
#       Index the BAM files
for i in $(cat SampleNames.txt)
do
    samtools index "$i"
done
#       Index the reference and ancestral sequences
echo "Indexing reference and ancestral sequences..."
find "$(pwd)" -name "*.fa.gz" -exec gzip -d {} \;
find "$(pwd)" -name "*.fa" -exec samtools faidx {} \;
#       Rename the inbreeding coefficients file
echo "Creating a list of inbreeding coefficients..."
mv test_F.txt InbreedingCoefficients.indF
#       Create a regions file
echo "Creating a regions file..."
for i in $(seq 12)
do
    echo "$i": >> regions.txt
done

#   Display final setup message
echo
echo "Test data can be found at $(pwd)"
echo "Please run 'source ~/.bash_profile' to complete installation"
