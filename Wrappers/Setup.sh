#!/bin/bash

set -e
set -u
set -o pipefail

#   Check to see if Git and Wget are installed
if `command -v git > /dev/null 2> /dev/null` && `command -v wget > /dev/null 2> /dev/null` && `command -v samtools > /dev/null 2> /dev/null`
then
    echo "Git, Wget, and SAMTools found"
else
    echo "Cannot find either Git, Wget, or SAMTools, please install and place in your PATH"
    exit 1
fi

#   Append to the angsd-wrapper file
OS=`uname -s`
if [[ "$OS" == "Darwin" ]]
then
    sed -i '' '37i\
    SOURCE='`pwd`'
    ' angsd-wrapper
elif [[ "$OS" == "Linux" ]]
then
    sed -i '37iSOURCE='`pwd`'' angsd-wrapper
else
    echo "Unknown operating system"
    exit 1
fi
echo alias "angsd-wrapper='`pwd -P`/angsd-wrapper'" >> ~/.bash_profile
source ~/.bash_profile

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

#   Download and set up the test data
cd ${ROOT}
cd ..
wget http://de.iplantcollaborative.org/dl/d/3A541C91-A66A-4651-949D-4E65028C4A2F/iplant.zip
unzip iplant.zip
cd iplant
#       Create a list of sample names
rm test_samples.txt
find `pwd` -name "*.bam" | sort > SampleNames.txt
#       Index the BAM files
for i in `cat SampleNames.txt`
do
    samtools index "$i"
done
#       Index the reference and ancestral sequences
find `pwd` -name "*.fa.gz" -exec gzip -d {} \;
find `pwd` -name "*.fa" -exec samtools faidx {} \;
#       Rename the inbreeding coefficients file
mv test_F.txt InbreedingCoefficients.txt
#       Create a regions file
for i in `seq 12`
do
    echo "$i": >> regions.txt
done

echo "Test data can be found at `pwd`"
