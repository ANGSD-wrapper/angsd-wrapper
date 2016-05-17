#!/bin/bash

set -e
set -u

#   Arguments
setup_routine=$1 # Which routine are we running
SOURCE=$2 # Where is ANGSD-wrapper?

#	Download and install SAMTools 1.3
function installSAMTools() {
	wget https://github.com/samtools/samtools/releases/download/1.3/samtools-1.3.tar.bz2 # Download SAMTools
	tar -xvjf samtools-1.3.tar.bz2 # Extract the tarball
	rm -f samtools-1.3.tar.bz2 # Get rid of the tarball
	cd samtools-1.3 # Change into the SAMTools directory
	./configure --prefix=$(pwd) # Configure the installation process, setting the install directory to be here
	make # Compile the code
	make install # Install SAMTools
	echo "export PATH=$(pwd):"'${PATH}' >> ~/.bash_profile # Add the path to bash_profile
}

#   Export the function
export -f installSAMTools

case "${setup_routine}" in
    "dependencies" )
        #   Check to see if Git and Wget are installed
        if ! $(command -v git > /dev/null 2> /dev/null); then echo "Please install Git and place in your PATH" >&2 ; exit 1; fi
        if ! $(command -v wget > /dev/null 2> /dev/null); then echo "Please install Wget and place in your PATH" >&2 ; exit 1; fi
        #   Let angsd-wrapper be run from anywhere
        echo alias "angsd-wrapper='`pwd -P`/angsd-wrapper'" >> ~/.bash_profile
        #   Make the 'dependencies' directory
        cd "${SOURCE}"
        mkdir dependencies
        cd dependencies
        ROOT=$(pwd)
        #   Check for SAMTools. If not found, install it
        if ! $(command -v samtools > /dev/null 2> /dev/null); then cd "${ROOT}"; installSAMTools; source ~/.bash_profile;cd "${ROOT}"; fi
        #   Install ngsF
        if [[ $(uname) == "Linux" ]] # Are we on Linux?
        then # If so, we can use the normal ngsF
            cd "${ROOT}"
            git clone https://github.com/fgvieira/ngsF.git
            cd ngsF
            git reset --hard c39b6ad35c8512d29f09dc4ffd7b6c30afcebd16
            make
            cd "${ROOT}"
        elif [[ $(uname) == "Darwin" ]] # If we're on Mac OS X
        then # We need a version with a compalition routine that works for Mac OS X
            cd "${ROOT}"
            git clone https://github.com/mojaveazure/ngsF.git
            cd ngsF
            bash install.sh
            cd "${ROOT}"
        else # I don't know how other BSDs works with ngsF, so yeah...
            echo "Failed to determine operating system. If not using a Windows-based machine, please file an issue and let us know!"
            exit 1
        fi
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
        echo
        #   Display final setup message
        echo "Please run 'source ~/.bash_profile' to complete installation"
        ;;
    "data" )
        #   Check for depenent programs
        if ! $(command -v wget > /dev/null 2> /dev/null); then echo "Please install Wget and place in your PATH" >&2 ; exit 1; fi
        if ! $(command -v samtools > /dev/null 2> /dev/null); then echo "Please install SAMTools and place in your PATH" >&2; exit 1; fi
        #   Download and set up the test data
        if [[ ${SOURCE} == '.' ]]; then SOURCE=$(pwd -P); fi
        cd "${SOURCE}"
        if [[ -d "Example_Data" ]]; then rm -rf Example_Data/; fi # Remove any old example datasets
        wget --no-check-certificate --output-document=Example_Data.tar.bz2 https://ndownloader.figshare.com/files/3667101
        tar -xvjf Example_Data.tar.bz2
        rm Example_Data.tar.bz2
        EXAMPLE_DIR="${SOURCE}/Example_Data"
        #       Change into the example data directory
        cd "${EXAMPLE_DIR}"
        #       Create lists of sample names
        echo "Creating sample lists..." >&2
        find "${EXAMPLE_DIR}"/Maize/ -name "*.bam" | sort > Maize/Maize_Samples.txt
        find "${EXAMPLE_DIR}"/Mexicana -name "*.bam" | sort > Mexicana/Mexicana_Samples.txt
        find "${EXAMPLE_DIR}"/Teosinte -name "*.bam" | sort > Teosinte/Teosinte_Samples.txt
        #   Make sure all inbreeding files are named "*.indF"
        for inbreeding in $(find ${EXAMPLE_DIR} -name "*Inbreeding*")
        do
            BASE=$(basename ${inbreeding} | cut -f 1 -d '.')
            DIR=$(dirname ${inbreeding})
            mv ${inbreeding} ${DIR}/${BASE}.indF
        done
        #       Index the reference and ancestral sequences
        echo "Indexing reference and ancestral sequences..." >&2
        cd Sequences
        find "${EXAMPLE_DIR}" -name "*.fa" -exec samtools faidx {} \;
        #   Index the BAM files
        echo "Indexing the BAM files..." >&2
        cd "${EXAMPLE_DIR}"/Maize # Maize samples
        for i in $(cat Maize_Samples.txt); do samtools index "$i"; done
        cd "${EXAMPLE_DIR}"/Mexicana # Mexicana samples
        for i in $(cat Mexicana_Samples.txt); do samtools index "$i"; done
        cd "${EXAMPLE_DIR}"/Teosinte # Teosinte samples
        for i in $(cat Teosinte_Samples.txt); do samtools index "$i"; done
        cd "${EXAMPLE_DIR}"
        #   Write the example configuration files
        echo "Writing example configuration files..." >&2
        source "${SOURCE}/Wrappers/writeConfigs.sh"
        CONFIG_DIR=$(writeConfigs ${EXAMPLE_DIR})
        #   Finished
        echo
        echo "Test data can be found at ${EXAMPLE_DIR}"
        echo "Example configuration files can be found at ${CONFIG_DIR}"
        ;;
esac
