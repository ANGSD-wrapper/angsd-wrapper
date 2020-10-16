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
	if ! $(command -v conda > /dev/null 2> /dev/null); then echo "Please install conda and place in your PATH" >&2 ; exit 1; fi
	
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
        cd "${ROOT}"
        git clone https://github.com/fgvieira/ngsF.git
        cd ngsF
        git reset --hard 807ca7216ab8c3fbd98e628ef1638177d5c752b9
        make
        cd "${ROOT}"

        #   Install ngsAdmix
        cd "${ROOT}"
        mkdir ngsAdmix
        cd ngsAdmix
        wget http://popgen.dk/software/download/NGSadmix/ngsadmix32.cpp
        g++ ngsadmix32.cpp -O3 -lpthread -lz -o NGSadmix
        cd "${ROOT}"

        #   Install ngsPopGen
        cd "${ROOT}"
        git clone https://github.com/mfumagalli/ngsPopGen.git
        cd ngsPopGen
        git reset --hard bbd73d5caa660f28111c69eefca3230ded4a97ac
        make
        cd "${ROOT}"
        echo

	#   Installing ANGSD and its dependencies
	conda env create -f /path/to/angsd-wrapper/Wrappers/angsd-wrapper.yml
	
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
        wget --no-check-certificate --output-document=Example_Data.tar.bz2 https://ndownloader.figshare.com/files/5282197
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
        for inbreeding in $(find ${EXAMPLE_DIR} -name "*Inbreeding.txt")
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
