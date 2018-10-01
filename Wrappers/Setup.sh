#!/bin/bash

set -e
set -u

#   Arguments
declare -a args=("$@")

setup_routine="${args[0]}" # Which routine are we running
SOURCE="${args[1]}" # Where is ANGSD-wrapper?
BASESOURCE="${args[2]}"

case "${setup_routine}" in
    "dependencies" )
        if [[ -x $(command -v singularity) ]]; then
                cd "${SOURCE}"/Wrappers
                if [[ ! -f carte731-angsd-wrapper-update-master-latest.simg ]]; then
                    echo -e "Installing CentOS-7 image for stable installation."
                    singularity pull shub://carte731/angsd-wrapper-update
                fi
                # ./angsd_singularity.simg "${SOURCE}" "${BASESOURCE}"
                ./carte731-angsd-wrapper-update-master-latest.simg "${SOURCE}" "${BASESOURCE}"
                echo -e "Angsd-Wrapper has been installed.\n"
        else
                echo -e "Please install or module load Singularity.\n"
                exit 1
        fi
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
