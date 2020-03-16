#!/bin/bash

set -e
set -o pipefail


#   A function to write Common_Config
function createCommon(){
    local path=$1 # Base path for sample list, inbreeding coefficients, and regions file
    local samples=$2 # Name of sample list
    local inbreeding=$3 # Name of inbreeding coefficients
    local ancestral=$4 # Full path to ancestral sequence
    local reference=$5 # Full path to reference sequence
    local project=$6 # Name of project
    local scratch=$7 # Full path to scratch directory
    local regions=${8} # Name of regions file
    local outdirectory=${9} # Where are we putting our
    local outname="${outdirectory}/${project}_Common_Config" # Create a name for the configuration file
    #   Write Common_Config
    echo -e "#!/bin/bash\nset -e \nset -o pipefail \n\n#   A simple script to hold common variables for ANGSD-wrapper\n\n#   Define a list of samples\nSAMPLE_LIST=${path}/${samples}\n\n#   Define a list of inbreeding coefficients\n#	This should end in '_.indF'\nSAMPLE_INBREEDING=${path}/${inbreeding}\n\n#   Ancestral sequence\nANC_SEQ=${ancestral}\n\n#   Reference sequence\nREF_SEQ=${reference}\n\n#   Name the project\nPROJECT=${project}\n\n#   Where do we put the outfiles?\n    #   Note, the final outdirectory will be\n    #   \${SCRATCH}/\${PROJECT}/<name_of_program/>\nSCRATCH=${scratch}\n\n#   Define the region being looked at\n#       Optional, but ANGSD is expensive to run without specifying regions to look at\nREGIONS=${path}/${regions}\n\n#   Set common parameters for all methods\nUNIQUE_ONLY=0\nMIN_BASEQUAL=20\nBAQ=1\nMIN_IND=1\nGT_LIKELIHOOD=2\nMIN_MAPQ=30\nN_CORES=32\nDO_MAJORMINOR=1\nDO_GENO=32\nDO_MAF=1\nDO_POST=1\n" > "${outname}"
    echo "${outname}" # Return the name of Common_Config
}

#   Export the function
export -f createCommon

#   A function to write a configuration file for Site Frequency Spectrum
function createSFS() {
    local common=$1 # Where is Common_Config?
    local project=$2 # What are we calling our configuration file?
    local outdirectory=$3 # Where are we putting our configuration file?
    local outname="${outdirectory}/${project}_Site_Frequency_Spectrum_Config" # Create a name for the configuration file
    #   Write the configuration file for SFS
    echo -e "#!/bin/bash\n\nset -e\nset -u\nset -o pipefail\n\n#   A simple script to hold variables for the Site Frequency Spectrum\n#   Are you using the Common_Config file?\n#       If so, where is it?\nCOMMON=${common}\n\n##############################################################################################\n#   If we aren't using the Common_Config file, specify these variables\n#   If Common_Config is specified, leave these blank\n#   Define a list of samples\nSAMPLE_LIST=\n\n#   Define a list of inbreeding coefficients\nSAMPLE_INBREEDING=\n\n#   Ancestral and Reference sequences\nANC_SEQ=\nREF_SEQ=\n\n#   Name the project\nPROJECT=\n\n#   Where do we put the outfiles?\n    #   Note, the final outdirectory will be\n    #   \${SCRATCH}/\${PROJECT}/SFS\nSCRATCH=\n\n#   Define the region being looked at\n#       Optional, but ANGSD is expensive to run without specifying regions to look at\nREGIONS=\n\n#   Parameters that are specified in Common_Config\nUNIQUE_ONLY=0\nMIN_BASEQUAL=20\nBAQ=1\nMIN_IND=1\nGT_LIKELIHOOD=2\nMIN_MAPQ=30\nN_CORES=32\nDO_MAJORMINOR=1\nDO_GENO=32\nDO_MAF=1\nDO_POST=1\n\n##############################################################################################\n\n#   Site Frequency Spectrum Parameters\n#       Listed below are the defaults, please modify for your samples\n#       Generate site allele frequencies\nDO_SAF=2\n#       Overwrite any previously generated results\nOVERRIDE=true\n\n#   For advanced users who want to change arguments used by ANGSD (i.e. which SAF method is used)\n#       Expected format is the same as how the flags and arguments are written on the command line:\n#       '-flag1 arg1 -flag2 arg2 ...' \nADVANCED_ARGS=''" > "${outname}"
}

#   Export the function
export -f createSFS

#   A function to write a configuration file for Thetas Estimations
function createThetas() {
    local common=$1 # Where is Common_Config
    local pest=$2 # Where is the pest file
    local project=$3 # What are we calling our configuration file?
    local outdirectory=$4 # Where are we putting our configuration file?
    local outname="${outdirectory}/${project}_Thetas_Config" # Create a name for the configuration file
    #   Write the configuration file for Thetas
    echo -e "#!/bin/bash\n\nset -e\nset -u\nset -o pipefail\n\n#   A simple script to hold variables for the Estimation of Thetas\n#   Are you using the Common_Config file?\n#       If so, where is it?\nCOMMON=${common}\n\n##############################################################################################\n#   If we aren't using the Common_Config file, specify these variables\n#   If Common_Config is specified, leave these blank\n#   Define a list of samples\nSAMPLE_LIST=\n\n#   Define a list of inbreeding coefficients\nSAMPLE_INBREEDING=\n\n#   Ancestral and Reference sequences\nANC_SEQ=\nREF_SEQ=\n\n#   Name the project\nPROJECT=\n\n#   Where do we put the outfiles?\n    #   Note, the final outdirectory will be\n    #   \${SCRATCH}/\${PROJECT}/Thetas\nSCRATCH=\n\n#   Define the region being looked at\n#       Optional, but ANGSD is expensive to run without specifying regions to look at\nREGIONS=\n\n#   Parameters that are specified in Common_Config\nUNIQUE_ONLY=0\nMIN_BASEQUAL=20\nBAQ=1\nMIN_IND=1\nGT_LIKELIHOOD=2\nMIN_MAPQ=30\nN_CORES=32\nDO_MAJORMINOR=1\nDO_MAF=1\n\n##############################################################################################\n\n#   Pest file\n#       This is the output from the site frequency spectrum wrapper\n#       This should end in '_DerivedSFS'\nPEST=${pest}\n\n#   Thetas Parameters\n#       Listed below are the defaults, please modify for your samples\nDO_SAF=2\nDO_THETAS=1\nOVERRIDE=true\nSLIDING_WINDOW=false\nWIN=50000\nSTEP=10000\n\n#   For advanced users who want to change arguments used by ANGSD (i.e. which SAF method is used)\n#       Expected format is the same as how the flags and arguments are written on the command line:\n#       '-flag1 arg1 -flag2 arg2 ...' \nADVANCED_ARGS=''" > "${outname}"
}

#   Export the function
export -f createThetas

#   A function to write a configuration file for Genotype Likelihoods
function createGenotypes() {
    local common="$1" # Where is Common_Config?
    local project="$2" # What are we calling our configuration file?
    local outdirectory="$3" # Where are we putting our configuration file?
    local outname="${outdirectory}/${project}_Genotypes_Config" # Create a name for the configuration file
    echo -e "#!/bin/bash\n\nset -e\nset -u\nset -o pipefail\n\n#   A simple script to hold variables for the Site Frequency Spectrum\n#   Are you using the Common_Config file?\n#       If so, where is it?\nCOMMON=${common}\n\n##############################################################################################\n#   If we aren't using the Common_Config file, specify these variables\n#   If Common_Config is specified, leave these blank\n#   Define a list of samples\nSAMPLE_LIST=\n\n#   Define a list of inbreeding coefficients\nSAMPLE_INBREEDING=\n\n#   Name the project\nPROJECT=\n\n#   Where do we put the outfiles?\n    #   Note, the final outdirectory will be\n    #   \${SCRATCH}/\${PROJECT}/GenotypeLikelihoods\nSCRATCH=\n\n#   Define the region being looked at\n#       Optional, but ANGSD is expensive to run without specifying regions to look at\nREGIONS=\n\n#   Set common parameters for all methods\n#       Use only uniquely-mapped reads\nUNIQUE_ONLY=0\n#       Set the minimum base quality\nMIN_BASEQUAL=20\n#       Calculate base alignment quality\nBAQ=1\n#       Set the minimum number of individuals required\nMIN_IND=1\n#       Calculate genotype likelihoods\nGT_LIKELIHOOD=2\n#       Set the minimum mapping quality for a base to be used\nMIN_MAPQ=30\n#       Set the number of threads to be used\nN_CORES=32\n#       Determine major and minor alleles\nDO_MAJORMINOR=1\n#       Call genotypes from genotype likelihoods\nDO_GENO=32\n#       Calculate allele frequencies\nDO_MAF=1\n#       Calculate the posterior probability\nDO_POST=1\n\n##############################################################################################\n\n#   Genotypes Parameters\n#       Listed below are the defaults, please modify for your samples\n#       Set the minimum posterior value for calling genotypes\nPOST_CUTOFF=0.95\n#       Set the maximum p-value for polymorphic sites\nSNP_PVAL=1e-6\n#       Output genotype likelihood frequency file\nDO_GLF=2\n\n#   For advanced users who want to change arguments used by ANGSD (i.e. which SAF method is used)\n#       Expected format is the same as how the flags and arguments are written on the command line:\n#       '-flag1 arg1 -flag2 arg2 ...' \nADVANCED_ARGS=''" > "${outname}"
}

#   Export the function
export -f createGenotypes

#   A function to write a configuration file for Admixture Analysis
function createAdmixture() {
    local common=$1 # Where is Common_Config?
    local likelihood=$2 # Where is the likelihood file
    local project=$3 # What are we calling our configuration file?
    local outdirectory=$4 # Where are we putting our configuration file?
    local outname="${outdirectory}/${project}_Admixture_Config" # Create a name for the configuration file
    #   Write out the configuration file for Admixture
    echo -e "#!/bin/bash\n\nset -e\nset -u\nset -o pipefail\n\n#   A simple script to hold the varialbes for the NGS Admixture\n#   Are you using the Common_Config file?\n#       If so, where is it?\nCOMMON=${common}\n\n##############################################################################################\n#   If we aren't using the Common_Config file, specify these variables\n#   If Common_Config is specified, leave these blank\n#   Name the project\nPROJECT=\n\n#   Where do we put the outfiles?\n    #   Note, the final outdirectory will be\n    #   \${SCRATCH}/\${PROJECT}/Admixture\nSCRATCH=\n\n#   Parameters that are specified in Common_Config\nN_CORES=32\n\n##############################################################################################\n\n#   The Likelihood file\n#       This is the .beagle.gz file from the Site Frequency Spectrum\nLIKELIHOOD=${likelihood}\n\n#   ngsAdmix Parameters\n#       Listed below are the defaults, please modify for your samples\nK=5\nMIN_MAF=0.05\nTOLERANCE=0.01\n" > "${outname}"
}

#   Export the function
export -f createAdmixture

#   A function to write a configuration file for Principal Component Analysis
function createPCA() {
    local common=$1 # Where is Common_Config?
    local project=$2 # What are we calling our configuration file?
    local outdirectory=$3 # Where are we putting our configuration file?
    local outname="${outdirectory}/${project}_Principal_Component_Analysis_Config" # Create a name for the configuration file
    #   Write out the configuration file for PCA
    echo -e "#!/bin/bash\n\nset -e\nset -u\nset -o pipefail\n\n#   A simple script to hold variables for the Principal Component Analysis\n#   Are you using the Common_Config file?\n#       If so, where is it?\nCOMMON=${common}\n\n##############################################################################################\n#   If we aren't using the Common_Config file, specify these variables\n#   If Common_Config is specified, leave these blank\n#   Define a list of samples\nSAMPLE_LIST=\n\n#   Name the project\nPROJECT=\n\n#   Where do we put the outfiles?\n    #   Note, the final outdirectory will be\n    #   \${SCRATCH}/\${PROJECT}/PCA\nSCRATCH=\n\n#   Region being looked at?\nREGIONS=\n##############################################################################################\n\n#   Principal Component Analysis Parameters\n#       Listed below are the defaults, please modify for your samples\nDO_MAF=2\nDO_MAJORMINOR=1\nDO_GENO=32\nDO_POST=1\nN_CORES=8\nNORM=0\nCALL=0\nGT_LIKELIHOOD=2\nN_SITES=100000\n\n#   For advanced users who want to change arguments used by ANGSD (i.e. which SAF method is used)\n#       Expected format is the same as how the flags and arguments are written on the command line:\n#       '-flag1 arg1 -flag2 arg2 ...' \nADVANCED_ARGS=''" > "${outname}"
}

#   Export the function
export -f createPCA

#   A function to write a configuration file for Abbababa
function createAbbababa() {
   local common=$1 # Where is Common_Config?
   local project=$2 # What are we calling our configuration file?
   local outdirectory=$3 # Where are we putting our configuration file?
   local outname="${outdirectory}/${project}_Abbababa_Config" # Create a name for the configuration file
   #   Write out the configuration file for Abbababa
   echo -e "#!/bin/bash\n\nset -e\nset -u\nset -o pipefail\n\n#   A simple script to hold variables for ABBA BABA\n#   Are you using the Common_Config file?\n#       If so, where is it?\nCOMMON=${common}\n\n##############################################################################################\n#   If we aren't using the Common_Config file, specify these variables\n#   If Common_Config is specified, leave these blank\n#   Define a list of samples\nSAMPLE_LIST=\n\n#   Name the project\nPROJECT=\n\n#   Where do we put the outfiles?\n#      Note, the final outdirectory will be\n#      \${SCRATCH}/\${PROJECT}/ABBABABA\nSCRATCH=\n\n#   Region being looked at?\n#       Optional, but ANGSD is expensive to run without specifying regions to look at\nREGIONS=\n\n#   Parameters that are specified in Common_Config\n#       Use only uniquely-mapped reads\nUNIQUE_ONLY=0\n#       Set the minimum base quality\nMIN_BASEQUAL=20\n#       Set the minimum number of individuals required\nMIN_IND=1\n#       Set the minimum mapping quality for a base to be used\nMIN_MAPQ=30\n#       Set the number of threads to be used\nN_CORES=32\\n\n##############################################################################################\n\n#   Fasta file to be used as an outgroup\n#       This can be generated using the ANC_SEQ function of angsd-wrapper\nOUTGROUP=\n\n#   ABBA BABA Parameters\n#       Listed below are the defaults, please modify for your samples\n#       Count allele frequencies\nDO_COUNTS=1\n#       Perform Abbababa analysis\nDO_ABBABABA=1\n#       Remove transitions from sample\nREMOVE_TRANS=0\n#       Set the size for each block\nBLOCKSIZE=1000\n\n#   For advanced users who want to change arguments used by ANGSD (i.e. which SAF method is used)\n#       Expected format is the same as how the flags and arguments are written on the command line:\n#       '-flag1 arg1 -flag2 arg2 ...' \nADVANCED_ARGS=''" > "${outname}"

}

#   Export the function
export -f createAbbababa


#   A function to write a configuration file for Ancestral Sequence
function createAncestral() {
    local common=$1 # Where is Common_Config?
    local project=$2 # What are we calling our configuration file?
    local outdirectory=$3 # Where are we putting our configuration file?
    local anc_bam=$4 # Where is our ancestral bam file?
    local outname="${outdirectory}/${project}_Ancestral_Sequence_Config" # Create a name for the configuration file
    #   Write out the configuration file for Abbababa
    echo -e "#!/bin/bash\n\nset -e\nset -u\nset -o pipefail\n\n#   A simple script to hold variables for generating an ancestral fasta file\n\n#   This script does NOT utilize the Common_Config file\n\n#   Where is the ancestral BAM file?\nANC_BAM=${anc_bam}\n\n#   What should we call the output file?\n#       Defaults to the same name as the ancestral BAM file\nOUT_NAME=\`basename '\${ANC_BAM}' .bam\`\n\n#   Where should we put the output file?\n#       Defaults to the same directory as the ancestral BAM file.\nOUT_DIR=\`dirname '\${ANC_BAM}'\`\n\n#   Full path to the output file\nOUT=\${OUT_DIR}/\${OUT_NAME}\n\n#   Ancestral Sequence Parameters\n#       Listed below are the defaults, please modify for your samples\n#       Extract FASTA sequence from BAM file\nDO_FASTA=1\n#       Count allele frequencies\n#       If DO_FASTA is 2, DO_COUNTS must be 1\n#       Otherwise, DO_COUNTS can be any other legal value\nDO_COUNTS=0\n\n#   For advanced users who want to change arguments used by ANGSD (i.e. which SAF method is used)\n#       Expected format is the same as how the flags and arguments are written on the command line:\n#       '-flag1 arg1 -flag2 arg2 ...' \nADVANCED_ARGS=''" > "${outname}"
}

#   Export the function
export -f createAncestral


#   A function to write a configuration file for FST
function createFST() {
    local common=$1 # Where is Common_Config?
    local project=$2 # What are we calling our configuration file?
    local outdirectory=$3 # Where are we putting our configuration file?
    local ancestral=$4 # Where is our ancestral sequence?
    local reference=$5 # Where is our reference sequence?
    local outname="${outdirectory}/${project}_FST_Config" # Create a name for the configuration file
    #   Write out the configuration file for Abbababa
    echo -e "#!/bin/bash\n\nset -e\nset -u\nset -o pipefail\n\n#   A simple script to hold variables for the FST\n#   Are you using the Common_Config file?\n#       If so, where is it?\nCOMMON=${common}\n\n##############################################################################################\n#   If we aren't using the Common_Config file, specify these variables\n#   If Common_Config is specified, leave these blank\n#   Ancestral and Reference sequences\nANC_SEQ=${ancestral}\nREF_SEQ=${reference}\n\n#   Name the project\nPROJECT=${project}\n\n#   Where do we put the outfiles?\n    #   Note, the final outdirectory will be\n    #   \${SCRATCH}/\${PROJECT}/Fst\nSCRATCH=\n\n#   Region being looked at?\nREGIONS=\n\n#   Parameters that are specified in Common_Config\n#       Use only uniquely-mapped reads\nUNIQUE_ONLY=0\n#       Set the minimum base quality\nMIN_BASEQUAL=20\n#       Calculate base alignment quality\nBAQ=1\n#       Calculate genotype likelihoods\nGT_LIKELIHOOD=2\n#       Set the minimum mapping quality for a base to be used\nMIN_MAPQ=30\n#       Set the number of threads to be used\nN_CORES=32\n#       Determine major and minor alleles\nDO_MAJORMINOR=1\n#       Call genotypes from genotype likelihoods\nDO_GENO=32\n#       Calculate allele frequencies\nDO_MAF=1\n#       Calculate the posterior probability\nDO_POST=1\n\n##############################################################################################\n\n#   What is group 1?\nGROUP_1=\n\n#   Sample list for group 1\nG1_SAMPLE_LIST=\n\n#   Inbreeding coefficients for group 1\nG1_INBREEDING=\n\n#   What is group 2?\nGROUP_2=\n\n#   Sample list for group 2\nG2_SAMPLE_LIST=\n\n#   Inbreeding coefficients for group 2\nG2_INBREEDING=\n\n#   FST Parameters\n#       Listed below are the defaults, please modify for your samples\n#       Generate site allele frequencies\nDO_SAF=2\n#       Set the minimum number of individuals required for group 1\nMIN_IND1=4\n#       Set the minimum number of individuals required for group 2\nMIN_IND2=4\n#       Overwrite any previously generated results\nOVERRIDE=true\n#       Calculate global Fst values\nGLOBAL=true\n#       Set the sliding window size\nWIN=1000\n#       Set the step size for sliding window analysis\nSTEP=500\n\n#   For advanced users who want to change arguments used by ANGSD (i.e. which SAF method is used)\n#       Expected format is the same as how the flags and arguments are written on the command line:\n#       '-flag1 arg1 -flag2 arg2 ...' \nADVANCED_ARGS=''" > "${outname}"
}

#   Export the function
export -f createFST


#   A function to write a configuration file for NGS_F
function createInbreeding() {
    local common=$1 # Where is Common_Config?
    local project=$2 # What are we calling our configuration file?
    local outdirectory=$3 # Where are we putting our configuration file?
    local outname="${outdirectory}/${project}_Inbreeding_Coefficients_Config" # Create a name for the configuration file
    #   Write out the configuration file for finding inbreeding coefficient using ngsF
    echo -e "#!/bin/bash\n\nset -e\nset -u\nset -o pipefail\n\n#   A simple script to hold variables for the Inbreeding Coefficients\n#   Are you using the Common_Config file?\n#       If so, where is it?\nCOMMON=${common}\n\n##############################################################################################\n#   If we aren't using the Common_Config file, specify these variables\n#   If Common_Config is specified, leave these blank\n#   Define a list of samples\nSAMPLE_LIST=\n\n#   Ancestral and Reference sequences\nANC_SEQ=\nREF_SEQ=\n\n#   Name the project\nPROJECT=\n\n#   Name the project\nPROJECT=\n\n#   Where do we put the outfiles?\n    #   Note, the final outdirectory will be\n    #   \${SCRATCH}/\${PROJECT}/ngsF\nSCRATCH=\n\n#   Region being looked at?\n#       Optional, but ANGSD is expensive to run without specifying regions to look at\nREGIONS=\n
#   Parameter that are specified in Common_Config\n#       Use only uniquely-mapped reads\nUNIQUE_ONLY=0\n#       Set the minimum base quality\nMIN_BASEQUAL=20\n#       Calculate genotype likelihoods\nGT_LIKELIHOOD=1\n#       Set the minimum mapping quality for a base to be used\nMIN_MAPQ=30\n#       Set the number of threads to be used\nN_CORES=32\n#       Determine major and minor alleles\nDO_MAJORMINOR=1\n#       Calculate allele frequencies\nDO_MAF=1\n\n##############################################################################################\n\n#   ngsF Parameters\n#       Listed below are the defaults, please modify for your samples\n#       Set the maximum p-value for polymorphic sites\nSNP_PVAL=1e-6\n#       Overwrite any previously generated results\nOVERRIDE=false\n#       Set the minimum root-mean-square deviation between to assume convergence\nMIN_EPSILON=1e-9\n#       Output genotype likelihood frequency file\nDO_GLF=3\n#       Set a seed value for creating approximate inbreeding coefficients\n#       Use the random number generator built into BASH\nSEED=\$RANDOM\n\n#   For advanced users who want to change arguments used by ANGSD (i.e. which SAF method is used)\n#       Expected format is the same as how the flags and arguments are written on the command line:\n#       '-flag1 arg1 -flag2 arg2 ...' \nADVANCED_ARGS=''" > "${outname}"
}

#   Export the function
export -f createInbreeding


#   A function to write all example configuration files for different projects
function configByProject() {
    local exampleDirectory=$1 # Where is the example data stored?
    local project=$2 # Which example dataset are we making configuration files for?
    local projectDirectory="${exampleDirectory}/${project}" # The example data is stored here
    local configurationDirectory=$3 # Where are we storing our configuration files?
    local ancestralSequence=$4 # Where is our ancestral sequence?
    local referenceSequence=$5 # Where is our reference sequence?
    local scratch="${exampleDirectory}/scratch" # Our scratch directory
    #   Write Common_Config and collect its path
    local commonPath=$(createCommon "${projectDirectory}" "${project}_Samples.txt" "${project}_Inbreeding.indF" "${ancestralSequence}" "${referenceSequence}" "${project}_Example" "${scratch}" "${project}_Regions.txt" "${configurationDirectory}")
    #   Write the configuration file for SFS
    createSFS "${commonPath}" "${project}_Example" "${configurationDirectory}"
    local pestPath="${scratch}/${project}_Example/SFS/${project}_Example_DerivedSFS.graph.me" # Where is our pest file stored?
    #   Write the configuration file for Thetas
    createThetas "${commonPath}" "${pestPath}" "${project}_Example" "${configurationDirectory}"
    #   Write the configuration file for Genotypes
    createGenotypes "${commonPath}" "${project}_Example" "${configurationDirectory}"
    local likelihoodPath="${scratch}/${project}_Example/GenotypeLikelihoods/${project}_Example_snps.beagle.gz" # Where is our likelihood file stored?
    #   Write the configuration file for Admixture
    createAdmixture "${commonPath}" "${likelihoodPath}" "${project}_Example" "${configurationDirectory}"
    #   Write the configuration file for PCA
    createPCA "${commonPath}" "${project}_Example" "${configurationDirectory}"
    #   Write the configuration file for Abbababa
    createAbbababa "${commonPath}" "${project}_Example" "${configurationDirectory}"
    #   Write the configuration file for Ancestral
    createAncestral "${commonPath}" "${project}_Example" "${configurationDirectory}" "${ancestralSequence}"
    #   Write the configuration file for FST
    createFST "${commonPath}" "${project}_Example" "${configurationDirectory}"
    #   Write the configuration file for Inbreeding
    createInbreeding "${commonPath}" "${project}_Example" "${configurationDirectory}"
}

#   Export the function
export -f configByProject

#   A function to manage the entire thing
function writeConfigs() {
    local exampleDirectory=$1 # Where are the example data stored?
    local configurationDirectory="${exampleDirectory}/Configuration_Files" # Where are we storing our example configuration files?
    mkdir -p ${configurationDirectory} # Make our configuration directory
    local ancestralSequence="${exampleDirectory}/Sequences/Tripsacum_TDD39103.fa" # Where is our ancestral sequence?
    local referenceSequence="${exampleDirectory}/Sequences/Zea_mays.AGPv3.30.dna_sm.chromosome.10.fa" # Where is our reference sequence?
    local -a projectNames=(Maize Mexicana Teosinte) # An array with the three different example datasets
    #   Write the example configuration files
    for sample in ${projectNames[@]}; do configByProject "${exampleDirectory}" "${sample}" "${configurationDirectory}" "${ancestralSequence}" "${referenceSequence}"; done
    #   Return the path for the configuration files
    echo "${configurationDirectory}"
}

#   Export the function
export -f writeConfigs
