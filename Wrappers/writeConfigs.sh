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
    echo -e "#!/bin/bash\nset -e \nset -o pipefail \n\n#   A simple script to hold common variables for ANGSD-wrapper\n\n#   Define a list of samples\nGROUP_SAMPLES=${path}/${samples}\n\n#   Define a list of inbreeding coefficients\n#	This should end in '_.indF'\nGROUP_INBREEDING=${path}/${inbreeding}\n\n#   Ancestral sequence\nANC_SEQ=${ancestral}\n\n#   Reference sequence\nREF_SEQ=${reference}\n\n#   Name the project\nPROJECT=${project}\n\n#   Where do we put the outfiles?\n    #   Note, the final outdirectory will be\n    #   \${SCRATCH}/\${PROJECT}/<name_of_program/>\nSCRATCH=${scratch}\n\n#   Define the region being looked at\n#       Optional, but ANGSD is expensive to run without specifying regions to look at\nREGIONS=${path}/${regions}\n\n#   Set common parameters for all methods\nUNIQUE_ONLY=0\nMIN_BASEQUAL=20\nBAQ=1\nMIN_IND=1\nGT_LIKELIHOOD=2\nMIN_MAPQ=30\nN_CORES=32\nDO_MAJORMINOR=1\nDO_GENO=32\nDO_MAF=1\nDO_POST=1\n" > "${outname}"
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
    echo -e "#!/bin/bash\n\nset -e\nset -u\nset -o pipefail\n\n#   A simple script to hold variables for the Site Frequency Spectrum\n#   Are you using the Common_Config file?\n#       If so, where is it?\nCOMMON=${common}\n\n##############################################################################################\n#   If we aren't using the Common_Config file, specify these variables\n#   If Common_Config is specified, leave these blank\n#   Define a list of samples\nSAMPLE_LIST=\n\n#   Define a list of inbreeding coefficients\nSAMPLE_INBREEDING=\n\n#   Ancestral and Reference sequences\nANC_SEQ=\nREF_SEQ=\n\n#   Name the project\nPROJECT=\n\n#   Where do we put the outfiles?\n    #   Note, the final outdirectory will be\n    #   \${SCRATCH}/\${PROJECT}/SFS\nSCRATCH=\n\n#   Define the region being looked at\n#       Optional, but ANGSD is expensive to run without specifying regions to look at\nREGIONS=\n\n#   Parameters that are specified in Common_Config\nUNIQUE_ONLY=0\nMIN_BASEQUAL=20\nBAQ=1\nMIN_IND=1\nGT_LIKELIHOOD=2\nMIN_MAPQ=30\nN_CORES=32\nDO_MAJORMINOR=1\nDO_GENO=32\nDO_MAF=1\nDO_POST=1\n\n##############################################################################################\n\n#   Site Frequency Spectrum Parameters\n#       Listed below are the defaults, please modify for your samples\nDO_SAF=2\nOVERRIDE=true\n" > "${outname}"
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
    echo -e "#!/bin/bash\n\nset -e\nset -u\nset -o pipefail\n\n#   A simple script to hold variables for the Estimation of Thetas\n#   Are you using the Common_Config file?\n#       If so, where is it?\nCOMMON=${common}\n\n##############################################################################################\n#   If we aren't using the Common_Config file, specify these variables\n#   If Common_Config is specified, leave these blank\n#   Define a list of samples\nSAMPLE_LIST=\n\n#   Define a list of inbreeding coefficients\nSAMPLE_INBREEDING=\n\n#   Ancestral and Reference sequences\nANC_SEQ=\nREF_SEQ=\n\n#   Name the project\nPROJECT=\n\n#   Where do we put the outfiles?\n    #   Note, the final outdirectory will be\n    #   \${SCRATCH}/\${PROJECT}/Thetas\nSCRATCH=\n\n#   Define the region being looked at\n#       Optional, but ANGSD is expensive to run without specifying regions to look at\nREGIONS=\n\n#   Parameters that are specified in Common_Config\nUNIQUE_ONLY=0\nMIN_BASEQUAL=20\nBAQ=1\nMIN_IND=1\nGT_LIKELIHOOD=2\nMIN_MAPQ=30\nN_CORES=32\nDO_MAJORMINOR=1\nDO_MAF=1\n\n##############################################################################################\n\n#   Pest file\n#       This is the output from the site frequency spectrum wrapper\n#       This should end in '_DerivedSFS'\nPEST=${pest}\n\n#   Thetas Parameters\n#       Listed below are the defaults, please modify for your samples\nDO_SAF=2\nDO_THETAS=1\nOVERRIDE=true\nSLIDING_WINDOW=false\nWIN=50000\nSTEP=10000\n" > "${outname}"
}

#   Export the function
export -f createThetas

#   A function to write a configuration file for Admixture Analysis
function createAdmixture() {
    local common=$1 # Where is Common_Config?
    local likelihood=$2 # Where is the likelihood file
    local project=$3 # What are we calling our configuration file?
    local outdirectory=$4 # Where are we putting our configuration file?
    local outname="${outdirectory}/${project}_Admixture_Config" # Create a name for the configuration file
    #   Write out the configuration file for Admixture
    echo -e "#!/bin/bash\n\nset -e\nset -u\nset -o pipefail\n\n#   A simple script to hold the varialbes for the NGS Admixture\n#   Are you using the Common_Config file?\n#       If so, where is it?\nCOMMON=${common}\n\n##############################################################################################\n#   If we aren't using the Common_Config file, specify these variables\n#   If Common_Config is specified, leave these blank\n#   Name the project\nPROJECT=\n\n#   Where do we put the outfiles?\n    #   Note, the final outdirectory will be\n    #   \${SCRATCH}/\${PROJECT}/Admixture\nSCRATCH=\n\n#   Parameters that are specified in Common_Config\nN_CORES=32\n\n##############################################################################################\n\n#   The Likelihood file\n#       This is the .beagle.gz file from the Site Frequency Spectrum\nLIKELIHOOD=${likelihood}\n\n#   ngsAdmix Parameters\n#       Listed below are the defaults, please modify for your samples\nK=5\nMIN_MAF=0.05\n" > "${outname}"
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
    echo -e "#!/bin/bash\n\nset -e\nset -u\nset -o pipefail\n\n#   A simple script to hold variables for the Principal Component Analysis\n#   Are you using the Common_Config file?\n#       If so, where is it?\nCOMMON=${common}\n\n##############################################################################################\n#   If we aren't using the Common_Config file, specify these variables\n#   If Common_Config is specified, leave these blank\n#   Define a list of samples\nSAMPLE_LIST=\n\n#   Name the project\nPROJECT=\n\n#   Where do we put the outfiles?\n    #   Note, the final outdirectory will be\n    #   \${SCRATCH}/\${PROJECT}/PCA\nSCRATCH=\n\n#   Region being looked at?\nREGIONS=\n##############################################################################################\n\n#   Principal Component Analysis Parameters\n#       Listed below are the defaults, please modify for your samples\nDO_MAF=2\nDO_MAJORMINOR=1\nDO_GENO=32\nDO_POST=1\nN_CORES=8\nNORM=0\nCALL=0\nGT_LIKELIHOOD=1\nN_SITES=100000\n" > "${outname}"
}

#   Export the function
export -f createPCA

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
    local likelihoodPath="${scratch}/${project}_Example/SFS/${project}_Example.beagle.gz" # Where is our likelihood file stored?
    #   Write the configuration file for Thetas
    createThetas "${commonPath}" "${pestPath}" "${project}_Example" "${configurationDirectory}"
    #   Write the configuration file for Admixture
    createAdmixture "${commonPath}" "${likelihoodPath}" "${project}_Example" "${configurationDirectory}"
    #   Write the configuration file for PCA
    createPCA "${commonPath}" "${project}_Example" "${configurationDirectory}"
}

#   Export the function
export -f configByProject

#   A function to manage the entire thing
function writeConfigs() {
    local exampleDirectory=$1 # Where are the example data stored?
    local configurationDirectory="${exampleDirectory}/Configuration_Files" # Where are we storing our example configuration files?
    mkdir -p ${configurationDirectory} # Make our configuration directory
    local ancestralSequence="${exampleDirectory}/Sequences/" # Where is our ancestral sequence?
    local referenceSequence="${exampleDirectory}/Sequences/" # Where is our reference sequence?
    local -a projectNames=(Maize Mexicana Teosinte) # An array with the three different example datasets
    #   Write the example configuration files
    for sample in ${projectNames[@]}; do configByProject "${exampleDirectory}" "${sample}" "${configurationDirectory}" "${ancestralSequence}" "${referenceSequence}"; done
    #   Return the path for the configuration files
    echo "${configurationDirectory}"
}

#   Export the function
export -f writeConfigs
