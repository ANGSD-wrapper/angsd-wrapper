#!/bin/bash

# set -e
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
    echo -e "#!/bin/bash\n \
set -e \n \
set -o pipefail \n \
\n \
#   A simple script to hold common variables for ANGSD-wrapper\n \
\n \
#   Define a list of samples\n \
SAMPLE_LIST=${path}/${samples}\n \
\n \
#   Define a list of inbreeding coefficients\n \
#	This should end in '_.indF'\n \
SAMPLE_INBREEDING=${path}/${inbreeding}\n \
\n \
#   Ancestral sequence\n \
ANC_SEQ=${ancestral}\n \
\n \
#   Reference sequence\n \
REF_SEQ=${reference}\n \
\n \
#   Name the project\n \
PROJECT=${project}\n \
\n \
#   Where do we put the outfiles?\n \
    #   Note, the final outdirectory will be\n \
    #   \${SCRATCH}/\${PROJECT}/<name_of_program/>\n \
SCRATCH=${scratch}\n \
\n \
#   Define the region being looked at\n \
#       Optional, but ANGSD is expensive to run without specifying regions to look at\n \
REGIONS=${path}/${regions}\n \
\n \
#   Set common parameters for all methods\n \
UNIQUE_ONLY=0\n \
MIN_BASEQUAL=20\n \
BAQ=1\n \
MIN_IND=1\n \
GT_LIKELIHOOD=2\n \
MIN_MAPQ=30\n \
N_CORES=32\n \
DO_MAJORMINOR=1\n \
DO_GENO=32\n \
DO_MAF=1\n \
DO_POST=1\n \
" > "${outname}"
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
    echo -e "#!/bin/bash\n \
\n \
set -e\n \
set -u\n \
set -o pipefail\n \
\n \
#   A simple script to hold variables for the Site Frequency Spectrum\n \
#   Are you using the Common_Config file?\n \
#       If so, where is it?\n \
COMMON=${common}\n \
\n \
##############################################################################################\n \
#   If we aren't using the Common_Config file, specify these variables\n \
#   If Common_Config is specified, leave these blank\n \
#   Define a list of samples\n \
SAMPLE_LIST=\n \
\n \
#   Define a list of inbreeding coefficients\n \
SAMPLE_INBREEDING=\n \
\n \
#   Ancestral and Reference sequences\n \
ANC_SEQ=\n \
REF_SEQ=\n \
\n \
#   Name the project\n \
PROJECT=\n \
\n \
#   Where do we put the outfiles?\n \
    #   Note, the final outdirectory will be\n \
    #   \${SCRATCH}/\${PROJECT}/SFS\n \
SCRATCH=\n \
\n \
#   Define the region being looked at\n \
#       Optional, but ANGSD is expensive to run without specifying regions to look at\n \
REGIONS=\n \
\n \
#   Parameters that are specified in Common_Config\n \
UNIQUE_ONLY=0\n \
MIN_BASEQUAL=20\n \
BAQ=1\n \
MIN_IND=1\n \
GT_LIKELIHOOD=2\n \
MIN_MAPQ=30\n \
N_CORES=32\n \
DO_MAJORMINOR=1\n \
DO_GENO=32\n \
DO_MAF=1\n \
DO_POST=1\n \
\n \
##############################################################################################\n \
\n \
#   Site Frequency Spectrum Parameters\n \
#       Listed below are the defaults, please modify for your samples\n \
#       Generate site allele frequencies\n \
DO_SAF=2\n \
#       Overwrite any previously generated results\n \
OVERRIDE=true\n \
\n \
#   For advanced users who want to change arguments used by ANGSD (i.e. which SAF method is used)\n \
#       Expected format is the same as how the flags and arguments are written on the command line:\n \
#       '-flag1 arg1 -flag2 arg2 ...' \n \
ADVANCED_ARGS=''" > "${outname}"
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
    echo -e "#!/bin/bash\n \
\n \
set -e\n \
set -u\n \
set -o pipefail\n \
\n \
#   A simple script to hold variables for the Estimation of Thetas\n \
#   Are you using the Common_Config file?\n \
#       If so, where is it?\n \
COMMON=${common}\n \
\n \
##############################################################################################\n \
#   If we aren't using the Common_Config file, specify these variables\n \
#   If Common_Config is specified, leave these blank\n \
#   Define a list of samples\n \
SAMPLE_LIST=\n \
\n \
#   Define a list of inbreeding coefficients\n \
SAMPLE_INBREEDING=\n \
\n \
#   Ancestral and Reference sequences\n \
ANC_SEQ=\n \
REF_SEQ=\n \
\n \
#   Name the project\n \
PROJECT=\n \
\n \
#   Where do we put the outfiles?\n \
    #   Note, the final outdirectory will be\n \
    #   \${SCRATCH}/\${PROJECT}/Thetas\n \
SCRATCH=\n \
\n \
#   Define the region being looked at\n \
#       Optional, but ANGSD is expensive to run without specifying regions to look at\n \
REGIONS=\n \
\n \
#   Parameters that are specified in Common_Config\n \
UNIQUE_ONLY=0\n \
MIN_BASEQUAL=20\n \
BAQ=1\n \
MIN_IND=1\n \
GT_LIKELIHOOD=2\n \
MIN_MAPQ=30\n \
N_CORES=32\n \
DO_MAJORMINOR=1\n \
DO_MAF=1\n \
\n \
##############################################################################################\n \
\n \
#   Pest file\n \
#       This is the output from the site frequency spectrum wrapper\n \
#       This should end in '_DerivedSFS'\n \
PEST=${pest}\n \
\n \
#   Thetas Parameters\n \
#       Listed below are the defaults, please modify for your samples\n \
DO_SAF=2\n \
DO_THETAS=1\n \
OVERRIDE=true\n \
SLIDING_WINDOW=false\n \
WIN=50000\n \
STEP=10000\n \
\n \
#   For advanced users who want to change arguments used by ANGSD (i.e. which SAF method is used)\n \
#       Expected format is the same as how the flags and arguments are written on the command line:\n \
#       '-flag1 arg1 -flag2 arg2 ...' \n \
ADVANCED_ARGS=''" > "${outname}"
}

#   Export the function
export -f createThetas

#   A function to write a configuration file for Genotype Likelihoods
function createGenotypes() {
    local common="$1" # Where is Common_Config?
    local project="$2" # What are we calling our configuration file?
    local outdirectory="$3" # Where are we putting our configuration file?
    local outname="${outdirectory}/${project}_Genotypes_Config" # Create a name for the configuration file
    echo -e "#!/bin/bash\n \
\n \
set -e\n \
set -u\n \
set -o pipefail\n \
\n \
#   A simple script to hold variables for the Site Frequency Spectrum\n \
#   Are you using the Common_Config file?\n \
#       If so, where is it?\n \
COMMON=${common}\n \
\n \
##############################################################################################\n \
#   If we aren't using the Common_Config file, specify these variables\n \
#   If Common_Config is specified, leave these blank\n \
#   Define a list of samples\n \
SAMPLE_LIST=\n \
\n \
#   Define a list of inbreeding coefficients\n \
SAMPLE_INBREEDING=\n \
\n \
#   Name the project\n \
PROJECT=\n \
\n \
#   Where do we put the outfiles?\n \
    #   Note, the final outdirectory will be\n \
    #   \${SCRATCH}/\${PROJECT}/GenotypeLikelihoods\n \
SCRATCH=\n \
\n \
#   Define the region being looked at\n \
#       Optional, but ANGSD is expensive to run without specifying regions to look at\n \
REGIONS=\n \
\n \
#   Set common parameters for all methods\n \
#       Use only uniquely-mapped reads\n \
UNIQUE_ONLY=0\n \
#       Set the minimum base quality\n \
MIN_BASEQUAL=20\n \
#       Calculate base alignment quality\n \
BAQ=1\n \
#       Set the minimum number of individuals required\n \
MIN_IND=1\n \
#       Calculate genotype likelihoods\n \
GT_LIKELIHOOD=2\n \
#       Set the minimum mapping quality for a base to be used\n \
MIN_MAPQ=30\n \
#       Set the number of threads to be used\n \
N_CORES=32\n \
#       Determine major and minor alleles\n \
DO_MAJORMINOR=1\n \
#       Call genotypes from genotype likelihoods\n \
DO_GENO=32\n \
#       Calculate allele frequencies\n \
DO_MAF=1\n \
#       Calculate the posterior probability\n \
DO_POST=1\n \
\n \
##############################################################################################\n \
\n \
#   Genotypes Parameters\n \
#       Listed below are the defaults, please modify for your samples\n \
#       Set the minimum posterior value for calling genotypes\n \
POST_CUTOFF=0.95\n \
#       Set the maximum p-value for polymorphic sites\n \
SNP_PVAL=1e-6\n \
#       Output genotype likelihood frequency file\n \
DO_GLF=2\n \
\n \
#   For advanced users who want to change arguments used by ANGSD (i.e. which SAF method is used)\n \
#       Expected format is the same as how the flags and arguments are written on the command line:\n \
#       '-flag1 arg1 -flag2 arg2 ...' \n \
ADVANCED_ARGS=''" > "${outname}"
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
    echo -e "#!/bin/bash\n \
\n \
set -e\n \
set -u\n \
set -o pipefail\n \
\n \
#   A simple script to hold the varialbes for the NGS Admixture\n \
#   Are you using the Common_Config file?\n \
#       If so, where is it?\n \
COMMON=${common}\n \
\n \
##############################################################################################\n \
#   If we aren't using the Common_Config file, specify these variables\n \
#   If Common_Config is specified, leave these blank\n \
#   Name the project\n \
PROJECT=\n \
\n \
#   Where do we put the outfiles?\n \
    #   Note, the final outdirectory will be\n \
    #   \${SCRATCH}/\${PROJECT}/Admixture\n \
SCRATCH=\n \
\n \
#   Parameters that are specified in Common_Config\n \
N_CORES=32\n \
\n \
##############################################################################################\n \
\n \
#   The Likelihood file\n \
#       This is the .beagle.gz file from the Site Frequency Spectrum\n \
LIKELIHOOD=${likelihood}\n \
\n \
#   ngsAdmix Parameters\n \
#       Listed below are the defaults, please modify for your samples\n \
K=5\n \
MIN_MAF=0.05\n \
TOLERANCE=0.01\n \
\n \
\n \
#   For advanced users who want to change arguments used by ANGSD (i.e. which SAF method is used)\n \
#       Expected format is the same as how the flags and arguments are written on the command line:\n \
#       '-flag1 arg1 -flag2 arg2 ...' \n \
ADVANCED_ARGS=''" > "${outname}"
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
    echo -e "#!/bin/bash\n \
\n \
set -e\n \
set -u\n \
set -o pipefail\n \
\n \
#   A simple script to hold variables for the Principal Component Analysis\n \
#   Are you using the Common_Config file?\n \
#       If so, where is it?\n \
COMMON=${common}\n \
\n \
##############################################################################################\n \
#   If we aren't using the Common_Config file, specify these variables\n \
#   If Common_Config is specified, leave these blank\n \
#   Define a list of samples\n \
SAMPLE_LIST=\n \
\n \
#   Name the project\n \
PROJECT=\n \
\n \
#   Where do we put the outfiles?\n \
    #   Note, the final outdirectory will be\n \
    #   \${SCRATCH}/\${PROJECT}/PCA\n \
SCRATCH=\n \
\n \
#   Region being looked at?\n \
REGIONS=\n \
##############################################################################################\n \
\n \
#   Principal Component Analysis Parameters\n \
#       Listed below are the defaults, please modify for your samples\n \
DO_MAF=2\n \
DO_MAJORMINOR=1\n \
DO_GENO=32\n \
DO_POST=1\n \
N_CORES=8\n \
NORM=0\n \
CALL=0\n \
GT_LIKELIHOOD=2\n \
N_SITES=100000\n \
\n \
#   For advanced users who want to change arguments used by ANGSD (i.e. which SAF method is used)\n \
#       Expected format is the same as how the flags and arguments are written on the command line:\n \
#       '-flag1 arg1 -flag2 arg2 ...' \n \
ADVANCED_ARGS=''" > "${outname}"
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
   echo -e "#!/bin/bash\n \
\n \
set -e\n \
set -u\n \
set -o pipefail\n \
\n \
#   A simple script to hold variables for ABBA BABA\n \
#   Are you using the Common_Config file?\n \
#       If so, where is it?\n \
COMMON=${common}\n \
\n \
##############################################################################################\n \
#   If we aren't using the Common_Config file, specify these variables\n \
#   If Common_Config is specified, leave these blank\n \
#   Define a list of samples\n \
SAMPLE_LIST=\n \
\n \
#   Name the project\n \
PROJECT=\n \
\n \
#   Where do we put the outfiles?\n \
#      Note, the final outdirectory will be\n \
#      \${SCRATCH}/\${PROJECT}/ABBABABA\n \
SCRATCH=\n \
\n \
#   Region being looked at?\n \
#       Optional, but ANGSD is expensive to run without specifying regions to look at\n \
REGIONS=\n \
\n \
#   Parameters that are specified in Common_Config\n \
#       Use only uniquely-mapped reads\n \
UNIQUE_ONLY=0\n \
#       Set the minimum base quality\n \
MIN_BASEQUAL=20\n \
#       Set the minimum number of individuals required\n \
MIN_IND=1\n \
#       Set the minimum mapping quality for a base to be used\n \
MIN_MAPQ=30\n \
#       Set the number of threads to be used\n \
N_CORES=32\\n \
\n \
##############################################################################################\n \
\n \
#   Fasta file to be used as an outgroup\n \
#       Defaults to ANC_SEQ if left undefined\n \
OUTGROUP=\n \
\n \
#   ABBA BABA Parameters\n \
#       Listed below are the defaults, please modify for your samples\n \
#       Count allele frequencies\n \
DO_COUNTS=1\n \
#       Perform Abbababa analysis\n \
DO_ABBABABA=1\n \
#       Remove transitions from sample\n \
REMOVE_TRANS=0\n \
#       Set the size for each block\n \
BLOCKSIZE=1000\n \
\n \
#   Use the last individual in the bam file as outgroup instead of a fasta file\n \
USE_LAST=\n \
\n \
#   Use only sites where the reads for the outgroup has the same base for all reads. Only works with USE_LAST=1\n \
ENHANCE=\n \
\n \
#   For advanced users who want to change arguments used by ANGSD (i.e. which SAF method is used)\n \
#       Expected format is the same as how the flags and arguments are written on the command line:\n \
#       '-flag1 arg1 -flag2 arg2 ...' \n \
ADVANCED_ARGS=''" > "${outname}"

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
    echo -e "#!/bin/bash\n \
\n \
set -e\n \
set -u\n \
set -o pipefail\n \
\n \
#   A simple script to hold variables for generating an ancestral fasta file\n \
\n \
#   This script does NOT utilize the Common_Config file\n \
\n \
#   Where is the ancestral BAM file?\n \
ANC_BAM=${anc_bam}\n \
\n \
#   Full path to the output file\n \
OUT=$(dirname "${anc_bam}")/$(basename "${anc_bam}" .bam).fa\n \
\n \
#   Ancestral Sequence Parameters\n \
#       Listed below are the defaults, please modify for your samples\n \
#       Extract FASTA sequence from BAM file\n \
DO_FASTA=1\n \
#       Count allele frequencies\n \
#       If DO_FASTA is 2, DO_COUNTS must be 1\n \
#       Otherwise, DO_COUNTS can be any other legal value\n \
DO_COUNTS=0\n \
\n \
#   For advanced users who want to change arguments used by ANGSD (i.e. which SAF method is used)\n \
#       Expected format is the same as how the flags and arguments are written on the command line:\n \
#       '-flag1 arg1 -flag2 arg2 ...' \n \
ADVANCED_ARGS=''" > "${outname}"
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
    echo -e "#!/bin/bash\n \
\n \
set -e\n \
set -u\n \
set -o pipefail\n \
\n \
#   A simple script to hold variables for the FST\n \
#   Are you using the Common_Config file?\n \
#       If so, where is it?\n \
COMMON=${common}\n \
\n \
##############################################################################################\n \
#   If we aren't using the Common_Config file, specify these variables\n \
#   If Common_Config is specified, leave these blank\n \
#   Ancestral and Reference sequences\n \
ANC_SEQ=${ancestral}\n \
REF_SEQ=${reference}\n \
\n \
#   Name the project\n \
PROJECT=${project}\n \
\n \
#   Where do we put the outfiles?\n \
    #   Note, the final outdirectory will be\n \
    #   \${SCRATCH}/\${PROJECT}/Fst\n \
SCRATCH=\n \
\n \
#   Region being looked at?\n \
REGIONS=\n \
\n \
#   Parameters that are specified in Common_Config\n \
#       Use only uniquely-mapped reads\n \
UNIQUE_ONLY=0\n \
#       Set the minimum base quality\n \
MIN_BASEQUAL=20\n \
#       Calculate base alignment quality\n \
BAQ=1\n \
#       Calculate genotype likelihoods\n \
GT_LIKELIHOOD=2\n \
#       Set the minimum mapping quality for a base to be used\n \
MIN_MAPQ=30\n \
#       Set the number of threads to be used\n \
N_CORES=32\n \
#       Determine major and minor alleles\n \
DO_MAJORMINOR=1\n \
#       Call genotypes from genotype likelihoods\n \
DO_GENO=32\n \
#       Calculate allele frequencies\n \
DO_MAF=1\n \
#       Calculate the posterior probability\n \
DO_POST=1\n \
\n \
##############################################################################################\n \
\n \
#   What is group 1?\n \
GROUP_1=\n \
\n \
#   Sample list for group 1\n \
G1_SAMPLE_LIST=\n \
\n \
#   Inbreeding coefficients for group 1\n \
G1_INBREEDING=\n \
\n \
#   What is group 2?\n \
GROUP_2=\n \
\n \
#   Sample list for group 2\n \
G2_SAMPLE_LIST=\n \
\n \
#   Inbreeding coefficients for group 2\n \
G2_INBREEDING=\n \
\n \
#   FST Parameters\n \
#       Listed below are the defaults, please modify for your samples\n \
#       Generate site allele frequencies\n \
DO_SAF=2\n \
#       Set the minimum number of individuals required for group 1\n \
MIN_IND1=4\n \
#       Set the minimum number of individuals required for group 2\n \
MIN_IND2=4\n \
#       Overwrite any previously generated results\n \
OVERRIDE=true\n \
#       Calculate global Fst values\n \
GLOBAL=true\n \
#       Set the sliding window size\n \
WIN=1000\n \
#       Set the step size for sliding window analysis\n \
STEP=500\n \
\n \
#   For advanced users who want to change arguments used by ANGSD (i.e. which SAF method is used)\n \
#       Expected format is the same as how the flags and arguments are written on the command line:\n \
#       '-flag1 arg1 -flag2 arg2 ...' \n \
ADVANCED_ARGS=''" > "${outname}"
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
    echo -e "#!/bin/bash\n \
\n \
set -e\n \
set -u\n \
set -o pipefail\n \
\n \
#   A simple script to hold variables for the Inbreeding Coefficients\n \
#   Are you using the Common_Config file?\n \
#       If so, where is it?\n \
COMMON=${common}\n \
\n \
##############################################################################################\n \
#   If we aren't using the Common_Config file, specify these variables\n \
#   If Common_Config is specified, leave these blank\n \
#   Define a list of samples\n \
SAMPLE_LIST=\n \
\n \
#   Ancestral and Reference sequences\n \
ANC_SEQ=\n \
REF_SEQ=\n \
\n \
#   Name the project\n \
PROJECT=\n \
\n \
#   Name the project\n \
PROJECT=\n \
\n \
#   Where do we put the outfiles?\n \
    #   Note, the final outdirectory will be\n \
    #   \${SCRATCH}/\${PROJECT}/ngsF\n \
SCRATCH=\n \
\n \
#   Region being looked at?\n \
#       Optional, but ANGSD is expensive to run without specifying regions to look at\n \
REGIONS=\n \

#   Parameter that are specified in Common_Config\n \
#       Use only uniquely-mapped reads\n \
UNIQUE_ONLY=0\n \
#       Set the minimum base quality\n \
MIN_BASEQUAL=20\n \
#       Calculate genotype likelihoods\n \
GT_LIKELIHOOD=1\n \
#       Set the minimum mapping quality for a base to be used\n \
MIN_MAPQ=30\n \
#       Set the number of threads to be used\n \
N_CORES=32\n \
#       Determine major and minor alleles\n \
DO_MAJORMINOR=1\n \
#       Calculate allele frequencies\n \
DO_MAF=1\n \
\n \
##############################################################################################\n \
\n \
#   ngsF Parameters\n \
#       Listed below are the defaults, please modify for your samples\n \
#       Set the maximum p-value for polymorphic sites\n \
SNP_PVAL=1e-6\n \
#       Overwrite any previously generated results\n \
OVERRIDE=false\n \
#       Set the minimum root-mean-square deviation between to assume convergence\n \
MIN_EPSILON=1e-9\n \
#       Output genotype likelihood frequency file\n \
DO_GLF=3\n \
#       Set a seed value for creating approximate inbreeding coefficients\n \
#       Use the random number generator built into BASH\n \
SEED=\$RANDOM\n \
\n \
#   For advanced users who want to change arguments used by ANGSD (i.e. which SAF method is used)\n \
#       Expected format is the same as how the flags and arguments are written on the command line:\n \
#       '-flag1 arg1 -flag2 arg2 ...' \n \
ADVANCED_ARGS=''" > "${outname}"
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
    for sample in ${projectNames[@]};
    do
        configByProject "${exampleDirectory}" "${sample}" "${configurationDirectory}" "${ancestralSequence}" "${referenceSequence}"
    done
    #   Return the path for the configuration files
    echo "${configurationDirectory}"
}

#   Export the function
export -f writeConfigs
