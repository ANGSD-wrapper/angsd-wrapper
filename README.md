#angsd-wrapper
=============

Wrapper scripts and documentation to make angsd more usable

Table of Contents:

0. [Usage](#usage)
1. [Overview](#overview)
2. [Configuration file guide](#config)
3. [ANGSD guide](#angsd)

##0. <a name="usage"></a>Usage

    $ bash scripts/ANGSD_SFS.sh scripts/SFS_example.conf

Scripts should be invoked with `bash` and executed from the toplevel directory so the paths are correct. 

##1. <a name="overview"></a>Overview
The scripts contained in this directory have been created to make the ANGSD workflow easier to use. It contains scripts from the [Ross-Ibarra](https://github.com/rossibarra/angsbigd) and [Morrell](https://github.com/MorrellLAB/angsbigd) labs. 

[ANGSD](https://github.com/ANGSD/angsd) is included in this repository as a [git subtree](https://hpc.uni.lu/blog/2014/understanding-git-subtree/). This allows the repository to keep up with the latest changes to ANGSD, and ensures a consistent version of ANGSD. If you need to pull the latest changes to angsd, use the following git command:

    git subtree pull --prefix=angsd --squash angsd master

Note, you will have to run `make` in the angsd directory to make the executables. 

##2. <a name="config"></a>Configuration file guide
###Config syntax 
Configuration files are interpreted by the script as pure bash. As such, they should follow bash syntax to avoid any errors. A full list of variables follows:
###SFS
####Required variables

- `UNIX_USER`: this variable fills in absolute paths for the rest of the config file and script. It should match the name of the user's home directory.
- `ANC_SEQ`: the path to the ancestral sequence file
- `REF_SEQ`: the path to the reference sequence file
- `TAXON`: the name of the data being analyzed. The script will look for files in the data directory with this name. These files include: `${TAXON}_samples.txt` and `${TAXON}_F.txt`. If these files are not present, the script will not work correctly. `${TAXON}_samples.txt` contains a list of paths to BAM files. Check the data folder for an example. `${TAXON}_F.txt` contains inbreeding coefficients for each of these samples. Check the data folder for an example.

####Optional variables

These variables can be included and will override the default values. 
- `DO_SAF` create SFS (default=2)
- `UNIQUE_ONLY` uniquely mapped reads (default=1)
- `MIN_BASEQUAL` minimum base quality (default=20)
- `BAQ` adjust qscores around indels (as SAMtools) (default=1)
- `MIN_IND` minimum number of individuals needed to use site (default=1)
- `GT_LIKELIHOOD` estimate genotype likelihoods (default=2)
- `MIN_MAPQ` minimum base mapping quality to use (default=30)
- `N_CORES` number of cores to use (default=32)
- `DO_MAJORMINOR` estimate major/minor alleles (default=1)
- `DO_MAF` calculate per site frequencies (default=1)
- `REGIONS` chromosome and region to use (default="1:")

##3. <a name="angsd"></a>ANGSD guide
Coming soon!