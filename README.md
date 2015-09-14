#angsd-wrapper
=============


Angsd-wrapper is a utility developed to aid in the analysis of next generation sequencing data. Users can calculate the following with this tool: 
- site frequency spectrum
- 2D site frequency spectrum and Fst estimations
- ABBA BABBA tests
- extract ancestral sequence from BAM files
- genotype likelihood estimations
- thetas estimations and neutrality test statistics
- inbreeding coefficient calculations
- runs ngsAdmix.  

Likelihood based approaches are used in ANGSD to calculate summary statistics for next generation sequencing data. The wrapper scripts and documentation aim to make angsd user friendly.

This is forked from [Arun Durvasula](https://github.com/arundurvasula/angsd-wrapper) in an attempt to simplify the angsd-wrapper package. To read about why I'm doing this, look no further than this [Gist](https://gist.github.com/mojaveazure/ce8c41440805be16c09c).

## Installing angsd-wrapper
### Here are the steps to install angsd-wrapper using the command-line:
First, you will need to download the angsd-wrapper directory from GitHub

> `git clone https://github.com/mojaveazure/angsd-wrapper.git`

Next you will go into the angsd-wrapper directory and make angsd-wrapper executable

> `cd angsd-wrapper/`

We will now setup angsd-wrapper

> `./angsd-wrapper setup please`

## To Use

### Basic usage
1. Run the `angsd-wrapper` script
2. Specify a wrapper
3. Modify the config file and tell the wrapper where the config file is
4. Run `angsd-wrapper <wrapper> <config_file>` as formatted

The default config files can be found in the `Configuration_Files directory`. You will need to modify them to suit your samples. The angsd-wrapper scripts will create a directory containing calculations from the specified wrapper for the samples that were run. The directory name will be whatever you specify in the configuration files script for the variable `SCRATCH=`.

If you run `angsd-wrapper` without any arguments, it will spit out a usage message.

### Configuration files

For each of the methods `angsd-wrapper` has, there is a config file for it. The configuration files hold variables used by the `wrappers`. This is where you need to modify and save the variables (i.e. specify filepaths of BAM files, FASTA files, sample lists, etc.) to suit your samples before running angsd-wrapper with a specified method.

### Dependencies
This package requires the following dependencies to work: [ANGSD](https://github.com/angsd/angsd), [ngsPopGen](https://github.com/mfumagalli/ngsPopGen), and [ngsF](https://github.com/fgvieira/ngsF) for various methods. In other words, angsd-wrapper depends on using the packages listed above to work. 

### Download and install supported versions of dependencies
Please run `angsd-wrapper setup please` to donwload and install the supported versions of dependencies. 

## Supported methods

- [SFS](https://github.com/arundurvasula/angsd-wrapper/wiki/Site-Frequency-Spectrum)
- [Thetas](https://github.com/arundurvasula/angsd-wrapper/wiki/Thetas)
- [2DSFS](https://github.com/arundurvasula/angsd-wrapper/wiki/2D-Site-Frequency-Spectrum) and [Fst](https://github.com/arundurvasula/angsd-wrapper/wiki/ngsTools-FST)
- [ABBA BABA](https://github.com/arundurvasula/angsd-wrapper/wiki/ABBA-BABA)
- [Ancestral Sequence](https://github.com/mojaveazure/angsd-wrapper/blob/master/Wrappers/Ancestral_Sequence.sh)
- [Genotypes](https://github.com/mojaveazure/angsd-wrapper/blob/master/Wrappers/Genotypes.sh)
- [ngsF](https://github.com/fgvieira/ngsF)

