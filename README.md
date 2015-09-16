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

To install angsd-wrapper, download from GitHub

> `git clone https://github.com/mojaveazure/angsd-wrapper.git`

Go into the angsd-wrapper directory

> `cd angsd-wrapper/`

Run the setup command

> `./angsd-wrapper setup please`

Finish the installation

> `source ~/.bash_profile`

## To Use

### A few requirements

- indexed BAM (.bai) or CRAM (.cram) file formats
- @HD header lines

`angsd-wrapper` does not work for BAM files and has specific header line requirements, so please index the BAM files being used and add the @HD header lines before running `angsd-wrapper`.

### Basic usage

To run angsd-wrapper, run

> `angsd-wrapper <wrapper> <config>`

Where < wrapper > is one of the available routines that angsd-wrapper can run and < config > is the relative path to the corresponding configuration file.

To see a list of available wrappers, run

> `angsd-wrapper`

### Configuration files

For each of the methods `angsd-wrapper` has, there is a config file for it. The configuration files hold variables used by the `wrappers`. This is where you need to modify and save the variables (i.e. specify filepaths of indexed BAM files/CRAM files,  FASTA files, sample lists, etc.) to suit your samples before running angsd-wrapper with a specified method.

The default config files can be found in the `Configuration_Files` directory. You will need to modify them to suit your samples. Please refer to the config files or the [wiki](https://github.com/arundurvasula/angsd-wrapper/wiki) to see what each variable is used for and how they should be specified. If you run `angsd-wrapper` without any arguments, it will spit out a usage message.

### Dependencies
This package requires the following dependencies to work: [ANGSD](https://github.com/angsd/angsd), [ngsPopGen](https://github.com/mfumagalli/ngsPopGen), and [ngsF](https://github.com/fgvieira/ngsF) for various methods. In other words, angsd-wrapper depends on using the packages listed above to work.

### Download and install supported versions of dependencies
Please run `./angsd-wrapper setup please` to donwload and install the supported versions of dependencies.

## Supported methods

- [SFS](https://github.com/arundurvasula/angsd-wrapper/wiki/Site-Frequency-Spectrum)
- [Thetas](https://github.com/arundurvasula/angsd-wrapper/wiki/Thetas)
- [2DSFS](https://github.com/arundurvasula/angsd-wrapper/wiki/2D-Site-Frequency-Spectrum) and [Fst](https://github.com/arundurvasula/angsd-wrapper/wiki/ngsTools-FST)
- [ABBA BABA](https://github.com/arundurvasula/angsd-wrapper/wiki/ABBA-BABA)
- [Ancestral Sequence](https://github.com/mojaveazure/angsd-wrapper/blob/master/Wrappers/Ancestral_Sequence.sh)
- [Genotypes](https://github.com/mojaveazure/angsd-wrapper/blob/master/Wrappers/Genotypes.sh)
- [ngsF](https://github.com/fgvieira/ngsF)

## To Do

 - Define requirements for BAM files (@HD/index)
 - Fix segfaults
 - Define variables in configuration files
 - Fix issue with `git pull`

