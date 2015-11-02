# ANGSD-wrapper

ANGSD-wrapper is a utility developed to aid in the analysis of next generation sequencing data. Users can do the following with this suite:
- Calculate a site frequency spectrum
- Calculate a 2D site frequency spectrum with corresponding Fst estimations
- Perform ABBA BABA tests
- Extract a FASTA sequence from BAM files
- Estimate genotype likelihoods
- Estimate thetas and various neutrality statistics
- Calculate per-individual inbreeding coefficient
- Find admixture proportions

Likelihood based approaches are used in ANGSD to calculate summary statistics for next generation sequencing data. The wrapper scripts and documentation aim to make angsd user friendly.

## Installing angsd-wrapper

To install angsd-wrapper, download from GitHub

```shell
git clone https://github.com/mojaveazure/angsd-wrapper.git
```

Go into the angsd-wrapper directory

```shell
cd angsd-wrapper/
```

Run the setup command

```shell
./angsd-wrapper setup please
```

Finish the installation

```shell
source ~/.bash_profile
```

## A note about BAM files

ANGSD requires BAM files as its input, and angsd-wrapper uses a list of BAM files to pass to ANGSD. These BAM files have a few requirements:

- The BAM files must have an '@HD' header line
- The BAM files must be indexed (.bai)

To see whether or not the BAM files have an '@HD' header line, run the following on your list of samples:
```shell
for sample in `cat ~/path/to/sample_list.txt`
do
    echo $sample
    samtools view -H $sample | head -1
done
```

If any samples start with '@SQ' instead of '@HD', ANGSD and angsd-wrapper will fail. This [Gist](https://gist.github.com/mojaveazure/d194c4705642eecf8437) will add an `@HD` header lines to your BAM files.

The index files must be generated after the BAM files. To index the BAM files using SAMTools, run the following on your sample list:

```shell
for sample in `cat ~/path/to/sample_list.txt`
do
    samtools index $sample
done
```

If you have GNU Parallel installed on your system, this process can be sped up:

```shell
cat ~/path/to/sample_list.txt | parallel samtools index {}
```

## Basic usage

To run ANGSD-wrapper, run

```shell
angsd-wrapper <wrapper> <config>
```

Where < wrapper > is one of the available routines that ANGSD-wrapper can run and < config > is the relative path to the corresponding configuration file.

To see a list of available wrappers, run

```shell
angsd-wrapper
```

### Configuration files

For each of the methods `angsd-wrapper` has, there is a config file for it. The configuration files hold variables used by the `wrappers`. This is where you need to modify and save the variables (i.e. specify filepaths of indexed BAM files/CRAM files,  FASTA files, sample lists, etc.) to suit your samples before running angsd-wrapper with a specified method.

The default config files can be found in the `Configuration_Files` directory. You will need to modify them to suit your samples. Please refer to the config files or the [wiki](https://github.com/arundurvasula/angsd-wrapper/wiki) to see what each variable is used for and how they should be specified. If you run `angsd-wrapper` without any arguments, it will spit out a usage message.

## Dependencies
This package requires the following dependencies to work:
 - [ANGSD](https://github.com/angsd/angsd)
 - [ngsPopGen](https://github.com/mfumagalli/ngsPopGen)
 - [ngsF](https://github.com/fgvieira/ngsF)
 - [ngsAdmix](http://www.popgen.dk/software/index.php/NgsAdmix)

These are downloaded and installed automatically when angsd-wrapper is [installed](https://github.com/mojaveazure/angsd-wrapper#installing-angsd-wrapper)

## Supported methods

 - [Site frequency spectrum (SFS)](https://github.com/arundurvasula/angsd-wrapper/wiki/Site-Frequency-Spectrum)
 - [Thetas estimations](https://github.com/arundurvasula/angsd-wrapper/wiki/Thetas)
 - [2D SFS](https://github.com/arundurvasula/angsd-wrapper/wiki/2D-Site-Frequency-Spectrum) and [Fst](https://github.com/arundurvasula/angsd-wrapper/wiki/ngsTools-FST)
 - ~~[Abbababa](https://github.com/arundurvasula/angsd-wrapper/wiki/ABBA-BABA)~~ Work in progress
 - [Ancestral sequence extractions](https://github.com/mojaveazure/angsd-wrapper/blob/master/Wrappers/Ancestral_Sequence.sh)
 - [Genotype likelihood estimations](https://github.com/mojaveazure/angsd-wrapper/blob/master/Wrappers/Genotypes.sh)
 - [ngsF](https://github.com/fgvieira/ngsF)
 - [Principal component analysis](https://github.com/arundurvasula/angsd-wrapper/wiki/Principle-Components-Analysis)
 - Admixture analysis
 - Inbreeding coefficients calculations

## To Do

 - ~~Define requirements for BAM files (@HD/index)~~ DONE!
 - Fix segfaults (In progress with the `ANGSD` team)
 - Define variables in configuration files
 - ~~Fix issue with `git pull`~~ DONE!

