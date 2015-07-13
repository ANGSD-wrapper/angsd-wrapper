#angsd-wrapper
=============

Wrapper scripts and documentation to make angsd easier to use.

Forked from [Arun Durvasula](https://github.com/arundurvasula/angsd-wrapper) in an attempt to simplify the angsd-wrapper package

To read about why I'm doing this, look no further than this [Gist](https://gist.github.com/mojaveazure/ce8c41440805be16c09c).

**NOTE: [ANGSD](https://github.com/angsd/angsd) requires BAM files to have the "@HD" header line. If your BAM files do not have this line, this [Gist](https://gist.github.com/mojaveazure/d194c4705642eecf8437) will add one.**

### To Use

Basic usage is to run the `ANGSD_Wrapper.sh` script, specify a wrapper, and tell it where the config file is. Default config files are found in the [`Configuration_Files`](https://github.com/mojaveazure/angsd-wrapper/tree/master/Configuration_Files) directory; please modify to suit your samples. Running `ANGSD_Wrapper.sh` without any arguments will spit out a usage message.

|**Call** | **Method**|
|_________|___________|
|SFS | Site Frequency Spectrum|


### Supported methods

- [SFS](https://github.com/arundurvasula/angsd-wrapper/wiki/Site-Frequency-Spectrum)
- [Thetas](https://github.com/arundurvasula/angsd-wrapper/wiki/Thetas)
- [2DSFS](https://github.com/arundurvasula/angsd-wrapper/wiki/2D-Site-Frequency-Spectrum)
- ~~[Fst](https://github.com/arundurvasula/angsd-wrapper/wiki/ngsTools-FST)~~
- [ABBA BABA](https://github.com/arundurvasula/angsd-wrapper/wiki/ABBA-BABA)
- [Ancestral Sequence](https://github.com/mojaveazure/angsd-wrapper/blob/master/Wrappers/Ancestral_Sequence.sh)
- [Genotypes](https://github.com/mojaveazure/angsd-wrapper/blob/master/Wrappers/Genotypes.sh)
