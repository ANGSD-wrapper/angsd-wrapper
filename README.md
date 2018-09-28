<<<<<<< HEAD
# angsd-wrapper
Utilities for analyzing next generation sequencing data.
=======
#angsd-wrapper
=============

__Update:__ The current version of ANGSD-wrapper is available at: https://github.com/mojaveazure/angsd-wrapper

## This version is now deprecated.

---
Wrapper scripts and documentation to make angsd easier to use.

### Installation
In order to install angsd-wrapper fully, you need to follow a few steps:

First, you need to make the ANGSD executable and the associated programs by copying this in your command line:

```
git clone https://github.com/arundurvasula/angsd-wrapper.git
cd angsd-wrapper/angsd
make
cd ../ngsPopGen
make
cd ../ngsF
make
```

Then you need to install the requisite R packages for using the R graphing application:

```
R
install.packages("shiny")
install.packages("ape")
install.packages("Hmisc")
install.packages("lattice")
source("http://bioconductor.org/biocLite.R")
biocLite("genomeIntervals")
install.packages("data.table")
```

### Documentation available on the [wiki](https://github.com/arundurvasula/angsd-wrapper/wiki).


### Contributing
Contributions are very welcome! You can help by submitting [issues](https://github.com/arundurvasula/angsd-wrapper/issues) on Github or by forking the repository and making changes yourself.

>>>>>>> angsd-wrapper/master
