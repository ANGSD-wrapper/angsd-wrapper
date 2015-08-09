#angsd-wrapper
=============

Wrapper scripts and documentation to make angsd easier to use.

### Installation
In order to install angsd-wrapper fully, you need to follow a few steps:

First, you need to make the ANGSD executable and the associated programs:

1. `$ git clone https://github.com/arundurvasula/angsd-wrapper.git`
2. `$ cd angsd-wrapper/angsd`
3. `$ make`
4. `$ cd ../ngsPopGen`
5. `$ make`
6. `$ cd ../ngsF`
7. `$ make`

Then you need to install the requisite R packages for using the R graphing application:

1. `$ R`
2. `> install.packages("shiny")`
3. `> install.packages("ape")`
4. `> install.packages("Hmisc")`
5. `> install.packages("lattice")`
6. `> source("http://bioconductor.org/biocLite.R")`
7. `> biocLite("genomeIntervals")`
8. `> install.packages("data.table")`

### Documentation available on the [wiki](https://github.com/arundurvasula/angsd-wrapper/wiki).


### Contributing
Contributions are very welcome! You can help by submitting [issues](https://github.com/arundurvasula/angsd-wrapper/issues) on Github or by forking the repository and making changes yourself.

