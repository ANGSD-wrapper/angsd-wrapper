#!/usr/bin/env Rscript

#   Start the Shiny graphics interface
args <- commandArgs(trailingOnly = TRUE)

#   Test to see if required packages are installed
pkgTest <- function(package) {
    if(package %in% rownames(installed.packages()) == FALSE) { # check to see if a packages is available
        install.packages(package) # if not, insall it
    }
}

#   A list of required pacakges
pkgList <- c("shiny", "lattice", "Hmisc", "ape", "data.table", "DT", "shinythemes")

#   Load the packages
batchInstall <- function(pkgList) {
    options(repos = c(CRAN = "http://cran.rstudio.com")) # set a repo mirror, we used RStudio just because
    for(dep in pkgList) {
        pkgTest(dep) # test to see if the package is installed
    }
    lapply(X = pkgList, FUN = library, character.only = TRUE) # load the packages to be used
}

#   Install packages from BioConductor
# biocLite is depricated in Bioconductor 3.10 and R 3.5 or greater. Need ot use BiocManager to install Bioconductor packages
bioInstall <- function() {
    #if("biocLite" %in% rownames(installed.packages()) == FALSE) {
        #source("http://bioconductor.org/biocLite.R")
    #}
    if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
   # BiocManager::install("BiocInstaller")
    #library(BiocInstaller)
    if("genomeIntervals" %in% rownames(installed.packages()) == FALSE) {
        #biocLite("genomeIntervals")
        BiocManager::install("genomeIntervals")
    }
    library(genomeIntervals)
}

#   Start Shiny
main <- function() {
    angsdwrapper <- args[1] # Where is ANGSD-wrapper located?
    dir.create(file.path(angsdwrapper, ".RLibs", fsep="/"), showWarnings = FALSE) # Create our directory for packages
    .libPaths(paste0(angsdwrapper, "/.RLibs")) # Enable this as a place to look for packages
    batchInstall(pkgList) # Install dependent packages
    bioInstall() # Install dependent packages needed from BioConductor
    setwd(paste0(angsdwrapper, "/shinyGraphing")) # Set our working directory
    runApp(launch.browser = TRUE) # Run Shiny graphing
}

main()
