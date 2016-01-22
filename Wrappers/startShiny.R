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
bioInstall <- function() {
    if("biocLite" %in% rownames(installed.packages()) == FALSE) {
        source("http://bioconductor.org/biocLite.R")
    }
    library(BiocInstaller)
    if("genomeIntervals" %in% rownames(installed.packages()) == FALSE) {
        biocLite("genomeIntervals")
    }
    library(genomeIntervals)
}

#   Start Shiny
main <- function() {
    angsdwrapper <- args[1] # Where is ANGSD-wrapper located?
    dir.create(file.path(angsd.wrapper, ".RLibs", fsep="/"), showWarnings = FALSE) # Create our directory for packages
    .libPaths(paste0(angsd.wrapper, "/.RLibs")) # Enable this as a place to look for packages
    batchInstall(pkgList) # Install dependent packages
    bioInstall() # Install dependent packages needed from BioConductor
    setwd(paste0(angsdwrapper, "/shinyGraphing")) # Set our working directory
    runApp(launch.browser = TRUE) # Run Shiny graphing
}

main()
