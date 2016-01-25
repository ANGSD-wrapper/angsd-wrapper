#!/usr/bin/env Rscript

# Collect arguments from the command line
args <- commandArgs(trailingOnly = TRUE)
angsd.wrapper <- args[1]
thetasOut <- args[2]
graph.me.thetas <- args[3]

# Set install directory for packages
dir.create(file.path(angsd.wrapper, ".RLibs", fsep="/"), showWarnings = FALSE)
.libPaths(paste0(angsd.wrapper, "/.RLibs"))

# Install dependent libraries
pkgTest <- function(package) {
    if(package %in% rownames(installed.packages()) == FALSE) { # check to see if a packages is available
        install.packages(package) # if not, insall it
    }
}

#   A list of required packages
pkgList <- c("data.table")

#   Load the packages
options(repos = c(CRAN = "http://cran.rstudio.com")) # set a repo mirror, we used RStudio just because
for(dep in pkgList) { # Iterate through each dependent package
    pkgTest(dep) # Test to see if it's installed
}
lapply(X = pkgList, FUN = library, character.only = TRUE) # load the packages to be used

# Define headers
thetas.headers <- c("(indexStart,indexStop)(firstPos_withData,lastPos_withData)(WinStart,WinStop)", "Chr","WinCenter", "tW", "tP", "tF", "tH", "tL", "Tajima", "fuf", "fud", "fayh", "zeng", "nSites")

# Read in data from .pestPG file
thetas <- fread(input = thetasOut, sep = "\t", col.names = thetas.headers)

# Filter invariant sites based on Tajima column
filter.thetas <- subset(x = thetas, Tajima != 0)

# Save to .thetas file
write.table(x = filter.thetas, file = graph.me.thetas, col.names = T, row.names = F)