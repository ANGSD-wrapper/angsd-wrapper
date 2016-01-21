#!/usr/bin/env Rscript

# Collect arguments from the command line
args <- commandArgs(trailingOnly = TRUE)
intersectOut <- args[1]
fstOut <- args[2]
mafs.g1 <- args[3]
mafs.g2 <- args[4]
graph.me.fst <- args[5]

# Install dependent libraries
pkgTest <- function(package) {
    if(package %in% rownames(installed.packages()) == FALSE) { # check to see if a packages is available
        install.packages(package) # if not, insall it
    }
}

#   A list of required packages
pkgList <- c("data.table", "dplyr")

#   Load the packages
options(repos = c(CRAN = "http://cran.rstudio.com")) # set a repo mirror, we used RStudio just because
for(dep in pkgList) { # Iterate through each dependent package
    pkgTest(dep) # Test to see if it's installed
}
lapply(X = pkgList, FUN = library, character.only = TRUE) # load the packages to be used

# Define headers
fst.headers <- c("A", "AB", "f", "FST", "Pvar")
intersect.headers <- c("Chr","position")

# Read in data from .fst file, shared.pos file, and two mafs files from 2 selected populations
# Subset data and select columns to pull
intersect <- fread(input = intersectOut, sep = "\t", col.names = intersect.headers)
fst <- fread(input = fstOut, sep = "\t", col.names = fst.headers)
fst <- subset(x = fst, select = c("Pvar"))
mafs1 <- fread(input = mafs.g1)
mafs1 <- subset(x = mafs1, select = c("chromo", "position", "major", "minor", "knownEM"))
mafs2 <- fread(input = mafs.g2)
mafs2 <- subset(x = mafs2, select = c("chromo", "position", "major", "minor", "knownEM"))

# Merge mafs1 and mafs2
mafs.1and2 <- merge(mafs1, mafs2, by = "position", suffixes = c(".mafs1", ".mafs2"))

# Combine two datasets
FST <- cbind(intersect, fst)

# Merge FST dataset and mafs.1and2 dataset and filter by columns
data.fst <- merge(mafs.1and2, FST, by = "position")
filter.knownEM <- subset(x = data.fst, knownEM.mafs1 + knownEM.mafs2 > 0)
filter.minor <- subset(x = filter.knownEM, minor.mafs2 == minor.mafs1 | minor.mafs2 == major.mafs1)
filter.major <- subset(x = filter.minor, major.mafs2 == major.mafs1 | major.mafs2 == minor.mafs1)

# Save to .fst file
write.table(x = filter.major, file = "maize.graph.me.fst", col.names = F, row.names = F)