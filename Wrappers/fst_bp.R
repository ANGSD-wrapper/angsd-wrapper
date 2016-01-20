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

#   A list of required pacakges
pkgList <- c("data.table", "dplyr")

#   Load the packages
options(repos = c(CRAN = "http://cran.rstudio.com")) # set a repo mirror, we used RStudio just because
for(dep in pkgList) { # Iterate through each dependent package
    pkgTest(dep) # Test to see if it's installed
}
lapply(X = pkgList, FUN = library, character.only = TRUE) # load the packages to be used

# Define headers
fst.headers <- c("A", "AB", "f", "FST", "Pvar")
intersect.headers <- c("Chr","bp")

# Read in data from .fst file and shared.pos file
intersect <- read.table(file = intersectOut, sep = "\t", col.names = intersect.headers)
fst <- read.table(file = fstOut, sep = "\t", col.names = fst.headers %>%
               select(c("Pvar")))
mafs1 <- read.table(file = mafs.g1 %>%
                 select(c(chromo, position, major, minor, knownEM)))
mafs2 <- read.table(file= = mafs.g2 %>%
                 select(c(chromo, position, major, minor, knownEM)))

# Combine two datasets
FST <- cbind(intersect, fst, mafs1, mafs2)

# Save to .fst file
write.table(x = FST, file = graph.me.fst, col.names = F, row.names = F)
