#!/usr/bin/env Rscript
library(data.table)

# Collect arguments from the command line
args <- commandArgs(trailingOnly = TRUE)
intersectOut <- args[1]
fstOut <- args[2]
graph.me.fst <- args[3]

# Define headers
fst.headers <- c("A", "AB", "f", "FST", "Pvar")
intersect.headers <- c("Chr","bp")

# Read in data from .fst file and shared.pos file
intersect <- fread(intersectOut, sep = "\t", col.names = intersect.headers)
fst <- read.table(fstOut, sep = "\t", col.names = fst.headers)

# Combine two datasets
FST <- cbind(intersect, fst)

# Save to .fst file
write.table(FST, graph.me.fst, col.names = F, row.names = F)
