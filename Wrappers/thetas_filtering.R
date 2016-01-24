#!/usr/bin/env Rscript

# Collect arguments from the command line
args <- commandArgs(trailingOnly = TRUE)
thetasOut <- args[1]
graph.me.thetas <- args[2]

# Define headers
thetas.headers <- c("(indexStart,indexStop)(firstPos_withData,lastPos_withData)(WinStart,WinStop)", "Chr","WinCenter", "tW", "tP", "tF", "tH", "tL", "Tajima", "fuf", "fud", "fayh", "zeng", "nSites")

# Read in data from .pestPG file
thetas <- fread(input = thetasOut, sep = "\t", col.names = thetas.headers)

# Filter invariant sites based on Tajima column
filter.thetas <- subset(x = thetas, Tajima != 0)

# Save to .thetas file
write.table(x = filter.thetas, file = graph.me.thetas, col.names = T, row.names = F)