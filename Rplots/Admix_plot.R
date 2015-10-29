#!/usr/bin/env Rscript

# An Admixture graphing utility for angsd-wrapper
#   Enable command line arguments
args <- commandArgs(trailingOnly = TRUE)

#   A function to plot Admixture
plot.admix <- function(admixFile) {
    #   Read in the data
    admix <- t(as.matrix(read.table(admixFile)))
    k <- regmatches(basename(admixFile, regexpr('([0-9]+)', basename(name))))
    barplot(admix, col = 1:3, space = 0, xlab = 'Individuals', ylab = 'Admixture', main = paste('Admixture assuming', k, 'populations', sep = " "))
}

#   Determine the number of populations for each sample
# find.k <- function(fileList) {
#     admixFiles <- read.table(fileList, col.names = 'Qopt Files', colClasses = 'character')
#     k <- sapply(admixFiles[, 1], function(name) regmatches(basename(name), regexpr('([0-9]+)', basename(name))))
#     admixFiles["K"] = k
#     return(admixFiles)
# }

main <- function() {
    inputList <- args[1]
    outputName <- paste0(args[2], '/', args[3], '_Admixture.pdf')
    admixList <- read.table(fileList, col.names = 'Qopt Files', colClasses = 'character')
    pdf(file = outputName, with = 6, height = 6)
    options(scipen=5)
    sapply(admixList[, 1], plot.admix)
    dev.off()
}

main()