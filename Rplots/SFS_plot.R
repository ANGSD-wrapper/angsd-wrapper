#!/usr/bin/env Rscript

# An SFS graphing utility for angsd-wrapper
#   Enable command line arguments
args <- commandArgs(trailingOnly = TRUE)

#   A function to plot the SFS
plot.sfs <- function(sfs) {
    #   Get initial data
    #   Do not transverse data frame!
    Derived <- as.matrix(read.table(sfs, header = FALSE))   #   read in table as matrix
    Derived <- as.data.frame(t(Derived))  #   coerced to matrix before switching rows and columns
    names(Derived) <- c("Allele_frequency")   #   assign column name to dataframe
    #   Get rid of first and last values for graphing purposes
    Derived <- Derived[3:nrow(x = Derived) - 1, ]
    #   Creating SFS graphs
    #   Bar plot
    barplot(Derived, xaxt = "n", xlab = "Derived Allele Frequency", ylab = "Proportion of SNPs", main = "Site Frequency Spectrum", col = "purple4",
            offset = 0, width = 1, las = 1)
    #   Need help generalizing 'at' parameter to fit any sample size input
    axis.x <- axis(1, at= 1:23, labels = numeric(), lwd = 2)
    axis.y <- axis(side = 2, at = , labels = numeric(), lwd = 1, lwd.ticks = 1, outer = FALSE, yaxt = "n")
    #   Histogram
    #   Make histogram excluding row 1 and 25, this needs to be fixed in the wrapper script later on
    hist(Derived, right = FALSE, xlab = "# of variants at site", ylab = "Proportion of SNPs",main = "Site Frequency Spectrum", col = "deepskyblue")
    #   Plot
    plot(Derived, type = "p", xlab = "# of variants at site", ylab = "Proportion of SNPs", main = "Site Frequency Spectrum", col = "navyblue", pch = 18)
}

#   Do the work here
main <- function() {
    inputData <- args[1]
    outputName <- paste0(args[2], '/', args[3], '_SFS.pdf')
    pdf(file = outputName, width = 6, height = 6)
    options(scipen=5)
    plot.sfs(sfs = inputData)
    dev.off()
}

main()
