#!/usr/bin/env Rscript

# An SFS graphing utility for angsd-wrapper
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
    barplot(Derived, xlab = "# of variants at site", ylab = "Proportion of SNPs", main = "Site Frequency Spectrum", col = "purple4")
    #   Histogram
    #   Make histogram excluding row 1 and 25, this needs to be fixed in the wrapper script later on
    hist(Derived, right = FALSE, xlab = "# of variants at site", ylab = "Proportion of SNPs",main = "Site Frequency Spectrum", col = "deepskyblue")
    #   Plot
    plot(Derived, type = "p", xlab = "# of variants at site", ylab = "Proportion of SNPs", main = "Site Frequency Spectrum", col = "navyblue", pch = 18)
}
