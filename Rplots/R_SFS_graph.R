# An SFS graphing utility for angsd-wrapper

#   Get initial data
#   Do not transverse data frame!
#   The following dataset is my Inversions dataset including 12 samples
setwd(dir = "/Users/chaochihliu/Dropbox/ANGSD Wrapper/SFS/")
Derived <- as.matrix(read.table("ANGSD_DerivedSFS", header = FALSE))   #   read in table as matrix
Derived <- as.data.frame(t(Derived))  #   coerced to matrix before switching rows and columns
names(Derived) <- c("Allele_frequency")   #   assign column name to dataframe

#   Creating SFS graphs
#   Bar plot
barplot(Derived[2:24,], xlab = "# of variants at site", ylab = "Proportion of SNPs", main = "Site Frequency Spectrum",
        col = "purple4")

#   Histogram
#   Make histogram excluding row 1 and 25, this needs to be fixed in the wrapper script later on
hist(Derived[2:24,], right = FALSE, xlab = "# of variants at site", 
     ylab = "Proportion of SNPs",main = "Site Frequency Spectrum", col = "deepskyblue")

#   Plot
plot(Derived[2:24,], type = "p", xlab = "# of variants at site",
     ylab = "Proportion of SNPs", main = "Site Frequency Spectrum", col = "navyblue", pch = 18)
  
