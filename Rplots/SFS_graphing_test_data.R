# An SFS graphing utility for angsd-wrapper

#   Get initial data
#   The following is the test dataset from Arun
#   Do not transverse data frame!
setwd(dir = "/Users/chaochihliu/Dropbox/ANGSD Wrapper/angsd-wrapper_test/SFS/")
Derived_test <- as.matrix(read.table("ANGSD_test_DerivedSFS", header = FALSE))   #   read in table as matrix
Derived_test <- as.data.frame(t(Derived_test))  #   coerced to matrix before switching rows and columns
names(Derived) <- c("Allele_frequency")   #   assign column name to dataframe

#   Creating SFS graphs
#   Bar plot
barplot(Derived_test[2:24,], xlab = "Sites", ylab = "Allele frequency", main = "Site Frequency Spectrum",
        col = "purple4")

#   Histogram
#   Make histogram excluding row 1 and 25, this needs to be fixed in the wrapper script later on
hist(Derived_test[2:24,], right = FALSE, xlab = "Sites", 
     ylab = "Allele frequency",main = "Site Frequency Spectrum", col = "deepskyblue")

#   Plot
plot(Derived_test[2:24,], type = "p", xlab = "Sites",
     ylab = "Allele frequency", main = "Site Frequency Spectrum", col = "navyblue", pch = 18)


