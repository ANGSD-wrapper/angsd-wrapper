#!/usr/bin/env Rscript

#   A Thetas graphing utility for angsd-wrapper
#   Enable command line arguments
args <- commandArgs(trailingOnly = TRUE)

#   Get initial data
pest.pg <- function(thetas) {
  #   read in table as data frame
  Pest <- read.table(thetas, header = FALSE)
  #   Assign column names to dataframe
  names(Pest) <- c("#", "Chromosome", "Wincenter", "tW", "tP", "tF", "tH", "tL", "Tajima's D", "fuf", "fud", "fayh", "zeng", "nSites")
  #   Subset data frames
  pest <- data.frame(Chromosome = Pest$Chromosome)
  pest$WinCenter <- Pest$Wincenter
  pest$nSites <- Pest$nSites
  pest$Watterson <- Pest$tW
  pest$Pairwise <- Pest$tP
  pest$Tajima <- Pest$"Tajima's D"
  return(pest)
}

#   Watterson Estimates by base pair
#   Creating new function called Watterson.chr
Watterson.chr <- function(chromosome, thetas) {
  #   subset and select columns to put into new data frame
  chromo <- subset(x = thetas, subset = Chromosome == chromosome, select = c("nSites", "Watterson"))
  totalWatterson <- sum(chromo$Watterson)
  divWatterson <- totalWatterson / chromo$nSites
  return(divWatterson)
}

#   Pairwise Estimates by base pair
#   Creating a new function
Pairwise.chr <- function(chromosome, thetas) {
  chromo <- subset(x = thetas, subset = Chromosome == chromosome, select = c("nSites", "Pairwise"))
  totalPairwise <- sum(chromo$Pairwise)
  divPairwise <- totalPairwise / chromo$nSites
  return(divPairwise)
}

#   Do the work here
main <- function() {
  #   Where is our data?
  inputData <- args[1]
  #   Create a name for the output file
  outputName <- past0(args[2], '/', args[3], '_Thetas.pdf')
  #   Read the data in as a dataframe
  pest <- pest.pg(thetas = inputData)
  #   Remove duplicate contigs
  chromoNames <- list(as.character(unique(pest$Chromosome)))
  #   Find the Watterson Estimates by basepair
  wattersons <- sapply(X = chromoNames[[1]], FUN = Watterson.chr, thetas = pest)
  #   Find Pairwise Estimates by base pair
  pairwises <- sapply(X = chromoNames[[1]], FUN = Pairwise.chr, thetas = pest)
  #   Estimated Watterson's Theta, Pairwise Theta, and TajD per base pair per region
  estTheta <- data.frame(Watterson = wattersons, Pairwise = pairwises, "Tajima's D" = pest$Tajima)
  return(wattersons, pairwises, estTheta)
}

# Not needed
# #   Summary Statistics for Watterson's Theta, Pairwise Theata, and Tajima's D
# summary(estTheta)

# Not needed
# #   Variance and Average of Watterson's Theta
# var(estTheta$Watterson)
# mean(estTheta$Watterson)

# Not needed
# #   Variance and Average of Pairwise Theta
# var(estTheta$Pairwise)
# mean(estTheta$Pairwise)

# Not needed
# #   Variance and Average of Tajima's D
# var(estTheta$Tajima)
# mean(estTheta$Tajima)

# Make these a function
#   Graphs: Histogram, Scatter Plot, and Barplot
#   Watterson's Theta
hist(Watterson, right = FALSE, xlab = "Watterson's Theta",
     main = "Frequency Distribution of Watterson's Theta", col = "lightsteelblue")

plot(x = pest$WinCenter, y = estTheta$Watterson, type = "p", xlab = "Position (bp)",
     ylab = "Watson's Theta", main = "Estimators of Watson's Theta",
     col = "lightsteelblue", pch = 19)

barplot(wattersons, col = sample(x= colors(distinct = TRUE),
                                 size = length(wattersons), replace = FALSE),
        main = "Estimated Watterson's Theta per base pair",
        xlab = "Region", ylab = "Estimated Theta", las = 2, names.arg = FALSE)

#   Pairwise Theta
hist(pairwises, right = FALSE, xlab = "Pairwise Theta",
     main = "Frequency Distribution of Pairwise Theta", col = "seagreen2")

plot(x = pest$WinCenter, y = estTheta$Pairwise, type = "p", xlab = "Position (bp)",
     ylab = "Pairwise Theta", main = "Estimators of Pairwise Theta",
     col = "seagreen2", pch = 19)

barplot(pairwises, col = sample(x= colors(distinct = TRUE), size = length(pairwises),
                                replace = FALSE), main = "Estimated Pairwise Theta per base pair",
        xlab = "Region", ylab = "Estimated Theta", las = 2, names.arg = FALSE)

#   Tajima's D
hist(estTheta$Tajima, right = FALSE, xlab = "Tajima's D",
     main = "Frequency Distribution of Tajima's D", col = "lightsalmon")

plot(x = pest$WinCenter, y = estTheta$Tajima.s.D, type = "p", xlab = "Position (bp)",
     ylab = "Tajima's D", main = "Estimators of Tajima's D",
     col = "lightsalmon", pch = 19)

barplot(estTheta$Tajima, col = sample(x= colors(distinct = TRUE), size = length(estTheta$T), replace = FALSE),
        main = "Estimated Tajima's D per base pair", xlab = "Region",
        ylab = "Estimated Tajima's D", las = 2, names.arg = FALSE)
