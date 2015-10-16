#!/usr/bin/env Rscript

#   A Thetas graphing utility for angsd-wrapper
#   Enable command line arguments
args <- commandArgs(trailingOnly = TRUE)

#   Get initial data
pest.pg <- function(thetas) {
  #   read in table as data frame
  Pest <- read.table(thetas, header = FALSE)
  #   Assign column names to dataframe
  names(Pest) <- c("#", "Chromosome", "Wincenter", "tW", "tP", "tF", "tH", "tL", "Tajima's D", "fuf", "fud",
                   "fayh", "zeng", "nSites")
  #   Subset data frames
  pest <- data.frame(Chromosome = Pest$Chromosome)
  pest$WinCenter <- Pest$Wincenter
  pest$nsites <- Pest$nSites
  pest$Watterson <- Pest$tW
  pest$Pairwise <- Pest$tP
  pest$Tajima <- Pest$"Tajima's D"
  return(pest)
}

#   Watterson Estimates by base pair
#   Creating new function called Watterson.chr
Watterson.chr <- function(chromosome, thetas) {
  #   subset and select columns to put into new data frame
  chromosome <- Pest$Chromosome
  chromo <- subset(x = thetas, select = c("nSites", "tW"))
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
}
