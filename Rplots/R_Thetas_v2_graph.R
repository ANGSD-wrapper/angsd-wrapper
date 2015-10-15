#!/usr/bin/env Rscript

#   A Thetas graphing utility for angsd-wrapper
#   Enable command line arguments
args <- commandArgs(trailingOnly = TRUE)

#   Get initial data
pest.pg <- function(thetas) {
  #   read in table as data frame
  Pest <- as.data.frame(thetas)
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
  #   Remove duplicate contigs
  chromoNames <- list(as.character(unique(pest$Chromosome)))
}

#   Watterson Estimates by base pair
#   Creating new function called Watterson.chr
Watterson.chr <- function(chromosome, thetas) {
  #   subset and select columns to put into new data frame
  chromosome <- Pest$Chromosome
  chromo <- subset(x = thetas, select = c("nSites", "tW"))
  totalWatterson <- sum(chromo$Watterson)
}

