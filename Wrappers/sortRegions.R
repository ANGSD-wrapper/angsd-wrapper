#!/usr/bin/Rscript

args <- commandArgs(trailingOnly = TRUE)

#   A function to read in the regions of the reference file
readFai <- function(infile){
    print(x = paste('Reading in fai file', infile))
    #   A vector containing the column names for Fai files
    fai.names <- c('SeqID', 'Length', 'Offset', 'Linebases', 'Linewidth')
    #   Read in the reference file, must be tab-delimited
    reference <- read.table(file = infile, header = FALSE, sep = '\t', as.is = TRUE)
    #   Name the columns of the reference Fai file
    names(x = reference) <- fai.names
    #   Return the reference regions file
    return(reference)
}

#   A function to split a string by a delimeter and
splitString <- function(input, delim = ' ', position = 1){
    #   Split the input string by the delimeter
    split <- strsplit(x = as.character(x = input), split = as.character(x = delim))
    #   Collect which part of the split string is being asked for
    string <- unlist(x = split)[position]
    #   Return the string of interest
    return(string)
}

#   A function to read in a provided regions file for ANGSD and SAMtools
readRegions <- function(infile){
    print(x = paste('Reading in regions file', infile))
    #   Read in the regions file, spilt by colon
    regions <- read.table(file = infile, header = FALSE, sep = ':', as.is = TRUE)
    #   If we only have contig names, no subregion
    if(ncol(x = regions) == 1){
        #   Give it a name, and be done
        names(x = regions) <- 'Contig'
    }
    #   If we have two columns, as the subregion is split by a '-' instead of a ':'
    else if(ncol(x = regions) == 2){
        #   Collect start and end positions
        start <- sapply(X = regions$V2, FUN = splitString, delim = '-', position = 1)
        end <- sapply(X = regions$V2, FUN = splitString, delim = '-', position = 2)
        #   Create a temporary data frame with the correct columns
        tmp <- data.frame(Contig = as.character(x = regions$V1), Start = as.numeric(x = start), End = as.numeric(x = end))
        #   Replace the original data frame with the temporary one
        regions <- tmp
    }
    #   If we don't have either 1 or 2
    else{
        #   Wrong regions file
        stop("Invalid regions file")
    }
    #   Return our regions data frame
    return(regions)
}

#   A function to sort by start positino
sortContig <- function(contig, regions){
    #   Find the positions where this contig occurs
    contig.positions <- which(x = as.character(x = regions$Contig) %in% contig)
    #   Order these positions by start
    ordered.positions <- contig.positions[order(regions[contig.positions, ]$Start)]
    #   Return the ordered positions
    return(ordered.positions)
}

#   A function to collapse regions data into one column
collapseRegion <- function(region){
    #   Paste the values from the row together
    collapsed <- paste0(region['Contig'], ':', region['Start'], '-', region['End'], collapse = '')
    #   Return the collapsed string
    return(collapsed)
}

#   A function to reverse a string
reverseString <- function(string){
    #   Split the string into individual characters
    split <- strsplit(x = string, split = '')[[1]]
    #   Reverse the split string
    reversed <- rev(x = split)
    #   Collapse the reversed string
    collapsed <- paste(rev, collapse = '')
    return(collapsed)
}

#   A function to write out a regions file
writeRegions <- function(regions, input.name){
    #   Create an output name based off the input name
    input.dir <- dirname(path = input.name)
    input.base <- basename(path = input.name)
    input.reversed <- reverseString(string = input.base)
    no.extension <- reverseString(string = unlist(x = strsplit(x = input.base, split = '.'))[2])
    output.name <- paste0(input.dir, '/', no.extension, '_sorted.txt')
    #   Write out the table
    print(x = paste('Writing sorted regions to', output.name))
    write.table(x = regions, file = output.name, quote = FALSE, row.names = FALSE, col.names = FALSE)
}

#   Run the program
main <- function(){
    #   Collect the arguments
    regions.file <- args[1]
    fai.file <- args[2]
    #   Read in the data
    regions <- readRegions(infile = regions.file)
    fai <- readFai(infile = fai.file)
    #   Sort the data
    ordered.positions <- unlist(x = lapply(X = fai$SeqID, FUN = sortContig, regions = regions))
    ordered.regions <- regions[ordered.positions, ]
    #   Write the sorted regions data to an output file
    if(ncol(x = regions) == 3){ # Collapse the regions file down to one column if necessary
        collapsed <- apply(X = regions, MARGIN = 1, FUN = collapseRegion)
        print(x = 'Collapsing regions')
        regions <- collapsed
    }
    writeRegions(regions = regions, input.name = regions.file)
}

main()
