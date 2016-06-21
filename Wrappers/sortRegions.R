#!/usr/bin/Rscript

args <- commandArgs(trailingOnly = TRUE)

#   A function to split a string by a delimeter and
splitString <- function(input, delim = ' ', position = 1){
    #   Split the input string by the delimeter
    split <- strsplit(x = as.character(x = input), split = as.character(x = delim))
    #   Collect which part of the split string is being asked for
    string <- unlist(x = split)[position]
    #   Return the string of interest
    return(string)
}

#   A function to reverse a string
reverseString <- function(string){
    #   Split the string into individual characters
    split <- strsplit(x = string, split = '')[[1]]
    #   Reverse the split string
    reversed <- rev(x = split)
    #   Collapse the reversed string
    collapsed <- paste(reversed, collapse = '')
    return(collapsed)
}

#   A function to get the basename of a string without any extension
baseName <- function(string){
    #   Get the base name
    base <- basename(path = string)
    #   Reverse the string and get the extension of the file
    base.reversed <- reverseString(string = base)
    base.extension <- paste0('.', reverseString(string = splitString(input = base.reversed, delim = '\\.', position = 1)))
    extension.position <- regexec(pattern = base.extension, text = base)[[1]][1]
    #   Split the base name for subsetting
    base.split <- strsplit(x = base, split = '')[[1]]
    #   Subset the base name to remove the extension and collapse into single character string
    base.extensionless <- paste0(base.split[1:extension.position - 1], collapse = '')
    return(base.extensionless)
}

#   A function to read in the regions of the reference file
readFai <- function(infile){
    cat('Reading in fai file', infile, '\n', sep = ' ')
    #   A vector containing the column names for Fai files
    fai.names <- c('SeqID', 'Length', 'Offset', 'Linebases', 'Linewidth')
    #   Read in the reference file, must be tab-delimited
    reference <- read.table(file = infile, header = FALSE, sep = '\t', as.is = TRUE)
    #   Name the columns of the reference Fai file
    names(x = reference) <- fai.names
    #   Return the reference regions file
    return(reference)
}

#   A function to read in a provided regions file for ANGSD and SAMtools
readRegions <- function(infile){
    cat('Reading in regions file', infile, '\n', sep = ' ')
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
    if(ncol(x = regions) == 3){
        ordered.positions <- contig.positions[order(regions[contig.positions, ]$Start)]
    }
    #   Return the ordered positions
    return(ordered.positions)
}

#   A function to collapse regions data into one column
collapseRegion <- function(region){
    #   Remove any whitespace from each entry
    contig <- gsub(pattern = "[[:space:]]", replacement = '', x = as.character(x = region[['Contig']]))
    start <- gsub(pattern = "[[:space:]]", replacement = '', x = as.character(x = region[['Start']]))
    end <- gsub(pattern = "[[:space:]]", replacement = '', x = as.character(x = region[['End']]))
    #   Paste the values from the row together
    collapsed <- paste0(contig, ':', start, '-', end, collapse = '')
    #   Return the collapsed string
    return(collapsed)
}

#   A function to write out a regions file
writeRegions <- function(regions, input.name){
    #   Create an output name based off the input name
    input.dir <- dirname(path = input.name)
    input.base <- baseName(string = input.name)
    output.name <- paste0(input.dir, '/', input.base, '_sorted.txt')
    #   Write out the table
    cat('Writing sorted regions to', output.name, '\n', sep = ' ')
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
    cat('Sorting the regions file', '\n', sep = ' ')
    ordered.positions <- unlist(x = lapply(X = unique(x = fai$SeqID), FUN = sortContig, regions = regions))
    ordered.regions <- regions[ordered.positions, ]
    row.names(x = ordered.regions) <- NULL
    #   Write the sorted regions data to an output file
    if(ncol(x = ordered.regions) == 3){ # Collapse the regions file down to one column if necessary
        cat('Collapsing regions', '\n', sep = ' ')
        collapsed <- data.frame(Contig = apply(X = ordered.regions, MARGIN = 1, FUN = collapseRegion))
        ordered.regions <- collapsed
    }
    writeRegions(regions = ordered.regions, input.name = regions.file)
}

main()
