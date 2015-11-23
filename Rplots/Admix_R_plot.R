# Read in admixture output data
# This works
admix <- fread("C:/Users/Chaochih/Dropbox/ANGSD_Wrapper/output_barley2_PH/output_barley2_PH.allK.qopt", header = FALSE, sep = "\n", na.strings = c("NA"))

admix <- fread("C:/Users/Chaochih/Dropbox/ANGSD_Wrapper/Inverted_2DSFS.qopt", header = FALSE, sep = "\n", na.strings = c("NA"))

admix <- read.table("C:/Users/Chaochih/Dropbox/ANGSD_Wrapper/Inverted_2DSFS.qopt", header = FALSE, fill = TRUE)

# This set works
# R reads in each row of values as a single column, need to figure out how to get it to read it in as multiple columns
admix.2pop <- admix[1:12, ]
admix.3pop <- admix[13:24, ]
admix.4pop <- admix[25:36, ]
admix.5pop <- admix[37:48, ]

setnames(admix.2pop, admix.headers.5)

# Convert data to matrix
admix.2dat <- t(as.data.frame(as.matrix(admix.2pop)))
admix.3dat <- t(as.data.frame(as.matrix(admix.3pop)))
admix.4dat <- t(as.data.frame(as.matrix(admix.4pop)))
admix.5dat <- t(as.data.frame(as.matrix(admix.5pop)))

# Graph admixture
barplot(admix.2dat,
        col=c("#006BA4","#FF800E","#A2C8EC", "#898989","#ABABAB",
              "#595959", "#5F9ED1","#CFCFCF","#FFBC79","#C85200"),
        space=0, width = 12,
        border=NA,  
        xlab="Individuals", 
        ylab="admixture proportion")

# Specify headers
# Will most likely have to set headers depending on the # of K ancestral populations
admix.headers.2 <- c("Population_1", "Population_2")
admix.headers.3 <- c("Population_1", "Population_2", "Population_3")
admix.headers.4 <- c("Population_1", "Population_2", "Population_3", "Population_4")
admix.headers.5 <- c("Population_1", "Population_2", "Population_3", "Population_4", "Population_5")
admix.headers.6 <- c("Population_1", "Population_2", "Population_3", "Population_4", "Population_5", "Population_6")
admix.headers.7 <- c("Population_1", "Population_2", "Population_3", "Population_4", "Population_5", "Population_6", "Population_7")
admix.headers.8 <- c("Population_1", "Population_2", "Population_3", "Population_4", "Population_5", "Population_6", "Population_7", "Population_8")
admix.headers.9 <- c("Population_1", "Population_2", "Population_3", "Population_4", "Population_5", "Population_6", "Population_7", "Population_8", "Population_9")
admix.headers.10 <- c("Population_1", "Population_2", "Population_3", "Population_4", "Population_5", "Population_6", "Population_7", "Population_8", "Population_9", "Population_10")