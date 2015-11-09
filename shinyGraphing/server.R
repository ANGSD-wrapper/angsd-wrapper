library(shiny)
library(genomeIntervals)
library(lattice)
library(Hmisc)
library(ape)
library(data.table)
options(shiny.maxRequestSize = -1)

# Assign headers to thetas, Fst and intersect data
thetas.headers <- c("(indexStart,indexStop)(firstPos_withData,lastPos_withData)(WinStart,WinStop)","Chr","WinCenter","tW","tP","tF","tH","tL","Tajima","fuf","fud","fayh","zeng","nSites")
fst.headers <- c("A", "AB", "f", "FST", "Pvar")
intersect.headers <- c("Chr","bp")
sfs.headers <- c("Allele Frequency")

#thetas <- fread("C:/Users/Chaochih/angsd-wrapper/shinyGraphing/BKN_Diversity.thetas.gz.pestPG", sep = "\t")

thetas <- fread("~/ANGSD-wrapper_beta_test/PH/angsd-wrapper/shinyGraphing/BKN_Diversity.thetas.gz.pestPG")

not.loaded <- TRUE

# Define server logic required to draw plots
shinyServer(

  function(input, output) {
    
    # Thetas input data
    dataInputThetas = reactive({
      data <- input$userThetas
      path <- as.character(data$datapath)
      thetas <- fread(input=path, sep="\t")
      setnames(thetas,thetas.headers)
      return(thetas)
    })
    
    # SFS input data
    dataInputSFS = reactive({
      data <- input$userSFS
      path <- as.character(data$datapath)
      #Derived <- as.matrix(read.table(input=path, header = FALSE))
      
      Derived <- as.matrix(read.table(path, header = FALSE))
      # Transverse data frame
      sfs <- as.data.frame(t(Derived))
      names(sfs) <- "Allele_frequency"
      alleles <- sfs$Allele_frequency[seq(3, nrow(sfs), by = 2)]
      return(alleles) 

    })

    # Fst input data
    dataInputFst = reactive({
      data <- input$userFst
      path <- as.character(data$datapath)
      fst <- read.table(file=path,
                        sep="\t",
                        col.names=fst.headers
      )
      return(fst)
    })

    # Intersect input data
    dataInputIntersect = reactive({
      data <- input$userIntersect
      path <- as.character(data$datapath)
      intersect <- read.table(file=path,
                              sep="\t",
                              col.names=intersect.headers
      )
      return(intersect)
    })

    # Admixture input data
    dataInputAdmix = reactive({
      data <- input$userAdmix
      path <- as.character(data$datapath)
      admix <- t(as.matrix(read.table(path)))
      return(admix)
    })
    
    # Admixture input data
    dataInputAdmix = reactive({
      data <- input$userAdmix2
      path <- as.character(data$datapath)
      admix <- t(as.matrix(read.table(path)))
      return(admix)
    })
    
    # ABBA BABA test input data
    dataInputABBABABA = reactive({
      data <- input$userABBABABA
      path <- as.character(data$datapath)
      ABBABABA <- read.table(path, sep="\t", header=T)
      return(ABBABABA)

    })

    # PCA plot input data
    dataInputPCA = reactive({
      data <- input$userPCA
      path <- as.character(data$datapath)
      PCA <- read.table(path, header=F)
      return(PCA)
    })

    # GFF file input
    gffInput = reactive({
      data <- input$userAnnotations
      path <- as.character(data$datapath)
      gff <- readGff3(path)
      return(gff)
    })
    
    
    # Thetas output data
    output$thetaChroms = renderUI({
      if(is.null(input$userThetas)){
        choices <- 10
      }
      else{
      thetas <- dataInputThetas()
      choices <- unique(thetas$Chr)
      }
      selectInput('thetaChrom', 'Chromosome to plot', 
                  choices,
                  multiple = TRUE)
    })
 
    # Create zoomable thetas plots
    ranges <- reactiveValues(x = NULL, y = NULL)
    
    # Create reactive plot on left side of page
    output$thetaPlot1 <- renderPlot({
      #  error handling code to provide a default dataset to graph
      thetas <- tryCatch({
        dataInputThetas()
      }, 
      error = function(err) {
        thetas <- read.table(file="BKN_Diversity.thetas.gz.pestPG",
                             sep="\t",
                             col.names=thetas.headers)
      }
      )

      thetas <- subset(thetas,Chr==input$thetaChrom)
      if(input$annotations){
        validate(need(input$userAnnotations, 'Need GFF file before clicking checkbox!'))
        gff <- gffInput()
        gff.gene <- subset(gff, type="gene")
        gff.df <- data.frame(gff.gene,annotation(gff))
        gff.df.chr <- subset(gff.df, seq_name==thetas$Chr[1])
        if(length(gff.df.chr$seq_name)==0){
          stop("Annotation does not match graphed region. Please make sure the first column of your GFF file matches the Chr column of the .pestPG file.")
        }
        gff.df.gene <- subset(gff.df.chr, type=="gene")
      }

      if(input$subset) {
        thetas.plot <- subset(thetas, WinCenter > input$WinCenterLow & WinCenter < input$WinCenterHigh)
      }
      else {
        thetas.plot <- thetas
      }

      # remove nsites=0
      thetas.plot <- subset(thetas.plot, nSites != 0)
      # remove data points with less than 50 sites. Calculate minimum from data?
      if(input$rm.nsites) {
        thetas.plot <- subset(thetas.plot, nSites > input$nsites)
      }
      #divide thetas by the number of sites in each window
      thetas.plot$tW <- thetas.plot$tW/thetas.plot$nSites
      thetas.plot$tP <- thetas.plot$tP/thetas.plot$nSites
      thetas.plot$tF <- thetas.plot$tF/thetas.plot$nSites
      thetas.plot$tH <- thetas.plot$tH/thetas.plot$nSites
      thetas.plot$tL <- thetas.plot$tW/thetas.plot$nSites

      data <- switch(input$thetaChoice,
                     "Watterson's Theta" = thetas.plot$tW,
                     "Pairwise Theta" = thetas.plot$tP,
                     "Fu and Li's Theta" = thetas.plot$tF,
                     "Fay's Theta" = thetas.plot$tH,
                     "Maximum likelihood (L) Theta" = thetas.plot$tL
      )
      if(input$annotations) {
        plot(thetas.plot$WinCenter,
             data, t="p", pch=19,col=rgb(0,0,0,0.5),
             xlab="Position (bp)",
             ylab=paste(input$thetaChoice,"Estimator Value"),
             main=paste("Estimators of theta along chromosome", thetas$Chr[1])
        )

        # Different representation of the data to plot
        rug(rect(gff.df.gene$X1, -1e2, 
                 gff.df.gene$X2, 0, 
                 col=rgb(0.18,0.55,0.8,0.75), 
                 border=NA))
        if(input$thetaLowess){lines(lowess(thetas.plot$WinCenter,data, f=0.1), 
                                    col="red")}
      }
      else {
        plot(thetas.plot$WinCenter,
             data, t="p", pch=19,col=rgb(0,0,0,0.5),
             xlab="Position (bp)",
             ylab=paste(input$thetaChoice,"Estimator Value"),
             main=paste("Estimators of theta along chromosome", thetas$Chr[1])
        )
        if(input$thetaLowess){lines(lowess(thetas.plot$WinCenter,data, f=0.1), 
                                    col="red")}
      }
    })
    
    # Create reactive plot that zooms in on right side
    output$thetaPlot2 <- renderPlot({
      # error handling code to provide a default dataset to graph
      thetas <- tryCatch({
        dataInputThetas()
      }, error = function(err) {
        thetas <- read.table(file="BKN_Diversity.thetas.gz.pestPG",
                             sep="\t",
                             col.names=thetas.headers)
      }
      )
      
      thetas <- subset(thetas,Chr==input$thetaChrom)
      if(input$annotations){
        validate(need(input$userAnnotations, 'Need GFF file before clicking checkbox!'))
        gff <- gffInput()
        gff.gene <- subset(gff, type="gene")
        gff.df <- data.frame(gff.gene,annotation(gff))
        gff.df.chr <- subset(gff.df, seq_name==thetas$Chr[1])
        if(length(gff.df.chr$seq_name)==0){
          stop("Annotation does not match graphed region. Please make sure the first column of your GFF file matches the Chr column of the .pestPG file.")
        }
        gff.df.gene <- subset(gff.df.chr, type=="gene")
      }
      
      if(input$subset) {
        thetas.plot <- subset(thetas, WinCenter > input$WinCenterLow & WinCenter < input$WinCenterHigh)
      }
      else {
        thetas.plot <- thetas
      }
      
      # remove nsites=0
      thetas.plot <- subset(thetas.plot, nSites != 0)
      # remove data points with less than 50 sites. Calculate minimum from data?
      if(input$rm.nsites) {
        thetas.plot <- subset(thetas.plot, nSites > input$nsites)
      }
      #divide thetas by the number of sites in each window
      thetas.plot$tW <- thetas.plot$tW/thetas.plot$nSites
      thetas.plot$tP <- thetas.plot$tP/thetas.plot$nSites
      thetas.plot$tF <- thetas.plot$tF/thetas.plot$nSites
      thetas.plot$tH <- thetas.plot$tH/thetas.plot$nSites
      thetas.plot$tL <- thetas.plot$tW/thetas.plot$nSites
      
      # Choose one of the following data sets to plot
      data <- switch(input$thetaChoice,
                     "Watterson's Theta" = thetas.plot$tW,
                     "Pairwise Theta" = thetas.plot$tP,
                     "Fu and Li's Theta" = thetas.plot$tF,
                     "Fay's Theta" = thetas.plot$tH,
                     "Maximum likelihood (L) Theta" = thetas.plot$tL
      )
      if(input$annotations) {
        plot(thetas.plot$WinCenter,
             data, t="p", pch=19,col=rgb(0,0,0,0.5),
             xlab="Position (bp)",
             ylab=paste(input$thetaChoice,"Estimator Value"),
             main=paste("Estimators of theta along chromosome", thetas$Chr[1])
        )
        
        rug(rect(gff.df.gene$X1, -1e2, 
                 gff.df.gene$X2, 0, 
                 col=rgb(0.18,0.55,0.8,0.75), 
                 border=NA))
        if(input$thetaLowess){lines(lowess(thetas.plot$WinCenter, data, f=0.1), 
                                    col="red")}
      }
      else {
        plot(thetas.plot$WinCenter,
             data, t="p", pch=19,col=rgb(0,0,0,0.5),
             xlab="Position (bp)",
             ylab=paste(input$thetaChoice,"Estimator Value"),
             main=paste("Estimators of theta along chromosome",
                        thetas$Chr[1]),
             xlim = ranges$x, ylim = ranges$y
        )
        if(input$thetaLowess){lines(lowess(thetas.plot$WinCenter,data, f=0.1), 
                                    col="red")}
      }
    })

    # Creating zoom function in thetasPlot1 and thetasPlot2
    observe({
      brush <- input$thetaPlot1_brush
      if (!is.null(brush)) {
        ranges$x <- c(brush$xmin, brush$xmax)
        ranges$y <- c(brush$ymin, brush$ymax)
        
      } else {
        ranges$x <- NULL
        ranges$y <- NULL
      }
    })
    
    # Creating reactive plot for selectionPlot1 and selectionPlot2
    ranges2 <- reactiveValues(x = NULL, y = NULL)
    
    # selectionPlot1
    output$selectionPlot1 <- renderPlot({
      # error handling code to provide a default dataset to graph
      thetas <- tryCatch({
        dataInputThetas()
      }, error = function(err) {
        thetas <- read.table(file="BKN_Diversity.thetas.gz.pestPG",
                             sep="\t",
                             col.names=thetas.headers
        )
      }
      )

      thetas <- subset(thetas,Chr==input$thetaChrom)

      if(input$subset) {
        thetas.plot <- subset(thetas, WinCenter > input$WinCenterLow & WinCenter < input$WinCenterHigh)
      }
      else {
        thetas.plot <- thetas
      }

      # remove nsites=0
      thetas.plot <- subset(thetas.plot, nSites != 0)

      data <- switch(input$selectionChoice,
                     "Tajima's D" = thetas.plot$Tajima,
                     "Fu and Li's F" = thetas.plot$fuf,
                     "Fu and Li's D" = thetas.plot$fud,
                     "Fay and Wu's H" = thetas.plot$fayh,
                     "Zeng's E" = thetas.plot$zeng
      )
      
      plot(thetas.plot$WinCenter,
           data, t="p", pch=19,col=rgb(0,0,0,0.5),
           xlab="Position (bp)",
           ylab=input$selectionChoice,
           main=paste("Neutrality test statistics along chromosome", thetas$Chr[1])
      )
      if(input$selectionLowess){lines(lowess(thetas.plot$WinCenter,data, f=0.1), col="red")}
    })
    
    # selectionPlot2
    output$selectionPlot2 <- renderPlot({
      # error handling code to provide a default dataset to graph
      thetas <- tryCatch({
        dataInputThetas()
      }, error = function(err) {
        thetas <- read.table(file="BKN_Diversity.thetas.gz.pestPG",
                             sep="\t",
                             col.names=thetas.headers
        )
      }
      )
      
      
      thetas <- subset(thetas,Chr==input$thetaChrom)
      
      if(input$subset) {
        thetas.plot <- subset(thetas, WinCenter > input$WinCenterLow & WinCenter < input$WinCenterHigh)
      }
      else {
        thetas.plot <- thetas
      }
      
      # remove nsites=0
      thetas.plot <- subset(thetas.plot, nSites != 0)
      
      data <- switch(input$selectionChoice,
                     "Tajima's D" = thetas.plot$Tajima,
                     "Fu and Li's F" = thetas.plot$fuf,
                     "Fu and Li's D" = thetas.plot$fud,
                     "Fay and Wu's H" = thetas.plot$fayh,
                     "Zeng's E" = thetas.plot$zeng
      )
      
      # Output plot of thetas
      plot(thetas.plot$WinCenter,
           data, t="p", pch=19,col=rgb(0,0,0,0.5),
           xlab="Position (bp)",
           ylab=input$selectionChoice,
           main=paste("Neutrality test statistics along chromosome", thetas$Chr[1]), 
           xlim = ranges2$x, ylim = ranges2$y
      )
      if(input$selectionLowess){lines(lowess(thetas.plot$WinCenter,data, f=0.1), col="red")}
    })
    
    # Creating zoom function in selectionPlot1 and selectionPlot2
    observe({
      brush <- input$selectionPlot1_brush
      if (!is.null(brush)) {
        ranges2$x <- c(brush$xmin, brush$xmax)
        ranges2$y <- c(brush$ymin, brush$ymax)
        
      } else {
        ranges2$x <- NULL
        ranges2$y <- NULL
      }
    })

    # Fst Plot
    output$fstChroms = renderUI({
      if(is.null(input$userIntersect)){
        choices <- 12
      }
      else{
        intersect <- dataInputIntersect()
        choices <- unique(intersect$Chr)
      }
      selectInput('fstChrom', 'Chromosome to plot', choices)
    })

    output$fstMin = renderUI({
      if(is.null(input$userIntersect)){
        min <- 10944
        max <- 27530233
      }
      else{
        intersect <- dataInputIntersect()
        min <- min(intersect$bp)
        max <- max(intersect$bp)
      }
      numericInput("intersectLow", "Base Start Position", value=min, min = min, max=max-1)
    })

    output$fstMax = renderUI({
      if(is.null(input$userIntersect)){
        min <- 10944
        max <- 27530233
      }
      else{
        intersect <- dataInputIntersect()
        min <- min(intersect$bp)
        max <- max(intersect$bp)
      }
      numericInput("intersectHigh", "Base Start Position", value=min+10000, min=min+1, max=max)
    })

    output$fstPlot <- renderPlot({
      #error handling code to provide a default dataset to graph
      fst <- tryCatch({
        dataInputFst()
      }, error = function(err) {
        fst <- fread("allo.indica.fst",
                             sep="\t")
      })
      intersect <- tryCatch({
        dataInputIntersect()
      }, error = function(err) {
        intersect <- fread("intersect.allo.indica_intergenic.txt",
                        sep="\t")
      })
      setnames(fst,fst.headers)
      setnames(intersect,intersect.headers)
      fst.intersect <- cbind(intersect,fst)
      fst.intersect <- subset(fst.intersect, Chr==input$fstChrom)
      fst.intersect <- subset(fst.intersect, FST>=0 & FST<=1)

      if(input$annotations){
        validate(need(input$userAnnotations, 'Need GFF file before clicking checkbox!'))
        gff <- gffInput()
        gff.gene <- subset(gff, type="gene")
        gff.df <- data.frame(gff.gene,annotation(gff))
        gff.df.gene <- subset(gff.df, type=="gene")
      }

      if(input$subset) {
        fst.plot <- subset(fst.intersect, bp >= input$intersectLow & bp <= input$intersectHigh)
      }
      else {
        fst.plot <- fst.intersect
      }

      if(input$annotations) {
        plot(fst.plot$bp,
             fst.plot$FST, t="p", pch=19,col=rgb(0,0,0,0.5),
             xlab="Position (bp)",
             ylab="Fst",
             main=paste("Fst along chromosome")
        )
        rug(rect(gff.df.gene$X1, 
                 -1e2, gff.df.gene$X2, 
                 0, 
                 col=rgb(0.18,0.55,0.8,0.75), 
                 border=NA))
        if(input$fstLowess){lines(lowess(fst.plot$bp,
                                         fst.plot$FST, 
                                         f=0.1), 
                                  col="red")}
      }
      else {
        plot(fst.plot$bp,
             fst.plot$FST, 
             t="p", 
             pch=19,
             col=rgb(0,0,0,0.5),
             xlab="Position (bp)",
             ylab=paste("Fst"),
             main=paste("Fst along chromosome", fst.plot$Chr[1])
        )
        if(input$fstLowess){lines(lowess(fst.plot$bp, fst.plot$FST, 
                                         f=0.1), 
                                  col="red")}
      }
    })

    output$SFSPlot <- renderPlot({
      sfs <- tryCatch({
        dataInputSFS()

      },error = function(err) {
        sfs <- scan("sfs_example.txt")
      })
#      subsfs <- sfs[-c(1,length(sfs))]/sum(sfs[-c(1,length(sfs))])
      
# Changed "Chromosome" for x axis label to "Derived Allele Frequency" ### Temporarily comment out Arun's version of the plot    
#      barplot(subsfs, xlab="Derived Allele Frequency", 
#              ylab="Proportion of SNPs", 
#              main="Site Frequency Spectrum",
#              names=1:length(sfs[-c(1,length(sfs))]), 
#              col="#A2C8EC", border=NA)
      sfs.bp <- barplot((alleles/sum(alleles)),
              xaxt = "n",
              xlab = "Derived Allele Frequency",
              ylab = "Proportion of SNPs",
              main = "Site Frequency Spectrum",
              offset = 0,
              xlim = NULL,
              ylim = NULL,
              axes = TRUE,
              names = 1:length(alleles),
              las = 1,
              pch = 18,
              xpd = TRUE,
              col = "magenta1")
      lab <- c(1:length(sfs.bp))
      axis(1, at = sfs.bp, labels = lab)

    })

    # Admixture plot 1
    output$admixPlot1 <- renderPlot({
      admix <- tryCatch({
        dataInputAdmix()
      },error = function(err) {
        admix <- t(as.matrix(read.table("ngsadmix_example.txt")))
      })
      barplot(admix, col=c("#006BA4","#FF800E","#A2C8EC",
                           "#898989","#ABABAB","#595959",
                           "#5F9ED1","#CFCFCF","#FFBC79","#C85200"),
              space=0, 
              border=NA, 
              xlab="Individuals", 
              ylab="admixture proportion")
    })
    
    # Admixture plot 2
    output$admixPlot2 <- renderPlot({
      admix <- tryCatch({
        dataInputAdmix()
        
      },error = function(err) {
        admix <- t(as.matrix(read.table("ngsadmix_example.txt")))
      })
      barplot(admix, col=c("#006BA4","#FF800E","#A2C8EC",
                           "#898989","#ABABAB","#595959",
                           "#5F9ED1","#CFCFCF","#FFBC79","#C85200"),
              space=0, 
              border=NA, 
              xlab="Individuals", 
              ylab="admixture proportion")
    })

    # ABBA BABA plot
    output$ABBABABATree <- renderPlot({
      ABBABABA <- tryCatch({
        dataInputABBABABA()
      }, error = function(err){
        ABBABABA <- read.table("abbababa.test", sep="\t", header=T)
      })
      d.current <- subset(ABBABABA, H2 == input$h2 & H3 == input$h3)
      tree <- read.tree(text=paste("(Outgroup,(", input$h3, ",(", input$h2, ",Taxon)));", sep=""))
      plot(tree, type = "cladogram", edge.width = 2, direction='downwards')

    })
    output$ABBABABAPlot <- renderPlot({
      ABBABABA <- tryCatch({
        dataInputABBABABA()
      }, error = function(err){
        ABBABABA <- read.table("abbababa.test", sep="\t", header=T)
      })
      d.current <- subset(ABBABABA, H2 == input$h2 & H3 == input$h3)
      mypanel.Dotplot <- function(x, y, ...) {
        panel.Dotplot(x,y,...)
        tips <- attr(x, "other")
        panel.abline(v=0, lty=3)
        trellis.par.set(mfrow=c(2,1))
        panel.arrows(x0 = tips[,1], y0 = y,
                     x1 = tips[,2], y1 = y,
                     length = 0.05, unit = "native",
                     angle = 90, code = 3)
      }
      Dotplot(factor(d.current$H1) ~ Cbind(d.current$Dstat, 
                                           d.current$Dstat-d.current$SE, 
                                           d.current$Dstat+d.current$SE), 
              col="blue", pch=20, panel = mypanel.Dotplot,
              xlab="D", ylab="Taxon", 
              title=paste("D statistic comparison where H2=", 
                          input$h2, " and H3=", input$h3, sep=""))
    })

    # Create zoomable PCA plots
    ranges3 <- reactiveValues(x = NULL, y = NULL)
    
    # Create PCAPlot1 you can select areas to zoom in on top
    output$PCAPlot1 <- renderPlot({
      PCA <- tryCatch({
        dataInputPCA()
      }, error = function(err) {
        PCA <- read.table("all.pop.covar", header=F)
      })
      eig <- eigen(PCA, symm=TRUE); # computes eigenvalues of matrices
      eig$val <- eig$val/sum(eig$val);
      PC <- as.data.frame(eig$vectors)
      plot(PC$V1, PC$V2, 
           pch=19, col=rgb(0,0,0,0.4), 
           xlab="PC1", ylab="PC2", 
           main="ngsCovar Results",
           asp=1)
    })
    
    # Create PCAPlot2 that will react to PCAPlot1 selected area
    output$PCAPlot2 <- renderPlot({
      PCA <- tryCatch({
        dataInputPCA()
      }, error = function(err) {
        PCA <- read.table("all.pop.covar", header=F)
      })
      eig <- eigen(PCA, symm=TRUE); # computes eigenvalues of matrices
      eig$val <- eig$val/sum(eig$val);
      PC <- as.data.frame(eig$vectors)
      plot(PC$V1, PC$V2, 
           pch=19, col=rgb(0,0,0,0.4), 
           xlab="PC1", ylab="PC2", 
           main="ngsCovar Results",
           xlim = ranges3$x, ylim = ranges3$y,
           asp=1)
    })
    
    # Creating zoom function in PCAPlot1 and PCAPlot2
    observe({
      brush <- input$PCAPlot1_brush
      if(!is.null(brush)) {
        ranges3$x <- c(brush$xmin, brush$xmax)
        ranges3$y <- c(brush$ymin, brush$ymax)
      } else {
        ranges3$x <- NULL
        ranges3$y <- NULL
      }
    })
  })