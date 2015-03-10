library(shiny)
library(genomeIntervals)
options(shiny.maxRequestSize = -1)
thetas.headers <- c("(indexStart,indexStop)(firstPos_withData,lastPos_withData)(WinStart,WinStop)","Chr","WinCenter","tW","tP","tF","tH","tL","Tajima","fuf","fud","fayh","zeng","nSites")

not.loaded <- TRUE

# Define server logic required to draw a histogram
shinyServer(
  
  
  function(input, output) {
  
  dataInputThetas = reactive({
    data <- input$userThetas
    path <- as.character(data$datapath)
    thetas <- read.table(file=path, 
                         sep="\t", 
                         col.names=thetas.headers
                         )
    return(thetas)
  })
  
  dataInputSFS = reactive({
    data <- input$userSFS
    path <- as.character(data$datapath)
    sfs <- exp(scan(path))
    return(sfs)
    
  })
  
  gffInput = reactive({
    data <- input$userAnnotations
    path <- as.character(data$datapath)
    gff <- readGff3(path)
    
    return(gff)
    
  })
    
  output$thetaPlot <- renderPlot({
    # error handling code to provide a default dataset to graph
    thetas <- tryCatch({
        dataInputThetas()
      }, error = function(err) {
        thetas <- read.table(file="BKN_Diversity.thetas.gz.pestPG", 
                             sep="\t", 
                             col.names=thetas.headers)
      }
      )
    if(input$annotations){
      validate(need(input$userAnnotations, 'Need GFF file before clicking checkbox!'))
      gff <- gffInput()
      gff.gene <- subset(gff, type="gene")
      gff.df <- data.frame(gff.gene,annotation(gff))
      gff.df.gene <- subset(gff.df, type=="gene")
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
    #thetas.plot <- subset(thetas.plot, nsites > 50) 
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
      rug(rect(gff.df.gene$X1, -1e2, gff.df.gene$X2, 0, col=rgb(0.18,0.55,0.8,0.75), border=NA))
      if(input$thetaLowess){lines(lowess(thetas.plot$WinCenter,data, f=0.1), col="red")}
    }
    else {
      plot(thetas.plot$WinCenter, 
           data, t="p", pch=19,col=rgb(0,0,0,0.5), 
           xlab="Position (bp)", 
           ylab=paste(input$thetaChoice,"Estimator Value"), 
           main=paste("Estimators of theta along chromosome", thetas$Chr[1])
      )
      if(input$thetaLowess){lines(lowess(thetas.plot$WinCenter,data, f=0.1), col="red")}
    }
  })
  
  output$selectionPlot <- renderPlot({
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
                   "Zeng's E" = thetas.plot$zengE
    )
    
    plot(thetas.plot$WinCenter, 
         data, t="p", pch=19,col=rgb(0,0,0,0.5), 
         xlab="Position (bp)", 
         ylab=input$selectionChoice, 
         main=paste("Neutrality test statistics along chromosome", thetas$Chr[1])
         )
    if(input$selectionLowess){lines(lowess(thetas.plot$WinCenter,data, f=0.1), col="red")}
  })
  
  output$SFSPlot <- renderPlot({
    sfs <- tryCatch({
      dataInputSFS()
      
    },error = function(err) {
      sfs <- exp(scan("sfs_example.txt"))
    })
    subsfs <- sfs[-c(1,length(sfs))]/sum(sfs[-c(1,length(sfs))])
    barplot(subsfs, xlab="Chromosomes", ylab="Proportion", main="Site Frequency Spectrum",names=1:length(sfs[-c(1,length(sfs))]), col="#A2C8EC", border=NA)
    
  })
})