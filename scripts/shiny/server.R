library(shiny)
library(genomeIntervals)
options(shiny.maxRequestSize = -1)
thetas.headers <- c("(indexStart,indexStop)(firstPos_withData,lastPos_withData)(WinStart,WinStop)","Chr","WinCenter","tW","tP","tF","tH","tL","Tajima","fuf","fud","fayh","zeng","nSites")

not.loaded <- TRUE
# thetas <- read.table(file="BKN_Diversity.thetas.gz.pestPG",
#                      sep="\t", 
#                      col.names=thetas.headers
#                      )
#gff <- readGff3("Zea_mays.AGPv3.23.chromosome.10.gff3.gz")

# Define server logic required to draw a histogram
shinyServer(
  function(input, output) {
  
  dataInput = reactive({
    data <- input$userThetas
    path <- as.character(data$datapath)
    thetas <- read.table(file=path, 
                         sep="\t", 
                         col.names=thetas.headers
                         )
    return(thetas)
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
        dataInput()
      }, error = function(err) {
        thetas <- read.table(file="BKN_Diversity.thetas.gz.pestPG", 
                             sep="\t", 
                             col.names=thetas.headers
                             )
      }
      )
    if(input$annotations){
      validate(need(input$userAnnotations, 'Need GFF file before clicking checkbox!'))
      gff <- gffInput()
#       tryCatch({
#         gffInput()
#       }, error = function(err) {
#         #gff <- readGff3("Zea_mays.AGPv3.23.chromosome.10.gff3.gz")
#         print("Need to provide a GFF file!")
#       }
#       )
      gff.gene <- subset(gff, type="gene")
    }

    if(input$subset) {
      thetas.plot <- subset(thetas, WinCenter > input$WinCenterLow & WinCenter < input$WinCenterHigh)
      }
      else {
        thetas.plot <- thetas
      }
    
#     if(input$annotations) {
#       gff.plot <- subset(gff.gene, gff.gene[,1] > input$WinCenterLow & gff.gene[,1] < input$WinCenterHigh)
#     }
    
    data <- switch(input$thetaChoice,
                   "Watterson's Theta" = thetas.plot$tW,
                   "Pairwise Theta" = thetas.plot$tP,
                   "Fu and Li's Theta" = thetas.plot$tF,
                   "Fay's Theta" = thetas.plot$tH,
                   "Maximum likelihood (L) Theta" = thetas.plot$tL
    )
    if(input$annotations) {
      plot(thetas.plot$WinCenter, 
           data, t="l", 
           xlab="Position (bp)", 
           ylab=paste(input$thetaChoice,"Estimator Value"), 
           main=paste("Estimators of theta along chromosome", thetas$Chr[1])#,
           #panel.first=rect(gff.gene[,1], -1e6, gff.gene[,2], 1e6, col=rgb(0,1,0,0.1), border=NA)
           )
      rug(rect(gff.gene[,1], -1e2, gff.gene[,2], 0, col=rgb(0,1,0,0.2), border=NA))
    }
    else {
      plot(thetas.plot$WinCenter, 
           data, t="l", 
           xlab="Position (bp)", 
           ylab=paste(input$thetaChoice,"Estimator Value"), 
           main=paste("Estimators of theta along chromosome", thetas$Chr[1])
      )
    }
  })
  
  output$selectionPlot <- renderPlot({
    # error handling code to provide a default dataset to graph
    thetas <- tryCatch({
        dataInput()
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
    
    data <- switch(input$selectionChoice,
                   "Tajima's D" = thetas.plot$Tajima,
                   "Fu and Li's F" = thetas.plot$fuf,
                   "Fu and Li's D" = thetas.plot$fud,
                   "Fay and Wu's H" = thetas.plot$fayh,
                   "Zeng's E" = thetas.plot$zengE
    )
    
    plot(thetas.plot$WinCenter, 
         data, t="l", 
         xlab="Position (bp)", 
         ylab=input$selectionChoice, 
         main=paste("Neutrality test statistics along chromosome", thetas$Chr[1])
         )
  })
})