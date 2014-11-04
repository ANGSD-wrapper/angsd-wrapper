library(shiny)
library(ggplot2)
thetas <- read.table(file="BKN_Diversity.thetas.gz.pestPG", sep="\t", col.names=c("(indexStart,indexStop)(firstPos_withData,lastPos_withData)(WinStart,WinStop)","Chr","WinCenter","tW","tP","tF","tH","tL","Tajima","fuf","fud","fayh","zeng","nSites"))

# Define server logic required to draw a histogram
shinyServer(
  function(input, output) {
  
  dataInput = reactive({
    data <- input$userThetas
    path <- as.character(data$datapath)
    thetas <- read.table(file=path, sep="\t", col.names=c("(indexStart,indexStop)(firstPos_withData,lastPos_withData)(WinStart,WinStop)","Chr","WinCenter","tW","tP","tF","tH","tL","Tajima","fuf","fud","fayh","zeng","nSites"))
    return(thetas)
  })
    
  output$thetaPlot <- renderPlot({
    # error handling code to provide a default dataset to graph
    thetas <- tryCatch({
      dataInput()
    }, error = function(err) {
      thetas <- read.table(file="BKN_Diversity.thetas.gz.pestPG", sep="\t", col.names=c("(indexStart,indexStop)(firstPos_withData,lastPos_withData)(WinStart,WinStop)","Chr","WinCenter","tW","tP","tF","tH","tL","Tajima","fuf","fud","fayh","zeng","nSites"))
    }
      )

    if(input$subset) {
      thetas.plot <- subset(thetas, WinCenter > input$WinCenterLow & WinCenter < input$WinCenterHigh)
    }
    else {
      thetas.plot <- thetas
    }
    
    data <- switch(input$thetaChoice,
                   "Watterson's Theta" = thetas.plot$tW,
                   "Pairwise Theta" = thetas.plot$tP,
                   "Fu and Li's Theta" = thetas.plot$tF,
                   "Fay's Theta" = thetas.plot$tH,
                   "Maximum likelihood (L) Theta" = thetas.plot$tL
    )
    plot(thetas.plot$WinCenter, data, t="l", xlab="Position (bp)", ylab=paste(input$thetaChoice,"Estimator Value"), main=paste("Estimators of theta along chromosome", thetas$Chr[1]))
  })
  
  output$selectionPlot <- renderPlot({
    # error handling code to provide a default dataset to graph
    thetas <- tryCatch({
      dataInput()
    }, error = function(err) {
      thetas <- read.table(file="BKN_Diversity.thetas.gz.pestPG", sep="\t", col.names=c("(indexStart,indexStop)(firstPos_withData,lastPos_withData)(WinStart,WinStop)","Chr","WinCenter","tW","tP","tF","tH","tL","Tajima","fuf","fud","fayh","zeng","nSites"))
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
    
    plot(thetas.plot$WinCenter, data, t="l", xlab="Position (bp)", ylab=input$selectionChoice, main=paste("Neutrality test statistics along chromosome", thetas$Chr[1]))
  })
})