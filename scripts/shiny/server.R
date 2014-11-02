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
    data <- switch(input$thetaChoice,
                   "tW" = thetas$tW,
                   "tP" = thetas$tP,
                   "tH" = thetas$tH,
                   "tL" = thetas$tL
                   )
    plot(thetas$WinCenter, data, t="l", xlab="Position (bp)", ylab=paste(input$thetaChoice,"Estimator Value"), main=paste("Estimators of theta along chromosome", thetas$Chr[1]))
#    ggplot(thetas, aes(WinCenter)) + 
#      geom_line(aes(y=data, colour="Theta Estimate", alpha=0.5))
  })
})