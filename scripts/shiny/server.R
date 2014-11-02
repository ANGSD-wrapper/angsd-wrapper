library(shiny)
library(ggplot2)
thetas <- read.table(file="~/Documents/Science/angsd-wrapper-shiny/BKN_Diversity.thetas.gz.pestPG", sep="\t", col.names=c("(indexStart,indexStop)(firstPos_withData,lastPos_withData)(WinStart,WinStop)","Chr","WinCenter","tW","tP","tF","tH","tL","Tajima","fuf","fud","fayh","zeng","nSites"))

# Define server logic required to draw a histogram
shinyServer(
  function(input, output) {
  
  # Expression that generates a histogram. The expression is
  # wrapped in a call to renderPlot to indicate that:
  #
  #  1) It is "reactive" and therefore should re-execute automatically
  #     when inputs change
  #  2) Its output type is a plot
  
  output$thetaPlot <- renderPlot({
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