library(shiny)

shinyUI(fluidPage(
  
  # Application title
  titlePanel("angsd-wrapper theta graph"),

  sidebarLayout(
    sidebarPanel(
#       checkboxGroupInput("thetaChoice", 
#                          label = h3("Choose estimators of theta to graph"), 
#                          choices = list("tW" = 1, 
#                                         "tP" = 2,
#                                         "tH" = 3,
#                                         "tL" = 4),
#                          selected = 1)
      
      selectInput("thetaChoice",
                  label = h6("Choose estimator of theta to graph"), 
                  choices = c("tW", "tP", "tH", "tL"),
                  selected = "tW"
                  ),
      fileInput('userThetas',
                label= h6('Choose thetas File')
      )
      ),
    
    # Show a plot of the thetas
      mainPanel(
        plotOutput("thetaPlot")
      )
    )
  )
)