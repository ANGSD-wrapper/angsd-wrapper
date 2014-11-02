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
      

        fileInput('userThetas',
                  label= 'Choose thetas File'
        ),
        
        selectInput("thetaChoice",
                    label = "Choose estimator of theta to graph", 
                    choices = c("tW", "tP", "tH", "tL"),
                    selected = "tW"
        ),
        checkboxInput("subset","Subset data?", value=FALSE),
        numericInput("WinCenterLow", "Base Start Position", value=0),
        numericInput("WinCenterHigh", "Base End Position", value=10000)
      ),
    
    # Show a plot of the thetas
      mainPanel(
        plotOutput("thetaPlot")
      )
    )
  )
)