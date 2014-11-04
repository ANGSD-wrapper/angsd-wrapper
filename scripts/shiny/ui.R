library(shiny)

shinyUI(fluidPage(
  
  # Application title
  titlePanel("angsd-wrapper graph"),

  sidebarLayout(
    sidebarPanel(
        fileInput('userThetas',
                  label= 'Choose thetas File'
        ),
        
        selectInput("thetaChoice",
                    label = "Choose estimator of theta to graph", 
                    choices = c("Watterson's Theta", "Pairwise Theta", "Fu and Li's Theta", "Fay's Theta", "Maximum likelihood (L) Theta"),
                    selected = "Watterson's Theta"
        ),
        selectInput("selectionChoice",
                    label = "Choose a neutrality test statistic to graph", 
                    choices = c("Tajima's D", "Fi and Li's D", "Fu and Li's F", "Fay and Wu's H", "Zeng's E"),
                    selected = "Tajima's D"
        ),
        hr(),
        checkboxInput("subset","Subset data?", value=FALSE),
        numericInput("WinCenterLow", "Base Start Position", value=0),
        numericInput("WinCenterHigh", "Base End Position", value=10000),
        hr(),
        checkboxInput("annotations","GFF annotations?", value=FALSE),
        fileInput('userAnnotations',
                  label= 'Choose GFF File'
        )
      ),
    # Show a plot of the thetas
      mainPanel(
        plotOutput("thetaPlot"),
        plotOutput("selectionPlot")
      )
    )
  )
)