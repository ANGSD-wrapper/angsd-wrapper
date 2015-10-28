library(shiny)

shinyUI(fluidPage(
  # Application title
  titlePanel("ANGSD-wrapper graph"),
  tabsetPanel(
    tabPanel(
      "Thetas",
      sidebarLayout(
        sidebarPanel(
          fileInput('userThetas',
                    label= "Choose '.pestPG' Thetas File"
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
          uiOutput('thetaChroms'),
          hr(),
          checkboxInput("thetaLowess","Theta Lowess", value=FALSE),
          checkboxInput("selectionLowess","Neutrality Test Statistic Lowess", value=FALSE),
          hr(),
          numericInput("WinCenterLow", "Base Start Position", value=0),
          numericInput("WinCenterHigh", "Base End Position", value=10000),
          checkboxInput("subset","Toggle subset data", value=FALSE),
          hr(),
          checkboxInput("rm.nsites", "Remove data where nSites < x", value=FALSE),
          numericInput("nsites", "x:",value=0),
          
          hr(),
          fileInput('userAnnotations',
                    label= 'Choose GFF File'
          ),
          checkboxInput("annotations","Toggle GFF annotations", value=FALSE)
        ),
        # Show a plot of the thetas
        mainPanel(
          plotOutput("thetaPlot"), 
          plotOutput("selectionPlot")
        )
      )
    ),
    tabPanel(
      "SFS",
      sidebarLayout(
        sidebarPanel(
          fileInput('userSFS',
                    label= "Choose '_DerivedSFS' SFS File"
          )
          
        ),
        mainPanel(
          plotOutput("SFSPlot")
        )
      )
    ),
    tabPanel(
      "ABBA BABA",
      sidebarLayout(
        sidebarPanel(
          fileInput('userABBABABA',
                    label= "Choose 'abbababa.txt' ABBABABA File"
          ),
          textInput('h2', label="H2", value='NA11993'),
          textInput('h3', label="H3", value='NA12763')
          
        ),
        mainPanel(
          plotOutput("ABBABABATree"),
          plotOutput("ABBABABAPlot"),
          helpText(a("https://youtu.be/unfzfe8f9NI", href="https://youtu.be/unfzfe8f9NI"))
        )
      )
    ),
    tabPanel(
      "Fst",
      sidebarLayout(
        sidebarPanel(
          fileInput('userFst',
                    label= "Choose '.fst' Fst File"
          ),
          fileInput('userIntersect',
                    label= 'Choose Intersect File'
          ),
          
          uiOutput('fstChroms'),
          
          hr(),
          checkboxInput("fstLowess","Fst Lowess", value=FALSE),
          hr(),
          
          uiOutput('fstMin'),
          uiOutput('fstMax'),
          
          checkboxInput("subset","Toggle subset data", value=FALSE),
          
          hr(),
          fileInput('userAnnotations',
                    label= 'Choose GFF File'
          ),
          checkboxInput("annotations","Toggle GFF annotations", value=FALSE)
          
        ),
        mainPanel(
          plotOutput("fstPlot")
        )
      )
    ),
    tabPanel(
      "PCA",
      sidebarLayout(
        sidebarPanel(
          fileInput('userPCA',
                    label= "Choose a '_PCA.covar' File"
          )
        ),
        mainPanel(
          plotOutput("PCAPlot"), 
          width=4
        )
      )
    ),
    tabPanel(
      "Admixture",
      sidebarLayout(
        sidebarPanel(
          fileInput('userAdmix',
                    label= "Choose '.qopt' admixture File"),
          numericInput("k", "Number of K ancestral populations to graph:", 1,
                       min = 1, max = 10)
          
        ),
        mainPanel(
          plotOutput("admixPlot")
          
        )
      )
    )
    #end of tabsetPanel
  )
)
)