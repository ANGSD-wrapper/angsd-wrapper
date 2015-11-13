library(shiny)

shinyUI(fluidPage(
  # Application title
  titlePanel("ANGSD-wrapper graph"),
  tabsetPanel(
    
    # Tab 1 - Thetas
    tabPanel("Thetas",
      sidebarLayout(
        sidebarPanel(
          headerPanel("Thetas Graphs", 
                      windowTitle = "Thetas Graphs"),
          # Choose input file
          fileInput('userThetas',
                    label= "Choose '.pestPG' Thetas File"
          ),
          
          # Watterson's theta
          selectInput("thetaChoice",
                      # Need to specify file to choose
                      label = "Choose estimator of theta to graph", 
                      choices = c("Watterson's Theta", "Pairwise Theta", "Fu and Li's Theta", "Fay's Theta", "Maximum likelihood (L) Theta"),
                      selected = "Watterson's Theta"
          ),
          
          # Tajima's D
          selectInput("selectionChoice",
                      label = "Choose a neutrality test statistic to graph", 
                      choices = c("Tajima's D", "Fi and Li's D", "Fu and Li's F", "Fay and Wu's H", "Zeng's E"),
                      selected = "Tajima's D"
          ),
          uiOutput('thetaChroms'),
          # Adding additional 'help' text
          tags$h4(class ="header",
                  tags$h4("Samples can be removed from graph by using cursor in text box and using backspace."),
                  tags$h4("You can also type parts of a sample name to search for it.")),
          hr(),
          checkboxInput("thetaLowess",
                        "Theta Lowess", 
                        value=FALSE),
          checkboxInput("selectionLowess",
                        "Neutrality Test Statistic Lowess", 
                        value=FALSE),
          hr(),
          tags$h4(class = "header",
                  tags$h4("There are two ways you can zoom in on the plots"),
                  tags$h4("1. You can click and drag to select area over the top graph and have the graph zoom in on the graph below"),
                  tags$h4("or"),
                  tags$h4("2. You can click the 'Toggle subset data' checkbox and enter in the interval you want to zoom in on")
                  ),
          numericInput("WinCenterLow", 
                       "Base Start Position", 
                       value=0),
          numericInput("WinCenterHigh", 
                       "Base End Position", 
                       value=10000),
          checkboxInput("subset",
                        "Toggle subset data", 
                        value=FALSE),
          hr(),
          checkboxInput("rm.nsites", 
                        "Remove data where nSites < x", 
                        value=FALSE),
          numericInput("nsites", 
                       "x:",
                       value=0),
          hr(),
          fileInput('userAnnotations',
                    label= "Choose '.gff' File"
          ),
          checkboxInput("annotations",
                        "Toggle GFF annotations", 
                        value=FALSE)
        ),
        # Show a plot of the thetas
        mainPanel(
          fluidRow(
            column(
              width = 12, height = 300,
              h3("Click and drag to select area to zoom on this plot"),
              plotOutput("thetaPlot1", 
                         dblclick = "thetaPlot1_dblclick", 
                         brush = brushOpts(id = "thetaPlot1_brush", 
                                           resetOnNew = TRUE))
              ),
            column(
              width = 12, height = 300,
              h3("Plot will zoom in on area selected above"),
              plotOutput("thetaPlot2")
            ),
            column(
              width = 3,
              verbatimTextOutput("thetaPlot2_clickinfo")
            ),
            column(
              width = 3,
              verbatimTextOutput("thetaPlot2_hoverinfo")
            ),
            column(
              width = 12, height = 300,
              h3("Click and drag to select area to zoom on this plot"), 
              plotOutput("selectionPlot1",
                       dblclick = "selectionPlot1_dblclick",
                       brush = brushOpts(id = "selectionPlot1_brush",
                                         resetOnNew = TRUE))
            ),
            column(
              width = 12, height = 300,
              h3("Plot will zoom in on area selected above"),
              plotOutput("selectionPlot2")
            )
          )
      )
      )
    ),
    
    # Tab 2 - SFS
    tabPanel(
      "SFS",
      sidebarLayout(
        sidebarPanel(
          headerPanel("Site Frequency Spectrum", 
                      windowTitle = "SFS Graph"),
          fileInput('userSFS',
                    label= "Choose '_DerivedSFS' SFS File"
          )
          
        ),
        mainPanel(
          plotOutput("SFSPlot")
        )
      ) 
    ),
    
    # Tab 3 - ABBA BABA
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
          helpText(a("https://youtu.be/unfzfe8f9NI", 
                     href="https://youtu.be/unfzfe8f9NI"))
        )
      )
    ),
    
    # Tab 4 - Fst
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
          checkboxInput("fstLowess",
                        "Fst Lowess", 
                        value=FALSE),
          hr(),
          
          uiOutput('fstMin'),
          uiOutput('fstMax'),
          
          checkboxInput("subset",
                        "Toggle subset data", 
                        value=FALSE),
          
          hr(),
          fileInput('userAnnotations',
                    label= "Choose '.gff' GFF File"
          ),
          checkboxInput("annotations",
                        "Toggle GFF annotations", 
                        value=FALSE)
          
        ),
        mainPanel(
          plotOutput("fstPlot")
        )
      )
    ),
    
    # Tab 5 - PCA
    tabPanel(
      "PCA",
      sidebarLayout(
        sidebarPanel(
          fileInput('userPCA',
                    label= "Choose a '_PCA.covar' File"
          ),
          tags$h5(class ="header",
                  tags$p("The zoom function works as follows:"),
                  tags$p("Click and drag over area you want to select on the top graph, the graph will zoom in on the bottom graph"))
        ),
        mainPanel(
          fluidRow(
            column(
              width = 12, height = 300,
              h3("Click and drag to select area on this plot"),
              plotOutput("PCAPlot1",
                         dblclick = "PCAPlot1_dblclick",
                         brush = brushOpts(id = "PCAPlot1_brush",
                                           resetOnNew = TRUE))
            ),
            column(
              width = 12, height = 300,
              h3("Plot will zoom in on area selected above",
                 plotOutput("PCAPlot2"))
            )
          )
        )
      )
    ),
    
    # Tab 6 - Admixture
    tabPanel(
      "Admixture",
      sidebarLayout(
        sidebarPanel(
          headerPanel('Admixture plots',
                      windowTitle = 'Admixture'),
          fileInput('userAdmix1',
                    label= "Choose '.qopt' admixture File",
                    multiple = TRUE),
          fileInput('userAdmix2',
                    label= "Choose '.qopt' admixture File")
#          numericInput("k", 
#                       "Number of K ancestral populations to graph:",
#                       1,
#                       min = 1, 
#                       max = 10)
        ),
        mainPanel(
          plotOutput("admixPlot1"),
          plotOutput("admixPlot2")
#          column(
#            width = 12, height = 200,
#            plotOutput("admixPlot1")
#          ),
#          column(
#            width = 12, height = 200,
#            plotOutput("admixPlot2")
          )
        )
      )
    )
    #end of tabsetPanel
  )
)