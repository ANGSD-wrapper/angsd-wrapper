library(shiny)
library(shinythemes)

shinyUI(fluidPage(
  # Theme
  theme = shinytheme("cerulean"),
  tags$style("body {background-color: #FFFFFF; }"),
  # Application title
  titlePanel("ANGSD-wrapper graph"),
  tabsetPanel(

    # Tab 1 - Thetas
    tabPanel("Thetas",
      sidebarLayout(
        sidebarPanel(
          headerPanel("Thetas Graphs", windowTitle = "Thetas Graphs"),
          # Choose input file
          fileInput('userThetas',
                    label= "Choose '.pestPG' Thetas File"
          ),

          # Watterson's theta
          selectInput("thetaChoice",
                      # Need to specify file to choose
                      label = "Choose estimator of theta to graph", 
                      choices = c("Watterson's Theta",
                                  "Pairwise Theta",
                                  "Fu and Li's Theta",
                                  "Fay's Theta",
                                  "Maximum likelihood (L) Theta"),
                      selected = "Watterson's Theta"
          ),
          
          # Tajima's D
          selectInput("selectionChoice",
                      label = "Choose a neutrality test statistic to graph", 
                      choices = c("Tajima's D",
                                  "Fi and Li's D",
                                  "Fu and Li's F",
                                  "Fay and Wu's H",
                                  "Zeng's E"),
                      selected = "Tajima's D"
          ),
          uiOutput('thetaChroms'),
          # Adding additional "help" text
          tags$h4(class ="header",
                  tags$h5("Samples can be removed from graph by using cursor in text box and using backspace."),
                  tags$h5("You can also type parts of a sample name to search for it.")),
          hr(),
          checkboxInput("thetaLowess",
                        "Theta Lowess", 
                        value=FALSE),
          checkboxInput("selectionLowess",
                        "Neutrality Test Statistic Lowess", 
                        value=FALSE),
          hr(),
          tags$h4(class = "header",
                  tags$h4("Interactive Zoom and Hover:"),
                  hr(),
                  tags$h5("There are two ways you can zoom in on the plots:"),
                  tags$h6("1. Click and drag on the upper graph to select an area to zoom in on."),
                  tags$h6("Lower graph will show zoomed area selected. If you hover over a point on the lower graph, a table of the 'x' and 'y' values will be shown."),
                  p("or"),
                  tags$h6("2. Click the 'Toggle subset data' checkbox and enter in the interval you want to zoom in on."),
                  hr()),
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
              h5("Click and drag to select area to zoom on this plot"),
              plotOutput("thetaPlot1",
                         dblclick = "thetaPlot1_dblclick",
                         brush = brushOpts("thetaPlot1_brush",
                                           resetOnNew = TRUE))
              ),
            column(
              width = 12, height = 300,
              h5("Plot will zoom in on area selected above"),
              plotOutput("thetaPlot2",
                         hover = hoverOpts("thetaPlot2_hover",
                                           delay = 50,
                                           nullOutside = FALSE)
            ),
            column(
              width = 6,
              verbatimTextOutput("thetaPlot2_hoverinfo")
            ),
            column(
              width = 12, height = 300,
              h5("Click and drag to select area to zoom on this plot"), 
              plotOutput("selectionPlot1",
                       dblclick = "selectionPlot1_dblclick",
                       brush = brushOpts(id = "selectionPlot1_brush",
                                         resetOnNew = TRUE))
            ),
            column(
              width = 12, height = 300,
              h5("Plot will zoom in on area selected above"),
              plotOutput("selectionPlot2",
                         hover = hoverOpts("selectionPlot2_hover",
                                           delay = 50,
                                           nullOutside = FALSE))
            ),
            column(
              width = 6,
              verbatimTextOutput("selectionPlot2_hoverinfo")
            )
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
          headerPanel("ABBA BABA",
                      windowTitle = "ABBA BABA Graph"),
          fileInput('userABBABABA',
                    label= "Choose 'abbababa.txt' ABBABABA File"
          ),
          tags$h5(class = "header",
                  tags$h5("Change the orientation of individuals/populations for your tree by typing name of H2 and H3 groups"),
                  tags$b("Note: the names of H2 and H3 must be an exact match with the names outputted in the table."),
                  hr()),
          textInput('h2', label="H2", value='NA11993'),
          textInput('h3', label="H3", value='NA12763'),
          hr(),
          tags$h4(class = "header",
                  tags$h5("The table includes the output ABBABABA file from ANGSD-wrapper with the addition of a 'Pvalue' column."),
                  tags$h5("1. Output of ABBABABA file from ANGSD-wrapper"),
                  tags$h5("2. Additional 'Pvalue' column calculated with the formula 2*pnorm(-abs(z-score))"),
                  hr(),
                  tags$h5("The table should be used as supplemental data to determine the most significant tree orientation. The table is searchable and can be organized by column."))
        ),
        mainPanel(
          plotOutput("ABBABABATree"),
          plotOutput("ABBABABAPlot"),
          dataTableOutput("ABBABABATable"),
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
          headerPanel("Fst",
                      windowTitle = "Fst Graph"),
          fileInput('userFst',
                    label= "Choose '.fst.graph.me' Fst File"
          ),
          tags$h4(class = "header",
                  tags$h4("Interactive Zoom and Hover:"),
                  hr(),
                  tags$h5("Click and drag on the upper graph to select an area to zoom in on."),
                  tags$h5("Lower graph will show zoomed area selected. If you hover over a point on the lower graph, a table of the 'x' and 'y' values will be shown.")),
          uiOutput('fstChroms'),
          hr(),
          checkboxInput("fstLowess",
                        "Fst Lowess",
                        value=FALSE),
          hr(),
          uiOutput('fstMin'),
          uiOutput('fstMax'),
          
          hr(),
          fileInput('userAnnotations',
                    label= "Choose '.gff' GFF File"
          ),
          checkboxInput("annotations",
                        "Toggle GFF annotations",
                        value=FALSE)

        ),
        mainPanel(
          fluidRow(
            column(
              width = 12, height = 300,
              h5("Click and drag to select area to zoom on this plot"),
              plotOutput("fstPlot1",
                         dblclick = "fstPlot1_dblclick",
                         brush = brushOpts("fstPlot1_brush",
                                           resetOnNew = TRUE))
            ),
            column(
              width = 12, height = 300,
              h5("Plot will zoom in on area selected above"),
              plotOutput("fstPlot2",
                         hover = hoverOpts("fstPlot2_hover",
                                           delay = 50,
                                           nullOutside = FALSE))
            ),
            column(
              width = 6, verbatimTextOutput("fstPlot2_hoverinfo")
            )
          )
        )
      )
    ),
    
    # Tab 5 - PCA
    tabPanel(
      "PCA",
      sidebarLayout(
        sidebarPanel(
          headerPanel("Principal Component Analysis",
                      windowTitle = "PCA Graph"),
          fileInput('userPCA',
                    label= "Choose a '_PCA.covar' File"
          ),
          tags$h4(class = "header",
                  tags$h4("Interactive Zoom and Hover:"),
                  hr(),
                  tags$h5("Click and drag on the upper graph to select an area to zoom in on."),
                  tags$h5("Lower graph will show zoomed area selected. If you hover over a point on the lower graph, a table of the 'x' and 'y' values will be shown."))
        ),
        mainPanel(
          fluidRow(
            column(
              width = 12, height = 300,
              h5("Click and drag to select area on this plot"),
              plotOutput("PCAPlot1",
                         dblclick = "PCAPlot1_dblclick",
                         brush = brushOpts("PCAPlot1_brush",
                                           delay = 50,
                                           resetOnNew = TRUE))
            ),
            column(
              width = 12, height = 300,
              h5("Plot will zoom in on area selected above",
                 plotOutput("PCAPlot2",
                            hover = hoverOpts("PCAPlot2_hover",
                                              nullOutside = FALSE)))
            ),
            column(
              width = 6,
              verbatimTextOutput("PCAPlot2_hoverinfo")
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
          fileInput('userAdmix',
                    label= "Choose '.qopt' admixture File",
                    multiple = TRUE)#,
#          checkboxGroupInput("checkGroup",
#                             label = )
        ),
        mainPanel(
          fluidRow(
            plotOutput("admixPlot")
          )
          )
        )
      )
    )
    #end of tabsetPanel
  )
)