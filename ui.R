library(shiny)
library(mlbench)
library(plotly)
library(shinythemes)
library(dplyr)
library(shinyFiles)
source("tab1_setup.R")
source("tab2_analyze.R")
# Load data
data(BreastCancer)

# Remove NAs
BreastCancer <- na.omit(BreastCancer)

# Remvove ID
BreastCancer <- BreastCancer[,-1]

# Store features and actual class in seprate variables
featureList <- colnames(BreastCancer)[-10]
class <- BreastCancer$Class

# Convert to numeric
BreastCancer[,1:9] <- apply(BreastCancer[,-10], 2, as.numeric)

# ui.R definition
ui <- fluidPage(
  # Set theme
  theme = shinytheme("spacelab"),
  tabsetPanel(id = "tabset",
              tabPanel(title = "setup", ui_tab1_setup()),
              tabPanel(title = "analyze", ui_tab2_analyze())
  ),
  
  # Some help text
  h2("Coupled events in plotly charts using Shiny"),
  h4("This Shiny app showcases coupled events using Plotly's ", tags$code("event_data()"), " function."),
  tags$ol(
    tags$li("The first chart showcases", tags$code("plotly_selected")),
    tags$li("The third chart showcases", tags$code("plotly_click"))
  ),
  
  # Vertical space
  tags$hr(),
  
  # # Feature selection
  # fixedRow(
  #   # column(3, selectInput(inputId = "featureInput1", label = "Select first feature", choices = featureList, selected = "Cell.Size")),
  #   column(6, selectInput(inputId = "featureInput2", label = "Select gene type", 
  #                         choices = c("all", "protein_coding", "non_coding"), selected = "all")),
  #   column(6, selectInput(inputId = "featureInput3", label = "Select alignment type", 
  #                         choices = c("peak summit in 2kb", "tss"), selected = "peak summit in 2kb"))
  # ),
  # 
  # # First row
  # fixedRow(
  #   column(6, plotlyOutput("PlotScatter2d", height = "600px")),
  #   column(6, plotlyOutput("PlotAggregateProfile", height = "600px"))),
  # 
  # tags$hr(),
  # tags$blockquote("First drag a selection box in the scatter plot to populate the barchart. Then select one of the bars in the barchat
  #   to populate the boxplot"),
  # 
  # 
  # # Second row
  # fixedRow(
  #   column(6, uiOutput("List1")),
  #   column(6, plotlyOutput("PlotSelectedProfiles", height = "600px"))
  # ),
  plotOutput("PlotTmp"),
  plotlyOutput("xy_values"),
  shinyFilesButton(id = "FilesLoadSet", label = "Find Files on Server", title = "Find Peaks to Annotate", multiple = F),
  fileInput(inputId = "UploadLoadSet", label = "Browse Local Files")
)