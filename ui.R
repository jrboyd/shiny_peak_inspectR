library(shiny)
library(mlbench)
library(plotly)
library(shinythemes)
library(dplyr)

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
  
  # Some help text
  h2("Coupled events in plotly charts using Shiny"),
  h4("This Shiny app showcases coupled events using Plotly's ", tags$code("event_data()"), " function."),
  tags$ol(
    tags$li("The first chart showcases", tags$code("plotly_selected")),
    tags$li("The third chart showcases", tags$code("plotly_click"))
  ),
  
  # Vertical space
  tags$hr(),
  
  # Feature selection
  fixedRow(
    # column(3, selectInput(inputId = "featureInput1", label = "Select first feature", choices = featureList, selected = "Cell.Size")),
    column(6, selectInput(inputId = "featureInput2", label = "Select gene type", 
                          choices = c("all", "protein_coding", "non_coding"), selected = "all")),
    column(6, selectInput(inputId = "featureInput3", label = "Select alignment type", 
                          choices = c("peak summit in 2kb", "tss"), selected = "peak summit in 2kb"))
  ),
  
  # First row
  fixedRow(
    column(6, plotlyOutput("Plot1", height = "600px")),
    column(6, plotlyOutput("Plot2", height = "600px"))),
  
  tags$hr(),
  tags$blockquote("First drag a selection box in the scatter plot to populate the barchart. Then select one of the bars in the barchat
    to populate the boxplot"),
  
  
  # Second row
  fixedRow(
    column(6, uiOutput("List1")),
    column(6, plotlyOutput("Plot3", height = "600px"))
  )
)