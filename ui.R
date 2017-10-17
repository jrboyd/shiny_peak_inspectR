library(shiny)
library(mlbench)
library(plotly)
library(shinythemes)
library(dplyr)
library(shinyFiles)
source("tab_setup.R")
source("tab_process.R")
source("tab_analyze.R")
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
              tabPanel(title = "1) setup", value = "1", ui_tab_setup()),
              tabPanel(title = "2) process", value = "2", ui_tab_process()),
              tabPanel(title = "3) analyze", value = "3", ui_tab_analyze())
  )
  

)