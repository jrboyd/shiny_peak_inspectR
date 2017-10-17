library(shiny)
library(shinyjs)
library(mlbench)
library(plotly)
library(shinythemes)
library(dplyr)
library(shinyFiles)
source("tab_setup.R")
source("tab_process.R")
source("tab_analyze.R")
source("shinyjs_tab.R")

# ui.R definition
ui <- fluidPage(
  #Activate shinyjs
  useShinyjs(), 
  #Exgend shinjs to support tabsetsPanel operations
  extendShinyjs(text = tab_jscode),
  inlineCSS(tab_css),
  # Set theme
  theme = shinytheme("spacelab"),
  tabsetPanel(id = "navbar",
              tabPanel(title = "1) setup", 
                       value = "tab1",
                       # value = "1", 
                       ui_tab_setup()),
              tabPanel(title = "2) process", 
                       value = "tab2",
                       # value = "2", 
                       ui_tab_process()),
              tabPanel(title = "3) analyze", 
                       value = "tab3",
                       # value = "3", 
                       ui_tab_analyze())
  )
)