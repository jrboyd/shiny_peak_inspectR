ui_tab_analyze = function(){
  tagList(h3("scatterplot of bigwig values"), 
          uiOutput("ProfilesLoaded"),
          uiOutput("GroupingsAvailable"),
          # plotlyOutput("xy_values"),
          fixedRow(
            column(width = 3,
                   radioButtons("RadioScatterplotType", label = "Plot Type", choices = c("volcano", "standard"), selected = "volcano"),
                   checkboxInput("CheckShowHelpers", label = "Show Plot Help", value = T)
            ),
            column(width = 9,
                   plotOutput("xy_values", 
                              click = "xy_click",
                              brush = brushOpts(
                                id = "xy_brush"
                              ))
            )
          ),
          h3("aggregated summary of selected values"),
          fixedRow(
            column(width = 3,
                   uiOutput("AggregateDisplayed")
            ),
            column(width = 9,
                   plotOutput("PlotAggregated")
            )
          ),
          h3("heatmap of all profiles"),
          h3("randomly selected individual profiles"),
          fixedRow(
            column(width = 3,
                   uiOutput("SelectIndividualRegionsDisplayed")
            ),
            column(width = 9,
                   plotOutput("PlotIndividualRegionsDisplayed", width = 360, height = 600)
            )
          ),
          DT::dataTableOutput("XY_Selected"))
  
}
