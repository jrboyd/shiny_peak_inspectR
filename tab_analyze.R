ui_tab_analyze = function(){
  tagList(h3("scatterplot of bigwig values"), 
          uiOutput("ProfilesLoaded"),
          uiOutput("GroupingsAvailable"),
          # plotlyOutput("xy_values"),
          plotOutput("xy_values", 
                     click = "xy_click",
                     brush = brushOpts(
                       id = "xy_brush"
                     )),
          h3("aggregated summary of selected values"),
          uiOutput("AggregateDisplayed"),
          plotOutput("PlotAggregated"),
          h3("heatmap of all profiles"),
          h3("randomly selected individual profiles"),
          sidebarLayout(
            sidebarPanel = 
              uiOutput("SelectIndividualRegionsDisplayed"),
            mainPanel = 
              plotOutput("PlotIndividualRegionsDisplayed", width = 360, height = 600)
          ),
          
          
          DT::dataTableOutput("XY_Selected"))
  
}
