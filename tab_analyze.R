ui_tab_analyze = function(){
  tagList(h3("scatterplot of bigwig values"), 
          uiOutput("ProfilesLoaded"),
          uiOutput("GroupingsAvailable"),
          plotlyOutput("xy_values"),
          h3("aggregated summary of selected values"),
          uiOutput("AggregateDisplayed"),
          plotOutput("PlotAggregated"),
          h3("heatmap of all profiles"),
          h3("randomly selected individual profiles"),
          DT::dataTableOutput("XY_Selected"))
          
}
