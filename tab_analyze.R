ui_tab_analyze = function(){
  tagList(h3("scatterplot of bigwig values"), 
          plotlyOutput("xy_values"),
          h3("aggregated summary of selected values"),
          plotOutput("PlotTmp"),
          h3("heatmap of all profiles"),
          h3("randomly selected individual profiles"))
}
