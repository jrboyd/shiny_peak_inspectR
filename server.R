# server.R definition
# source('JG_runx_intersect.R')
server <- function(input, output){
  get_selected_plot_df = function(){
    event.data <- event_data("plotly_selected", source = "scatter")
    event_data("plotly_selected", source = "scatter")
    p_df = plotly_data(main_plot)
    #gene_type hardcoded as level, ugh
    grp_counts = table(factor(p_df$gene_type))
    addons = c(0, cumsum(grp_counts)[-length(grp_counts)]) + 1
    # If NULL dont do anything
    if(is.null(event.data) == T) return(NULL)
    idx = addons[event.data$curveNumber + 1] + event.data$pointNumber
    plot_df[idx,]
  }
  # plot_df = reactive({
  #   if(input$featureInput2 == "all"){
  #     plot_df = scat_prof  
  #   }else{
  #     plot_df = scat_prof[gene_type == input$featureInput2]
  #   }
  #   plot_df
  # })
  
  # Observes the second feature input for a change
  observeEvent(input$featureInput2,{
    
    if(input$featureInput2 == "all"){
      inner_df = scat_prof
    }else{
      inner_df = scat_prof[gene_type == input$featureInput2]
    }
    
    # Add column names
    # colnames(plot_df) <- c("x", "y", "Class")
    
    # Do a plotly contour plot to visualize the two featres with
    # the number of malignant cases as size
    # Note the use of 'source' argument
    output$Plot1 <- renderPlotly({
      print("plot1")
      
      # print(inner_df)
      bg = sub("#", "", rgb(t(col2rgb("snow2")/255)))
      set.seed(0)
      plot_df <<- inner_df#[sample(nrow(inner_df))]
      p = ggplot(plot_df, 
                 aes(x = RUNX1, y = RUNX2, color = gene_type)) + 
        geom_point(aes(label1 = gene_name, label2 = seqnames)) + 
        scale_color_manual(values = c("protein_coding" = "#e41a1c", "non_coding" = "#377eb8")) +
        coord_cartesian(xlim = c(0, max(scat_prof$RUNX1)), ylim = c(0, max(scat_prof$RUNX2))) +
        labs(title = "FE at strongest peak within 2kb of tss")
      # print(class(p))
      
      ply = ggplotly(p, tooltip = c("label1", "label2"), source = "scatter") %>%
        layout(title = "the dat",
               dragmode =  "select",
               plot_bgcolor = bg)
      main_plot <<- ply
      ply
      
    })
    
    # Create a contour plot of the number of malignant cases
    
    
    # Assign to parent environment
    # plot_df <<- plot_df
  })
  observeEvent(input$featureInput3,{
    output$Plot2 <- renderPlotly({
      # Get subset based on selection
      event.data <- event_data("plotly_selected", source = "scatter")
      # If NULL dont do anything
      if(is.null(event.data) == T) return(NULL)
      if(length(event.data) == 0) return(NULL)
      gid = get_selected_plot_df()$gene_id
      
      if(input$featureInput3 == "tss"){
        agg_prof = prof[gene_id %in% gid, .(y = mean(y)),by = .(x_tss_dist, group, index_se, gene_type)]
        p = ggplot(agg_prof) + geom_line(aes(x = x_tss_dist, y = y, color = group))
      }else{
        agg_prof = prof[gene_id %in% gid, .(y = mean(y)),by = .(x_summit_dist, group, index_se, gene_type)]
        p = ggplot(agg_prof) + geom_line(aes(x = x_summit_dist, y = y, color = group))
      }
      p = p + facet_grid(. ~ gene_type) + labs(x = paste("distance to", input$featureInput3), y = "FE")
      
      # p = ggplot(prof[gene_id %in%gid])
      # if(input$featureInput3 == "tss"){
      #   p = p + geom_line(aes(x = x_tss_dist, y = y, color = group))
      # }else{
      #   p = p + geom_line(aes(x = x_summit_dist, y = y, color = group))
      # }
      # p = p + facet_grid(gene_name ~ group) + annotate("line", x = c(-2000,2000), y = rep(ym, 2)) + 
      #   labs(x = "relative position to peak summit", y = "FE") + 
      #   theme(strip.text.y = element_text(size = 8, colour = "black", angle = 0))
      ggplotly(p)
    })
  })
  

  
  # Coupled event 1
  output$Plot3 <- renderPlotly({
    # Get subset based on selection
    event.data <- event_data("plotly_selected", source = "scatter")
    # If NULL dont do anything
    if(is.null(event.data) == T) return(NULL)
    if(is.null(input$List1) == T) return(NULL)
    gid = input$List1
    # gid = plot_df[event.data$pointNumber]$gene_id
    p = ggplot(prof[gene_name %in% gid])#[1:min(8, length(gid))]])
    if(input$featureInput3 == "tss"){
      p = p + geom_line(aes(x = x_tss_dist, y = y, color = group))
    }else{
      p = p + geom_line(aes(x = x_summit_dist, y = y, color = group))
    }
    p = p + facet_grid(gene_name ~ group) + annotate("line", x = c(-2000,2000), y = rep(ym, 2)) + 
      labs(x = paste("distance to", input$featureInput3), y = "FE") + 
      theme(strip.text.y = element_text(size = 8, colour = "black", angle = 0))
    ggplotly(p)
  })
  
  output$List1 <- renderUI({
    event.data <- event_data("plotly_selected", source = "scatter")
    # If NULL dont do anything
    if(is.null(event.data) == T) return(NULL)
    
    gid = get_selected_plot_df()$gene_name
    
    selectInput("List1", 
                label = "Genes to plot", 
                choices = sort(gid), 
                selected = gid[sample(length(gid), min(8, length(gid)))], 
                multiple = T)
  })
  
  
}