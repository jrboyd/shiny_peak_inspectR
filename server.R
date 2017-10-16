# server.R definition
# source('JG_runx_intersect.R')
source("functions_process_bw.R")
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


shinyFiles2load = function(shinyF, roots){
  root_path = roots[shinyF$root]
  rel_path = paste0(unlist(shinyF$files), collapse = "/")
  file_path = paste0(root_path, "/", rel_path)
  return(file_path)
}

shinyFiles2save = function(shinyF, roots){
  root_path = roots[shinyF$root]
  rel_path = paste0(unlist(shinyF$name), collapse = "/")
  file_path = paste0(root_path, "/", rel_path)
  return(file_path)
}

roots_load <<- c("../peak_annotatR/intersectR_beds/", 
                 dir("/slipstream/galaxy/uploads/working/qc_framework", pattern = "^output", full.names = T))
names(roots_load) <- basename(roots_load)

server <- function(input, output){
  
  shinyFileChoose(input, 'FilesLoadSet', 
                  roots= roots_load, 
                  filetypes=c("bed", "txt", "Peak"))
  
  output$xy_values = renderPlotly({
    prof_dt = profiles_dt()
    x = prof_dt[sample == bw_toplot[1]]
    y = prof_dt[sample == bw_toplot[2]]
    x_summit_val = x[, .(xval = FE[which(FE == max(FE, na.rm = T))[1]]), by = .(hit)][order(as.numeric(hit))]
    y_summit_val = y[, .(yval = FE[which(FE == max(FE, na.rm = T))[1]]), by = .(hit)][order(as.numeric(hit))]
    xy_val = merge(x_summit_val, y_summit_val)
    xy_val$hit = as.integer(xy_val$hit)
    feat_gr = isolate(features_gr())
    feat_dt = as.data.table(feat_gr)
    feat_dt$id = as.integer(feat_dt$id)
    setkey(feat_dt, id)
    xy_val = cbind(xy_val, feat_dt[.(xy_val$hit),])
    xy_val$group = "neither"
    bw_toplot = c("MCF7_bza_H3K4AC", "MCF7_bza_H3K4ME3")
    print(xy_val)
    xy_val[get(bw_toplot[1]) & get(bw_toplot[2]), group := "both"]
    xy_val[get(bw_toplot[1]) & !get(bw_toplot[2]), group := bw_toplot[1]]
    xy_val[!get(bw_toplot[1]) & get(bw_toplot[2]), group := bw_toplot[2]]
    set.seed(0)
    #STOPPED HERE
    # p = ggplot(xy_val[sample(1:nrow(xy_val))]) + geom_point(aes(x = xval, y = yval, col = group))
    p = ggplot(xy_val) + geom_point(aes(x = xval, y = yval, col = group))
    ggplotly(p) %>% 
      layout(title = "the dat",
             dragmode =  "select")
    
    #scoring strategy
  })
  
  bw_toplot = c("H3K4AC", "H3K4ME3")
  grp_tocolor = c("MCF7_bza_H3K4AC", "MCF7_bza_H3K4ME3")
  # profiles_dt = reactiveVal(NULL, "profiles_dt")
  
  profiles_dt = reactive({NULL})
  #should be an observe event to process
  #   if(is.null(features_gr())) return(NULL)
  #   bw_file = "/slipstream/galaxy/production/galaxy-dist/static/UCSCtracks/breast/10A_progression/AF-MCF10AT1_RUNX1_pooled_FE.bw"
  #   MCF7bza_bws = dir("/slipstream/galaxy/production/galaxy-dist/static/UCSCtracks/breast/MCF7_drug_treatments_pooled_inputs", pattern = "MCF7_bza_.+_FE.bw", full.names = T)
  #   names(MCF7bza_bws) = sub("/slipstream/galaxy/production/galaxy-dist/static/UCSCtracks/breast/MCF7_drug_treatments_pooled_inputs/MCF7_bza_", "", MCF7bza_bws) %>%
  #     sub(pattern = "_FE.bw", replacement = "")
  #   fgr = features_gr()
  #   
  #   feature_size = ceiling(quantile(width(fgr), .75) / 1000) * 1000 #TODO feature size - handling features several times larger than average
  #   mids = start(fgr) + floor(width(fgr) / 2)
  #   start(fgr) = mids - feature_size / 2
  #   end(fgr) = mids + feature_size / 2 - 1
  #   win_size = 50 #TODO win_size
  #   if(is.null(fgr$id)) fgr$id = 1:length(fgr)#TODO id
  #   if(exists("out_dt")) remove(out_dt, envir = globalenv())
  #   for(tp in bw_toplot){
  #     bw_file = MCF7bza_bws[tp]
  #     bw_gr = fetch_windowed_bw(bw_file = bw_file, win_size = win_size, qgr = fgr)
  #     bw_dt = bw_gr2dt(bw_gr, qgr = fgr, win_size = win_size)
  #     bw_dt$sample = tp
  #     
  #     if(exists("out_dt")){
  #       out_dt = rbind(out_dt, bw_dt)
  #     }else{
  #       out_dt = bw_dt
  #     }
  #   }    
  #   return(out_dt)
  # })
  
  output$PlotTmp = renderPlot({
    if(is.null(profiles_dt())) return(NULL)
    win_size = 50 #TODO win_size
    gg_bw_banded_quantiles(profiles_dt(), win_size = win_size)
  })
  
  # features_gr = reactiveVal(NULL, "features_gr")
  features_file = reactiveVal(NULL, "features_file")
  features_name = reactiveVal(NULL, "features_name")
  
  features_gr = reactive({
    if(is.null(features_file())) return(NULL)
    GRanges(fread(features_file()))
  })
  
  output$SetPreview = DT::renderDataTable(
    {
      
      fgr = features_gr()
      if(is.null(fgr)){
        return(data.frame())
      }else{
        return(DT::datatable(as.data.frame(fgr), 
                             options = list(
                               pageLength = 5)))
      }
      
    })
  
  #temporary to drive reactivity
  observeEvent(profiles_dt(), {
    print(head(features_gr()))
  })
  
  observeEvent(input$FilesLoadSet, {
    print("server find file")
    file_path = shinyFiles2load(input$FilesLoadSet, roots_load)
    features_file(file_path)
    features_name(basename(file_path))
  })
  observeEvent(input$UploadLoadSet, {
    print("local find file")
    features_file(input$UploadLoadSet$datapath)
    features_name(input$UploadLoadSet$name)
  })
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
    output$PlotScatter2d <- renderPlotly({
      # print("plot1")
      
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
    output$PlotAggregateProfile <- renderPlotly({
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
      
      ggplotly(p) %>% 
        layout(title = "the dat",
               dragmode =  "select")
    })
  })
  
  
  
  # Coupled event 1
  output$PlotSelectedProfiles <- renderPlotly({
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