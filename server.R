source("server_setup.R")
source("module_filter_modal.R")
source("module_annotate_modal.R")
server <- function(input, output, session){
  ###initialize
  js$disableTab("tab2")
  js$disableTab("tab3")
  
  shinyFileChoose(input, 'FilesLoadSet', 
                  roots= roots_load_set, 
                  filetypes=c("bed", "txt", "Peak"))
  
  shinyFileChoose(input, 'FilesLoadBigwig', 
                  roots= roots_load_bw, 
                  filetypes=c("bigwig", "bw"))
  
  output$ProfilesLoaded = renderUI({
    sample_names = unique(profiles_dt()$sample)
    x = sample_names[min(1, length(sample_names))]
    y = sample_names[min(2, length(sample_names))]
    fixedRow(
      column(width = 6,
             selectInput(inputId = "x_variable", "X-sample", choices = sample_names, selected = x)
      ),
      column(width = 6,
             selectInput(inputId = "y_variable", "Y-sample", choices = sample_names, selected = y)
      )
    )
  })
  
  output$GroupingsAvailable = renderUI({
    fgr = features_gr_filtered()
    chrms = unique(seqnames(fgr))
    mdat = elementMetadata(fgr)
    if(!is.null(mdat$id)) mdat$id = NULL
    if(ncol(mdat) > 0){
      col_classes = sapply(1:ncol(mdat), function(i)class(mdat[,i]))
      #groups composed of just T and F can be compared to create new sets
      logical_groups = sort(colnames(mdat)[which(col_classes == "logical")])
      #groups of factors can only be analyzed individually
      factor_names = colnames(mdat)[which(col_classes != "logical")]
      names(factor_names) = factor_names
      factor_groups = list(chromosomes = chrms)
      factor_groups = append(factor_groups, lapply(factor_names, function(x){
        unique(mdat[[x]])
      }))
    }else{
      logical_groups = character()
      factor_groups = list(chromosomes = chrms)
    }
    factor_groups = lapply(factor_groups, as.character)
    if(length(logical_groups) > 0){
      #create a conditional selector
      #the factor grouping chromosome (GRanges seqnames) is always assumed present
      #no radio button will be shown if logicals aren't present
      tagList(
        radioButtons(inputId = "GroupingType", label = "Groupings Available", choices = c("logical derived", "predefined"), selected = "logical derived"),
        conditionalPanel(
          condition = "input.GroupingType == 'logical derived'",
          selectInput("SelectLogicalGrouping", label = "Select Logical Groups", choices = logical_groups, multiple = T, selectize = T, selected = logical_groups[1:max(1, min(2, length(logical_groups)))])
        ),
        conditionalPanel(
          condition = "input.GroupingType == 'predefined'",
          selectInput("SelectFactorGrouping", label = "Select Factor Group", choices = names(factor_groups), selected = names(factor_groups)[1], multiple = F, selectize = F)
        )
      )
    }else{
      tagList(
        (radioButtons(inputId = "GroupingType", label = "Groupings Available", choices = c("predefined"), selected = "predefined")),
        selectInput("SelectFactorGrouping", label = "Select Factor Group", 
                    choices = names(factor_groups), selected = names(factor_groups)[1], 
                    multiple = F, selectize = F)
      )
      
    }
    
  })
  
  #set xy data when appropriate
  observe({
    if(is.null(input$GroupingType)) return(NULL)
    prof_dt = profiles_dt()
    #TODO selector for this
    # samples_loaded = unique(prof_dt$sample)
    x_variable = input$x_variable
    y_variable = input$y_variable
    if(is.null(x_variable) || is.null(y_variable)) return(NULL)
    x = prof_dt[sample == x_variable]
    y = prof_dt[sample == y_variable]
    x_summit_val = x[, .(xval = FE[which(FE == max(FE, na.rm = T))[1]]), by = .(hit)][order(as.numeric(hit))]
    y_summit_val = y[, .(yval = FE[which(FE == max(FE, na.rm = T))[1]]), by = .(hit)][order(as.numeric(hit))]
    xy_val = merge(x_summit_val, y_summit_val)
    xy_val$hit = as.integer(xy_val$hit)
    feat_gr = isolate(features_gr_filtered())
    feat_dt = as.data.table(feat_gr)
    feat_dt$id = as.integer(feat_dt$id)
    setkey(feat_dt, id)
    xy_val = cbind(xy_val, feat_dt[.(xy_val$hit),])
    
    if(input$GroupingType == "logical derived"){
      logical_groups = input$SelectLogicalGrouping
      
      switch(paste0("l", length(logical_groups)),
             l0 = {
               xy_val$plotting_group = "no plotting_group set"
             },
             l1 = {
               print(1)
               xy_val$plotting_group = "negative"
               xy_val[get(logical_groups[1]), plotting_group := logical_groups[1]]
             },
             l2 = {
               print(2)
               xy_val$plotting_group = "neither"
               xy_val[get(logical_groups[1]) & get(logical_groups[2]), plotting_group := "both"]
               xy_val[get(logical_groups[1]) & !get(logical_groups[2]), plotting_group := logical_groups[1]]
               xy_val[!get(logical_groups[1]) & get(logical_groups[2]), plotting_group := logical_groups[2]]
             },
             l3 = {
               print(3)
               xy_val$plotting_group = "none"
               xy_val[get(logical_groups[1]) & get(logical_groups[2]) & get(logical_groups[3]), plotting_group := "all"]
               xy_val[!get(logical_groups[1]) & get(logical_groups[2]) & get(logical_groups[3]), plotting_group := paste(logical_groups[2:3], collapse = " & ")]
               xy_val[get(logical_groups[1]) & !get(logical_groups[2]) & get(logical_groups[3]), plotting_group := paste(logical_groups[c(1,3)], collapse = " & ")]
               xy_val[get(logical_groups[1]) & get(logical_groups[2]) & !get(logical_groups[3]), plotting_group := paste(logical_groups[1:2], collapse = " & ")]
               xy_val[!get(logical_groups[1]) & !get(logical_groups[2]) & get(logical_groups[3]), plotting_group := logical_groups[3]]
               xy_val[get(logical_groups[1]) & !get(logical_groups[2]) & !get(logical_groups[3]), plotting_group := logical_groups[1]]
               xy_val[!get(logical_groups[1]) & get(logical_groups[2]) & !get(logical_groups[3]), plotting_group := logical_groups[2]]
             },
             {
               print("too many groups")
               xy_val$plotting_group = "too many groups"
             }
      )
    }else if(input$GroupingType == "predefined"){
      fgrp_name = input$SelectFactorGrouping
      if(fgrp_name == "chromosomes"){
        fgrp = xy_val$seqnames
      }else{
        fgrp = xy_val[[fgrp_name]]
      }
      xy_val$plotting_group = fgrp
      
    }else{
      stop(paste("bad GroupingType", input$GroupingType))
    }
    xy_plot_dt(xy_val)
  })
  
  output$xy_values = renderPlot({
    plotted_dt = visible_dt()
    if(is.null(plotted_dt)) return(NULL)
    #isolate because these trigger updates of xy_plot_dt() which already triggers reactivity
    x_variable = isolate(input$x_variable)
    y_variable = isolate(input$y_variable)
    p = ggplot(plotted_dt) + 
      geom_point(aes(x = xval, y = yval, col = plotting_group)) +
      labs(x = x_variable, y = y_variable, title = "Max FE in regions")
    p
  })
  
  bw_toplot = c("H3K4AC", "H3K4ME3")
  grp_tocolor = c("MCF7_bza_H3K4AC", "MCF7_bza_H3K4ME3")
  # profiles_dt = reactiveVal(NULL, "profiles_dt")
  
  profiles_dt = reactiveVal({NULL})
  
  #the sampled and unsampled data being plotted in the scatterplot
  xy_plot_dt = reactiveVal({NULL})
  
  seed = reactiveVal(0)
  
  #the visible scatterplot data after sampling
  visible_dt = reactive({
    xy_val = xy_plot_dt()
    if(is.null(xy_val)) return(NULL)
    n_displayed = min(2000, nrow(xy_val))
    set.seed(seed())
    xy_val[sample(x = 1:nrow(xy_val), size = n_displayed)]
  })
  
  #the data selected by brushing
  selected_dt = reactive({
    xy_dt = xy_plot_dt()
    if(is.null(xy_dt)) return(NULL)
    if(is.null(input$xy_brush)){
      sel_dt = xy_dt
    }else{
      sel_dt = brushedPoints(xy_dt, input$xy_brush)
    }
    return(sel_dt)
  })
  
  selected_profiles = reactive({
    p_dt = profiles_dt()
    setkey(p_dt, hit)
    
    hit_tp = selected_dt()$hit
    class(hit_tp) = class(p_dt$hit)
    p_dt[.(hit_tp)]
  })
  
  output$AggregateDisplayed = renderUI({
    sample_names = unique(profiles_dt()$sample)
    selectInput(inputId = "SelectAggDisplayed", label = "Aggregates Displayed", choices = sample_names, selected = c(input$x_variable, input$y_variable), selectize = T, multiple = T)
  })
  
  output$PlotAggregated = renderPlot({
    if(is.null(profiles_dt())) return(NULL)
    if(is.null(input$SelectAggDisplayed)) return(NULL)
    to_disp = input$SelectAggDisplayed
    win_size = 50 #TODO win_size
    p_dt = selected_profiles()
    if(nrow(p_dt) == 0) return(NULL)
    showNotification("agg plot", type = "message", duration = 2)
    ggplot_list = lapply(to_disp, function(sample_grp){
      p = gg_bw_banded_quantiles(p_dt[sample == sample_grp], win_size = win_size)
      p = p + labs(title = sample_grp)
      if(sample_grp != to_disp[length(to_disp)]) p = p + guides(fill = "none")
      return(p)
    })
    grid.arrange(grobs = ggplot_list, nrow = 1)
    
  })
  
  features_file = reactiveVal(NULL, "features_file")
  features_name = reactiveVal(NULL, "features_name")
  
  features_gr = reactiveVal(NULL)
  observeEvent(features_file(), {
    if(is.null(features_file())) return(NULL)
    dt = fread(features_file())
    colnames(dt)[1:3] = c("seqnames", "start", "end")
    if(is.null(dt$id)) dt$id = 1:nrow(dt)
    features_gr(GRanges(dt))
  })
  
  get_filtering_DF = function(){    
    if(is.null(features_gr())) return(NULL)
    as.data.frame(features_gr())
  }
  set_filtering_DF = function(new_df){
    features_gr(GRanges(new_df))
  }
  get_file_path = function(){
    features_file()
  }
  server_filterModal(input, output, session, 
                     get_filtering_DF = get_filtering_DF, 
                     set_filtering_DF = set_filtering_DF,
                     get_file_path = get_file_path)
  
  
  get_annotateing_DF = function(){    
    if(is.null(features_gr())) return(NULL)
    as.data.frame(features_gr())
  }
  set_annotateing_DF = function(new_df){
    features_gr(GRanges(new_df))
  }
  server_annotateModal(input, output, session, 
                     get_annotateing_DF = get_annotateing_DF, 
                     set_annotateing_DF = set_annotateing_DF,
                     roots_reference = roots_load_set)
  
  
  features_gr_filtered = reactive({
    fgr = features_gr()[input$SetPreview_rows_all]
    if(length(fgr) == 0) fgr = features_gr()
    fgr
  })
  
  output$SetPreview = DT::renderDataTable(
    {
      fgr = features_gr()
      if(is.null(fgr)){
        return(data.frame())
      }else{
        return(DT::datatable(as.data.frame(fgr),
                             # filter = list(position = "top", clear = TRUE, plain = F),
                             options = list(
                               pageLength = 5), rownames = F))
      }
    })
  
  observeEvent(input$SetPreview_rows_all, {
    print(length(input$SetPreview_rows_all))
  })
  
  observeEvent(input$FilesLoadSet, {
    print("server find file")
    file_path = shinyFiles2load(input$FilesLoadSet, roots_load_set)
    features_file(file_path)
    features_name(basename(file_path))
  })
  observeEvent(input$UploadLoadSet, {
    print("local find file")
    features_file(input$UploadLoadSet$datapath)
    features_name(input$UploadLoadSet$name)
  })
  
  bigwigSelected = reactiveVal()
  bigwigAdded = reactiveVal(data.frame(filename = character(), filepath = character()))
  
  observeEvent(input$FilesLoadBigwig, {
    file_path = shinyFiles2load(input$FilesLoadBigwig, roots_load_bw)
    updateTextInput(session, inputId = "TxtAddBigWig", value = basename(file_path))
    bigwigSelected(file_path)
  })
  
  observeEvent(input$BtnAddBigiwg, {
    if(is.null(bigwigSelected())) return(NULL)
    tmp = bigwigAdded()
    tmp = rbind(tmp, data.frame(filename = input$TxtAddBigWig, filepath = bigwigSelected()))
    bigwigAdded(tmp)
  })
  
  observeEvent(input$BtnFinishSetup, {
    js$enableTab("tab2")
    updateTabsetPanel(session = session, inputId = "navbar", selected = 'tab2')
  })
  
  observeEvent(input$navbar, {
    if(input$navbar == "tab2"){
      if(is.null(features_gr()) && nrow(bigwigAdded()) == 0){
        showNotification(ui = "No setup detected. Loading some example data!", duration = 10, id = "Note_ExDataWarning", type = "warning")
        print("no setup detected, loading some example data!")
        #setting features_file, features_name, and bigwigAdded is sufficient for valid setup
        features_file(paste0(bed_path, "/MCF7_bza_500ext.bed"))
        features_name("MCF7_bza_500ext")
        MCF7bza_bws = dir("/slipstream/galaxy/production/galaxy-dist/static/UCSCtracks/breast/MCF7_drug_treatments_pooled_inputs", pattern = "MCF7_bza_.+_FE.bw", full.names = T)
        names(MCF7bza_bws) = sub("/slipstream/galaxy/production/galaxy-dist/static/UCSCtracks/breast/MCF7_drug_treatments_pooled_inputs/MCF7_bza_", "", MCF7bza_bws) %>%
          sub(pattern = "_FE.bw", replacement = "")
        example_bw = data.frame(filename = names(MCF7bza_bws), filepath = MCF7bza_bws)
        bigwigAdded(example_bw)
      }
    }
  })
  
  
  
  output$AddedBigWigs = DT::renderDataTable(
    DT::datatable(bigwigAdded(), rownames = F))
  
  output$BedLength = renderText({
    fgr = features_gr_filtered()
    return(length(fgr))
  })
  
  output$BedSummary = DT::renderDataTable({
    fgr = features_gr_filtered()
    nchr = length(unique(seqnames(fgr)))
    col_classes = sapply(1:ncol(elementMetadata(fgr)), function(i){
      class(elementMetadata(fgr)[[i]])
    })
    col_counts = sapply(1:ncol(elementMetadata(fgr)), function(i){
      length(unique(elementMetadata(fgr)[[i]]))
    })
    mat = cbind(c("chr", colnames(elementMetadata(fgr))),
                c(nchr, col_counts))
    colnames(mat) = c("factor", "n")
    DT::datatable(mat)
  })
  
  output$BigWigSummary = DT::renderDataTable({
    bw_info = bigwigAdded()
    bw_sizes = file.size(as.character(bw_info$filepath))
    bw_sizes = sapply(bw_sizes, function(x)utils:::format.object_size(x, "auto"))
    bw_info$filepath = NULL
    bw_info$size = bw_sizes
    DT::datatable(bw_info, rownames = F)
  })
  
  observeEvent(input$BtnFinishProcess, {
    showNotification("Processing has begun!", type = "message")
    fgr = features_gr_filtered()
    bw_toplot = bigwigAdded()
    feature_size = ceiling(quantile(width(fgr), .75) / 1000) * 1000 #TODO feature size - handling features several times larger than average
    mids = start(fgr) + floor(width(fgr) / 2)
    start(fgr) = mids - feature_size / 2
    end(fgr) = mids + feature_size / 2 - 1
    win_size = 50 #TODO win_size
    fgr = fgr[order(fgr$id)]
    bed_dir = digest(fgr)
    win_dir = win_size
    cache_path = paste(bw_cache_path, bed_dir, win_dir, sep = "/")
    dir.create(cache_path, showWarnings = F, recursive = T)
    if(is.null(fgr$id)) fgr$id = 1:length(fgr)#TODO id
    if(exists("out_dt")) remove(out_dt, envir = globalenv())
    for(i in 1:nrow(bw_toplot)){
      bw_file = as.character(bw_toplot$filepath[i])
      bw_name = as.character(bw_toplot$filename[i])
      cache_file = paste0(basename(bw_file), "_", file.size(bw_file), ".RData") %>% 
        paste(cache_path, ., sep = "/")
      if(file.exists(cache_file)){
        showNotification(paste("Loading cached regions for", bw_name), type = "message", duration = 10)
        load(cache_file)
      }else{
        showNotification(paste("Fetching regions for", bw_name), type = "message", duration = 10)
        bw_gr = fetch_windowed_bw(bw_file = bw_file, win_size = win_size, qgr = fgr)
        bw_dt = bw_gr2dt(bw_gr, qgr = fgr, win_size = win_size)
        bw_dt$sample = bw_name
        save(bw_dt, file = cache_file)
      }
      if(exists("out_dt")){
        out_dt = rbind(out_dt, bw_dt)
      }else{
        out_dt = bw_dt
      }
    }
    profiles_dt(out_dt)
    js$enableTab("tab3")
    updateTabsetPanel(session = session, inputId = "navbar", selected = 'tab3')
  })
  
  output$XY_Selected <- DT::renderDataTable({
    #Add reactivity for selection
    selected_dt()
  })
  
  output$SelectIndividualRegionsDisplayed <- renderUI({
    profs = selected_profiles()
    if(nrow(profs) == 0) return(NULL)
    showNotification("update select UI")
    gid = profs$hit
    if(length(gid) > 50){
      gid = sample(gid, 50)
    }
    selectInput("SelectIndividual", 
                label = "Genes to plot", 
                choices = sort(gid), 
                selected = gid[sample(length(gid), min(8, length(gid)))], 
                multiple = T)
  })
  
  output$PlotIndividualRegionsDisplayed = renderPlot({
    hits = input$SelectIndividual
    if(is.null(hits))return(NULL)
    if(length(hits) == 0)return(NULL)
    profs = selected_profiles()
    if(nrow(profs) == 0) return(NULL)
    sel_profs = profs[hit %in% hits]
    if(nrow(sel_profs) == 0) return(NULL)
    ggplot(sel_profs) + geom_line(aes(x = x-1000, y = FE, color = sample)) + facet_grid(hit ~ .) +
      labs(x = "bp from center")
  })
  
  observeEvent(input$stop, {
    print("debug stop")
  })
  
  #examples
  observeEvent(input$ExampleMCF7_bza, {
    showNotification(ui = "This data is for 4 histone marks from the MCF7 cell line treated with bezadoxifene.", duration = 10, id = "Note_ExMCF7_bza", type = "warning")
    #setting features_file, features_name, and bigwigAdded is sufficient for valid setup
    features_file(paste0(bed_path, "/MCF7_bza_500ext.bed"))
    features_name("MCF7_bza_500ext")
    MCF7bza_bws = dir("/slipstream/galaxy/production/galaxy-dist/static/UCSCtracks/breast/MCF7_drug_treatments_pooled_inputs", pattern = "MCF7_bza_.+_FE.bw", full.names = T)
    names(MCF7bza_bws) = sub("/slipstream/galaxy/production/galaxy-dist/static/UCSCtracks/breast/MCF7_drug_treatments_pooled_inputs/MCF7_bza_", "", MCF7bza_bws) %>%
      sub(pattern = "_FE.bw", replacement = "")
    example_bw = data.frame(filename = names(MCF7bza_bws), filepath = MCF7bza_bws)
    bigwigAdded(example_bw)
  })
  observeEvent(input$ExampleKasumi, {
    showNotification(ui = "This data is Dan Trombly's ChIPseq in Kasumi cell lines. For Kaleem.", duration = 10, id = "Note_ExKasumi", type = "warning")
    #setting features_file, features_name, and bigwigAdded is sufficient for valid setup
    features_file(paste0(bed_path, "/Kasumi_bivalency.bed"))
    features_name("Kasumi_AE_bivalency")
    ex_bws = dir("/slipstream/galaxy/production/galaxy-dist/static/UCSCtracks/aml/hg38/kasumi/DT_aml-eto_hg38/", pattern = "_FE.bw", full.names = T)
    names(ex_bws) = basename(ex_bws) %>% sub("Kasumi1_", "", .) %>%
      sub(pattern = "_pooled_FE.bw", replacement = "")
    example_bw = data.frame(filename = names(ex_bws), filepath = ex_bws)
    bigwigAdded(example_bw)
  })
  observeEvent(input$ExampleBivalency, {
    showNotification(ui = "Bivalency within 2kb of TSSes.", duration = 10, id = "Note_ExBiv", type = "warning")
    #setting features_file, features_name, and bigwigAdded is sufficient for valid setup
    features_file(paste0(bed_path, "/TSS_serial_bivalency_2kb_ext.bed"))
    features_name("TSS_2kb_ext")
    ex_bws = c(
      "/slipstream/galaxy/uploads/working/qc_framework/output_bivalency_redo_patients_H7/H7_H3K4ME3_pooled/H7_H3K4ME3_pooled_FE.bw",
      "/slipstream/galaxy/uploads/working/qc_framework/output_bivalency_redo_patients_H7/H7_H3K27ME3_pooled/H7_H3K27ME3_pooled_FE.bw",
      "/slipstream/galaxy/uploads/working/qc_framework/output_bivalency_redo_patients_H7/patients_GOOD1-H3K4ME3_pooled/patients_GOOD1-H3K4ME3_pooled_FE.bw",
      "/slipstream/galaxy/uploads/working/qc_framework/output_bivalency_redo_patients_H7/patients_GOOD1-H3K27ME3_pooled/patients_GOOD1-H3K27ME3_pooled_FE.bw",
      "/slipstream/galaxy/uploads/working/qc_framework/output_bivalency_redo_patients_H7/patients_POOR5-H3K4ME3_pooled/patients_POOR5-H3K4ME3_pooled_FE.bw",
      "/slipstream/galaxy/uploads/working/qc_framework/output_bivalency_redo_patients_H7/patients_POOR5-H3K27ME3_pooled/patients_POOR5-H3K27ME3_pooled_FE.bw",
      "/slipstream/galaxy/production/galaxy-dist/static/UCSCtracks/breast/MCF7_drug_treatments_pooled_inputs/MCF7_ctrl_H3K4ME3_FE.bw",
      "/slipstream/galaxy/production/galaxy-dist/static/UCSCtracks/breast/MCF7_drug_treatments_pooled_inputs/MCF7_ctrl_H3K27ME3_FE.bw",
      "/slipstream/galaxy/production/galaxy-dist/static/UCSCtracks/breast/MCF10A_drug_treatments_pooled_inputs/MCF10A_ctrl_H3K4ME3_FE.bw",
      "/slipstream/galaxy/production/galaxy-dist/static/UCSCtracks/breast/MCF10A_drug_treatments_pooled_inputs/MCF10A_ctrl_H3K27ME3_FE.bw"
    )
    names(ex_bws) = basename(ex_bws) %>% sub("Kasumi1_", "", .) %>%
      sub(pattern = "_pooled_FE.bw", replacement = "") %>%
      sub(pattern = "ctrl_", replacement = "") %>%
      sub(pattern = "OOR", replacement = "") %>%
      sub(pattern = "OOD", replacement = "") %>%
      sub(pattern = "_FE.bw", replacement = "") %>%
      sub(pattern = "patients_", replacement = "") %>%
      sub(pattern = "-", replacement = "_")
    example_bw = data.frame(filename = names(ex_bws), filepath = ex_bws)
    bigwigAdded(example_bw)
  })
  
}