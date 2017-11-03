source("server_setup.R")
source("module_filter_modal.R")
source("module_annotate_modal.R")
source("module_rename_modal.R")
source("module_scatter_plot.R")
source("module_example_data.R")
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
  
  profiles_dt = reactiveVal({NULL})
  selected_dt = reactiveVal(data.table())
  
  set_selected = function(new_selected_dt){
    selected_dt(new_selected_dt)
  }
  
  
  server_scatter_plot(input, output, session, 
                        get_features = features_gr_filtered, 
                        get_profiles = profiles_dt, 
                        set_selected = set_selected)
  
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
    if(is.null(dt$id)){
      showNotification("No id column detected.  Using row number.  Either create your own or add a gene based id with Annotate button below table.", type = "warning", duration = 10)
      dt$id = paste0("region_", 1:nrow(dt))
    }
    if(class(dt$id) != "character") dt$id = as.character(dt$id)
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
                               scrollX = T,
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
  bigwigAdded = reactiveVal(data.frame(filename = character(), filepath = character(), stringsAsFactors = F))
  
  observeEvent(input$FilesLoadBigwig, {
    file_path = shinyFiles2load(input$FilesLoadBigwig, roots_load_bw)
    updateTextInput(session, inputId = "TxtAddBigWig", value = basename(file_path))
    bigwigSelected(file_path)
  })
  
  observeEvent(input$BtnAddBigwig, {
    if(is.null(bigwigSelected())) return(NULL)
    tmp = bigwigAdded()
    tmp = rbind(tmp, data.frame(filename = input$TxtAddBigWig, filepath = bigwigSelected()))
    updateTextInput(session, "TxtAddBigWig", value = "")
    bigwigAdded(tmp)
  })
  
  observeEvent(input$BtnFinishSetup, {
    bws = bigwigAdded()
    features = features_gr_filtered()
    if(is.null(bws)){
      showNotification("bigwigs NULL, can't proceed", type = "error")
      return()
    }
    if(nrow(bws) == 0){
      showNotification("no bigwigs have been added, can't proceed", type = "error")
      return()
    }
    if(is.null(features)){
      showNotification("features have not been loaded, can't proceed", type = "error")
      return()
    }
    if(length(features) == 0){
      showNotification("no valid features, can't proceed", type = "error")
      return()
    }
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
        example_bw = data.frame(filename = names(MCF7bza_bws), filepath = MCF7bza_bws, stringsAsFactors = F)
        bigwigAdded(example_bw)
      }
    }
  })
  
  
  
  output$AddedBigWigs = DT::renderDataTable(
    DT::datatable(bigwigAdded(), 
                  rownames = F, 
                  selection = "single", 
                  options = list(
                    scrollX = T,
                    pageLength = 5)))
  
  observeEvent(input$AddedBigWigs_rows_selected, {
    showNotification(as.character(input$AddedBigWigs_rows_selected))
    print(input$AddedBigWigs_rows_selected)
  })
  observeEvent(input$BtnRemoveBigWig, {
    if(is.null(input$AddedBigWigs_rows_selected)) return()
    showNotification(as.character(input$AddedBigWigs_rows_selected))
    bigwigAdded(bigwigAdded()[-input$AddedBigWigs_rows_selected,])
    
  })
  observeEvent(input$BtnRenameBigWig, {
    if(is.null(input$AddedBigWigs_rows_selected)) return()
    showNotification(as.character(input$AddedBigWigs_rows_selected))
    showModal(renameModal(as.character(bigwigAdded()[input$AddedBigWigs_rows_selected, ]$filename)))
  })
  
  server_renameModal(input, output, session, 
                     set_name = function(new_name){
                       df = bigwigAdded()
                       df[input$AddedBigWigs_rows_selected, ]$filename = new_name
                       bigwigAdded(df)
                     })
  
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
    DT::datatable(mat, 
                  options = list(
                    scrollX = T,
                    pageLength = 5))
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
  
  server_example_data(input, output, session, 
                      set_features_file = features_file, 
                      set_features_name = features_name, 
                      set_bigwig = bigwigAdded)
  
}