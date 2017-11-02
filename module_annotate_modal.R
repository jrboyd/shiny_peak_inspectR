
annotateModal <- function(failed = FALSE) {
  modalDialog(
    span('(Please annotate the selected sample)'),
    DT::dataTableOutput("DTPeaksAnnotate", width = "auto"),
    if (failed)
      div(tags$b("Hey man, this didn't work! Not sure why.", style = "color: red;")),
    
    footer = tagList(fluidRow(
      column(width = 8,
             uiOutput("DTPeaksAnnotateElements")),
      column(width = 4,
             actionButton("BtnCancelAnnotate", "Cancel"),
             actionButton("BtnConfirmAnnotate", "Confirm")
      )
    )
    ),
    size = "l",
    title = "Annotateing"
  )
}

#assumed colnames for MACS2 peak files
peak_cn = c("seqnames", "start", "end", "id", "score", "strand", "FE", "p-value", "q-value", "summit_pos")
get_col_classes.df = function(df){
  sapply(1:ncol(df), function(i)class(df[[i]]))
}
get_col_classes = function(file, skipFirst = F){
  df = read.table(file, stringsAsFactors = F, header = F, nrows = 10, skip = ifelse(skipFirst, 1, 0))
  get_col_classes.df(df)
}
test_for_header = function(peak_file){
  withFirst = get_col_classes(peak_file, skipFirst = F)
  skipFirst = get_col_classes(peak_file, skipFirst = T)
  #definitely no header
  if(all(withFirst == skipFirst)){
    return(F)
  }
  if(all(withFirst == "character")){
    return(T)
  }
  warning(paste("can't determine if", peak_file, "has a header, guess not."))
  return(F)
}
load_peak_wValidation = function(peak_file, with_notes = F){
  has_header = test_for_header(peak_file)
  df = read.table(peak_file, stringsAsFactors = F, header = has_header)
  col_classes = get_col_classes.df(df)
  if(!has_header){
    if(ncol(df) == length(peak_cn)){
      if(with_notes){
        showNotification("assuming file is narrowPeak.", type = "warning")
      }else{
        print("assuming file is narrowPeak.")
      }
      colnames(df) = peak_cn  
    }else{
      if(with_notes){
        showNotification("file not narrowPeak, loading as minimal bed file.", type = "warning")
      }else{
        print("file not narrowPeak, loading as minimal bed file.")
      }
      bed_cn = peak_cn[1:4]
      nc = min(ncol(df), length(bed_cn))
      colnames(df)[1:nc] = peak_cn[1:nc]
    }
  }else{#try to make colnames GRanges compatible
    if(all(col_classes[1:5] == c("character", "integer", "integer", "integer", "character"))){
      if(with_notes){
        showNotification("file looks like saved GRanges.", type = "warning")
      }else{
        print("file looks like saved GRanges.")
      }
      colnames(df)[1:5] = c("seqnames", "start", "end", "width", "strand")
    }else if(all(col_classes[1:ncol(df)] == c("character", "integer", "integer", 
                                              "character", "integer", "character", 
                                              "numeric", "numeric", "numeric", "integer")[1:ncol(df)])){
      if(with_notes){
        showNotification("file looks like bed or encode peak", type = "warning")
      }else{
        print("file looks like bed or encode peak")
      }
    }else{
      colnames(df)[1:3] = c("seqnames", "start", "end")
      if(with_notes){
        showNotification("forced to assume first 3 columns are minimal bed, might break.", type = "warning")
      }else{
        print("forced to assume first 3 columns are minimal bed, might break.")
      }
    }
  }
  print(head(df))
  print(paste0(nrow(df), " total rows..."))
  invisible(df)
}

server_annotateModal = function(input, output, session, get_annotateing_DF, set_annotateing_DF, roots_reference, load_file = load_peak_wValidation){
  AnnotateingDF = reactiveVal(NULL)
  
  observeEvent(input$BtnAnnotateSet, {
    if(is.null(get_annotateing_DF())){
      showNotification("No data selected.", type = "error")
      return()
    }
    #dataTable instance used to annotate set data after adding
    AnnotateingDF(get_annotateing_DF())
    showModal(annotateModal())
  })
  observeEvent(input$BtnCancelAnnotate, {
    AnnotateingDF(NULL)
    removeModal()
  })
  output$DTPeaksAnnotate = DT::renderDataTable({
    df = AnnotateingDF()
    if(is.null(df)) return(NULL)
    DT::datatable(df, 
                  filter = list(position = "top", clear = TRUE, plain = F),
                  options = list(
                    scrollX = T,
                    pageLength = 10), rownames = F)
  })
  
  shinyFileChoose(input, 'FilesAnnotationReference', 
                  roots= roots_reference, 
                  filetypes=c("bed", "gff", "gtf"))
  
  observeEvent(AnnotateingDF(), {
    df = AnnotateingDF()
    # showNotification("check update DTPeaksAnnotateElements")
    if(is.null(df)) return(NULL)
    output$DTPeaksAnnotateElements = renderUI({
      tagList(
        shinyFilesButton("FilesAnnotationReference", label = "Find Reference", multiple = F),
        radioButtons("RadioAnnotateFeature", choices = "promoter", selected = "promoter"),
        selectizeInput(inputId = "SelectAnnotateReference"),
        
        numericInput("NumAnnotateMaxDistance", "Max Annotation Distance", value = nrow(df), min = 0, max = 10^9, step = 10000),
        actionButton("BtnAnnotateWithGtf", label = "Annotate")
      )
    })
  })
  
  observeEvent(input$BtnAnnotateWithGtf, {
    showNotification("BtnAnnotateWithGtf", type = "message")
    df = AnnotateingDF()
    df = df[sample(input$DTPeaksAnnotate_rows_all, size = input$NumAnnotateNumberOfRegions), ]
    AnnotateingDF(df)
  })
  observeEvent(input$BtnConfirmAnnotate, {
    set_annotateing_DF(AnnotateingDF()[input$DTPeaksAnnotate_rows_all,])
    AnnotateingDF(NULL)
    removeModal()
    # showNotification("Confirm annotate.", type = "message")
  })
}
