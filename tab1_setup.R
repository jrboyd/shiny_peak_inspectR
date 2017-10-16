ui_tab1_setup = function(){
  tagList(h3("Pick a bed file (ideally after annotating with peak_annotatR)"), 
          shinyFilesButton(id = "FilesLoadSet", label = "Find Files on Server", title = "Find Peaks to Annotate", multiple = F),
          fileInput(inputId = "UploadLoadSet", label = "Browse Local Files"),
          DT::dataTableOutput("SetPreview"),
          h3("Pick bigwigs to visualize at bed intervals"))
}
