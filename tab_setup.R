ui_tab_setup = function(){
  tagList( 
          h3("Pick a bed file (ideally after annotating with peak_annotatR)"),
          shinyFilesButton(id = "FilesLoadSet", label = "Find Files on Server", title = "Find Peaks to Annotate", multiple = F),
          fileInput(inputId = "UploadLoadSet", label = "Browse Local Files"),
          DT::dataTableOutput("SetPreview"),
          h3("Pick bigwigs to visualize at bed intervals"),
          shinyFilesButton(id = "FilesLoadBigwig", label = "Find bigwig on Server", title = "Find Peaks to Annotate", multiple = F),
          textInput("TxtAddBigWig", label = "Bigwig name"),
          actionButton(inputId = "BtnAddBigiwg", label = "Add Bigwig"),
          DT::dataTableOutput("AddedBigWigs"),
          actionButton("BtnFinishSetup", label = "Finish Setup")
  )
}
