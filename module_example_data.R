server_example_data = function(input, output, session, set_features_file, set_features_name, set_bigwig){
  #examples
  observeEvent(input$ExampleMCF7_bza, {
    showNotification(ui = "This data is for 4 histone marks from the MCF7 cell line treated with bezadoxifene.", duration = 10, id = "Note_ExMCF7_bza", type = "warning")
    #setting set_features_file, set_features_name, and set_bigwig is sufficient for valid setup
    set_features_file(paste0(bed_path, "/MCF7_bza_500ext.bed"))
    set_features_name("MCF7_bza_500ext")
    MCF7bza_bws = dir("/slipstream/galaxy/production/galaxy-dist/static/UCSCtracks/breast/MCF7_drug_treatments_pooled_inputs", pattern = "MCF7_bza_.+_FE.bw", full.names = T)
    names(MCF7bza_bws) = sub("/slipstream/galaxy/production/galaxy-dist/static/UCSCtracks/breast/MCF7_drug_treatments_pooled_inputs/MCF7_bza_", "", MCF7bza_bws) %>%
      sub(pattern = "_FE.bw", replacement = "")
    example_bw = data.frame(filename = names(MCF7bza_bws), filepath = MCF7bza_bws, stringsAsFactors = F)
    set_bigwig(example_bw)
  })
  observeEvent(input$ExampleKasumi, {
    showNotification(ui = "This data is Dan Trombly's ChIPseq in Kasumi cell lines. For Kaleem.", duration = 10, id = "Note_ExKasumi", type = "warning")
    #setting set_features_file, set_features_name, and set_bigwig is sufficient for valid setup
    set_features_file(paste0(bed_path, "/Kasumi_bivalency.bed"))
    set_features_name("Kasumi_AE_bivalency")
    ex_bws = dir("/slipstream/galaxy/production/galaxy-dist/static/UCSCtracks/aml/hg38/kasumi/DT_aml-eto_hg38/", pattern = "_FE.bw", full.names = T)
    names(ex_bws) = basename(ex_bws) %>% sub("Kasumi1_", "", .) %>%
      sub(pattern = "_pooled_FE.bw", replacement = "")
    example_bw = data.frame(filename = names(ex_bws), filepath = ex_bws, stringsAsFactors = F)
    set_bigwig(example_bw)
  })
  observeEvent(input$ExampleBivalency, {
    showNotification(ui = "Bivalency within 2kb of TSSes.", duration = 10, id = "Note_ExBiv", type = "warning")
    #setting set_features_file, set_features_name, and set_bigwig is sufficient for valid setup
    set_features_file(paste0(bed_path, "/TSS_serial_bivalency_2kb_ext.bed"))
    set_features_name("TSS_2kb_ext")
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
    example_bw = data.frame(filename = names(ex_bws), filepath = ex_bws, stringsAsFactors = F)
    set_bigwig(example_bw)
  })
}