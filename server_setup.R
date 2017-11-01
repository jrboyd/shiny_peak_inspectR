bw_cache_path = "~/ShinyApps/shiny_peak_data/cached_profiles"
bed_path = "~/ShinyApps/shiny_peak_data/beds/"
names(bed_path) = "intersectR"
source("functions_process_bw.R")

shinyFiles2load = function(shinyF, roots){
  root_path = roots[shinyF$root]
  rel_path = paste0(unlist(shinyF$files), collapse = "/")
  file_path = paste0(root_path, "/", rel_path)
  return(file_path)
}

require(magrittr)
user_roots = dir("/slipstream/home/", full.names = T) %>% paste0(. , "/ShinyData")
user_roots = subset(user_roots, dir.exists(user_roots))
names(user_roots) = dirname(user_roots) %>% basename()
qcframework_load <<- dir("/slipstream/galaxy/uploads/working/qc_framework", pattern = "^output", full.names = T)
names(qcframework_load) <- basename(qcframework_load)
roots_load_set = c(bed_path, user_roots, qcframework_load)

# names(roots_load_set) <- basename(roots_load_set)
roots_load_bw <<- c("/slipstream/galaxy/production/galaxy-dist/static/UCSCtracks/", 
                    dir("/slipstream/galaxy/uploads/working/qc_framework", pattern = "^output", full.names = T))
names(roots_load_bw) <- basename(roots_load_bw)
roots_load_bw = c(user_roots, roots_load_bw)




