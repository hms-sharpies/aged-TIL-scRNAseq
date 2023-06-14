# Integrate 
# Amy Huang
# 12/2/2020


# libraries ----

library(tidyverse)
library(here)
library(Seurat)
library(myfunctions)
source(here("scripts/seurat_workflow/parse_config.R"))
source(here("scripts/seurat_workflow/update_record.R"))


# directory paths and begin record ----

output_dir <- here("results", config$output_dir_name)
obj_dir <- file.path(output_dir, "objects")
record_path <- here("results", config$output_dir_name, "record.txt")
update_record(record_path, "\nRun script 04_integrate.R")


# read in data ----

seurat_singlets <- readRDS(file.path(obj_dir, "seurat_singlets.rds"))

# Normalization ----

if (config$normalization_method == "lognormalize") {
  
  print("not yet implemented, sorry :(")
  
} else if (config$normalization_method == "sctransform") {
  
  seurat_singlets <- map(seurat_singlets, SCTransform)
  
}


# combination method ----

if (config$integration_method == "merge") {
  
  seurat_obj <- merge(
    seurat_singlets[[1]], 
    y = seurat_singlets[2:4], 
    add.cell.ids = names(seurat_singlets), 
    merge.data = TRUE
  )
  saveRDS(
    seurat_obj, 
    file.path(obj_dir, paste0("seurat_", config$integration_method, ".rds"))
  )
  update_record(
    record_path, 
    "objects/seurat_merged.rds \tSeurat normalized (RNA), merged object"
  )
  
} else if (integration_method == "integrate") {
  
  seurat_singlets <- map(seurat_singlets, set_assay, "SCT")
  integrate_seurat_list <- function(seurat_list) {
    sc_features <- SelectIntegrationFeatures(object.list = seurat_list, 
                                             nfeatures = 3000)
    seurat_list <- PrepSCTIntegration(object.list = seurat_list, 
                                      anchor.features = sc_features)
    sc_anchors <- FindIntegrationAnchors(object.list = seurat_list, 
                                         normalization.method = "SCT", 
                                         anchor.features = sc_features)
    rm(seurat_list, sc_features)
    IntegrateData(anchorset = sc_anchors, 
                  normalization.method = "SCT")
  }
  seurat_obj <- integrate_seurat_list(seurat_singlets) 
  saveRDS(seurat_obj, 
          file.path(obj_dir, paste0("seurat_", integration_method, ".rds")))
  update_record(
    record_path, 
    "objects/seurat_integrated.rds \tSeurat normalized (RNA), integrated object"
  )
  
}
