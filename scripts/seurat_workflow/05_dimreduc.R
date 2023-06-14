# Dimensionality Reduction
# Amy Huang
# 11/24/2020

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
update_record(record_path, "\nRun script 05_dimreduc.R")

# read in data ----

seurat_obj <- readRDS(
  file.path(obj_dir, paste0("seurat_", config$integration_method, ".rds"))
)

# dimred workflow ----

if (config$clustering_assay == "SCT") {
  VariableFeatures(seurat_obj[["SCT"]]) <- rownames(seurat_obj[["SCT"]]@scale.data)
}

seurat_obj <- seurat_obj %>% 
  set_assay(config$clustering_assay) %>% 
  RunPCA(assay = config$clustering_assay, 
         reduction.name = paste0(config$clustering_assay, "_pca"), 
         reduction.key = paste0(config$clustering_assay, "PC_")) %>% 
  RunUMAP(reduction = paste0(config$clustering_assay, "_pca"), dims = 1:30) 

saveRDS(seurat_obj, file.path(obj_dir, "seurat_obj.rds"))
update_record(
  record_path, 
  paste0("objects/seurat_obj.rds \t",
         "Overwrite with ", config$integration_method, "d Seurat object, dimred complete")
)
