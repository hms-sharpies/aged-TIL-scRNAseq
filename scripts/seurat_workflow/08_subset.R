# Subset out groups of interest
# Amy Huang
# 1/4/2021

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
update_record(record_path, "\nRun script 08_subset.R")

# read in data ----

seurat_obj <- readRDS(file.path(obj_dir, "seurat_obj.rds"))

# subset ----

DefaultAssay(seurat_obj) <- "RNA"
clustering_res <- 0.3
clustering_ident <- paste0(config$clustering_assay, "_snn_res.", clustering_res)
Idents(seurat_obj) <- clustering_ident

subset_clusters <- c("7", "8", "9")
seurat_subset <- subset(seurat_obj, idents = subset_clusters, invert = TRUE)

# clean up the subset

seurat_subset@assays[[config$clustering_assay]] <- NULL
clustering_sols <- sprintf("%s_snn_res.%s", config$clustering_assay, seq(0.1, 0.9, 0.1))
for (i in clustering_sols) {
  seurat_subset@meta.data[[i]] <- NULL
}
seurat_subset@reductions[[paste0(config$clustering_assay, "_pca")]] <- NULL
seurat_subset@reductions[["umap"]] <- NULL

# save subset

saveRDS(seurat_subset, file.path(obj_dir, "seurat_subset.rds"))
update_record(record_path, paste0(
  "objects/seurat_subset.rds \t", 
  "Seurat object, subset with clustering res ", clustering_res, 
  " and clusters ", reduce(subset_clusters, paste)
))

seurat_obj <- seurat_subset
