# Clustering
# Amy Huang
# 11/24/2020

# libraries ----
library(tidyverse)
library(here)
library(Seurat)
library(myfunctions)
library(clustree)
source(here("scripts/seurat_workflow/parse_config.R"))
source(here("scripts/seurat_workflow/update_record.R"))

# directory paths and begin record ----

output_dir <- here("results", config$output_dir_name)
obj_dir <- file.path(output_dir, "objects")
record_path <- here("results", config$output_dir_name, "record.txt")
update_record(record_path, "\nRun script 06_clustering.R")

cluster_dir <- file.path(output_dir, "clustering")
mkdir_if(cluster_dir)
update_record(
  record_path, 
  "clustering/ \tDirectory for clustering results"
)

# read in data ----

seurat_obj <- readRDS(file.path(obj_dir, "seurat_obj.rds"))

# clustering workflow ----

seurat_obj <- FindNeighbors(seurat_obj, 
                            reduction = paste0(config$clustering_assay, "_pca"), 
                            dims = 1:30)

resolutions <- seq(0.1, 0.9, 0.1)
for (r in resolutions) {
  seurat_obj <- FindClusters(seurat_obj, resolution = r)
}

# save object ----

saveRDS(seurat_obj, file.path(obj_dir, "seurat_obj.rds"))

update_record(
  record_path, 
  "objects/seurat_obj.rds \tOverwrite with clustered Seurat object"
)

# cluster viz ----

cluster_names <- str_subset(names(seurat_obj@meta.data), "res")

map(
  cluster_names, 
  ~(seurat_obj %>% 
    set_assay("RNA") %>% 
    DimPlot(group.by = ..1, reduction = "umap", label = TRUE) + 
    ggtitle(..1) + 
    scale_color_manual(values = alison_palette)) %>% 
    ggsave(file.path(cluster_dir, paste0("umap_", ..1, ".png")), 
           plot = ., 
           height = 7, width = 7)
)
update_record(
  record_path, 
  "clustering/umap_{sample}.png \tPlot UMAP at different clustering resolutions"
)

# split cluster viz ----

map(
  cluster_names, 
  ~(seurat_obj %>% 
      set_assay("RNA") %>% 
      DimPlot(group.by = ..1, split.by = "orig.ident", 
              reduction = "umap", label = TRUE, ncol = 2) + 
      ggtitle(..1) + 
      scale_color_manual(values = alison_palette)) %>% 
    ggsave(file.path(cluster_dir, paste0("umap_split_orig.ident_", ..1, ".png")), 
           plot = ., 
           height = 7, width = 7)
)
update_record(
  record_path, 
  "clustering/umap_split_orig.ident_{sample}.png \tPlot UMAP at different clustering resolutions, split by sample identity"
)

# clustree ----

clustree(seurat_obj, prefix = paste0(config$clustering_assay, "_snn_res."))
ggsave(file.path(cluster_dir, "clustree.png"), height = 7, width = 7)
update_record(
  record_path, 
  "clustering/clustree.png \tPlot clustree"
)
