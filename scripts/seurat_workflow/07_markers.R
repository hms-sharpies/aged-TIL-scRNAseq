# Markers
# Amy Huang
# 11/24/2020

# libraries ----
library(tidyverse)
library(here)
library(Seurat)
library(myfunctions)
library(cowplot)
library(viridis)
source(here("scripts/seurat_workflow/parse_config.R"))
source(here("scripts/seurat_workflow/update_record.R"))

# directory paths and begin record ----

output_dir <- here("results", config$output_dir_name)
obj_dir <- file.path(output_dir, "3_objects")
record_path <- here("results", config$output_dir_name, "record.txt")
update_record(record_path, "\nRun script 07_markers.R")

marker_dir <- file.path(output_dir, "markers")
mkdir_if(marker_dir)
update_record(
  record_path, 
  "markers/ \tDirectory for markers from first iteration"
)

# read in data ----

seurat_obj <- readRDS(file.path(obj_dir, "seurat_obj.rds"))

# rna feature projection ----

markers <- read_tsv(here("data/markers", config$markers_filename))
seurat_obj <- set_assay(seurat_obj, "RNA")

rna_feature_proj_dir <- file.path(marker_dir, "rna_feature_proj")
mkdir_if(rna_feature_proj_dir)
update_record(record_path, "markers/rna_feature_proj/ \tDirectory for UMAP projections of marker genes")

map2(
  markers[[2]], 
  markers[[1]], 
  ~(FeaturePlot(seurat_obj, feature = ..1, reduction = "umap") + 
      ggtitle(sprintf("%s - %s", ..2, ..1)) + 
      scale_color_viridis()) %>% 
    ggsave(file.path(rna_feature_proj_dir, paste0(..2, "-", ..1, ".png")), 
           plot = ., 
           width = 7, height = 7)
)
update_record(
  record_path, 
  "markers/rna_feature_proj/*.png \tUMAP projections of marker genes"
)

rna_feature_proj_split_dir <- file.path(marker_dir, "rna_lymphocyte_feature_proj_split_age")
mkdir_if(rna_feature_proj_split_dir)
update_record(record_path, "markers/rna_lymphocyte_feature_proj_split_age/ \tDirectory for UMAP projections of marker genes, split by age")

add_title_and_scale <- function(ggprot, title1, title2) {
  ggprot + 
    ggtitle(sprintf("%s - %s", title2, title1)) + 
    scale_color_viridis()
}
map2(
  markers[[2]], 
  markers[[1]], 
  ~map(FeaturePlot(seurat_obj, feature = ..1, reduction = "umap", 
                split.by = "age", combine = FALSE),  
       add_title_and_scale, ..1, ..2) %>% 
    plot_grid(plotlist = .) %>% 
    ggsave(
      file.path(rna_feature_proj_split_dir, paste0(..2, "-", ..1, ".png")), 
      plot = ., 
      width = 10, height = 5)
)
update_record(
  record_path, 
  "markers/rna_lymphocyte_feature_proj_split_age/*.png \tUMAP projections of marker genes, split by sample"
)

# marker tables ----

# conserved markers----

DefaultAssay(seurat_obj) <- "RNA"
clustering_res <- c(0.3, 0.5, 0.8)

for (res in clustering_res) {
  clustering_ident <- paste0(config$clustering_assay, "_snn_res.", res)
  Idents(seurat_obj) <- clustering_ident
  groups <- seurat_obj@meta.data[[clustering_ident]] %>% unique() %>% 
    set_names()
  
  conserved_markers <- map_dfr(
    groups, 
    ~possibly(FindConservedMarkers, otherwise = tibble())(
      seurat_obj, 
      ..1, 
      grouping.var = "age"
    ) %>% as_tibble(rownames = "gene"), 
    .id = "cluster"
  )
  write_csv(conserved_markers, 
            file.path(marker_dir, paste0("res.", res, 
                                         "_conserved_markers.csv")))
  update_record(
    record_path, 
    paste0(
      "markers/res.", res, "_conserved_markers.csv \t", 
      "table of conserved genes for ", res, " clustering resolution"
    )
  )
}

# markers for each cluster ----

for (res in clustering_res) {
  cluster_markers <- seurat_obj %>% 
    set_idents(clustering_ident) %>% 
    possibly(FindAllMarkers, otherwise = tibble())(
      assay = "RNA", 
      logfc.threshold = 0, 
      only.pos = FALSE
    )
  write_csv(cluster_markers, file.path(marker_dir, paste0("res.", res, 
                                                          "_cluster_markers.csv")))
  update_record(
    record_path, 
    paste0(
      "markers/res.", res, "_conserved_markers.csv \t", 
      "table of marker genes for ", res, " clustering resolution"
    )
  )
}

# differential age markers, all clusters ----

age_markers <- seurat_obj %>% 
  set_idents("age") %>% 
  FindAllMarkers(
    assay = "RNA", 
    logfc.threshold = 0
  )
write_csv(age_markers, file.path(marker_dir, paste0("res.", res, 
                                                    "_age_markers.csv")))
update_record(
  record_path, 
  paste0(
    "markers/res.", res, "_age_markers.csv \t", 
    "table of marker genes for each age"
  )
)

# differential age markers, by cluster ----

for (res in clustering_res) {
  clustering_ident <- paste0(config$clustering_assay, "_snn_res.", res)
  seurat_obj <- set_idents(seurat_obj, clustering_ident)
  groups <- seurat_obj@meta.data[[clustering_ident]] %>% unique() %>% 
    set_names()

  diff_age_markers_by_cluster <- map_dfr(
    groups, 
    ~FindMarkers(
      seurat_obj, 
      assay = "RNA", 
      ident.1 = "O", 
      group.by = "age", 
      subset.ident = ..1, 
      logfc.threshold = 0
      ) %>% 
      rownames_to_column(var = "gene"), 
    .id = "cluster"
  )
  
  write_csv(diff_age_markers_by_cluster, 
            file.path(marker_dir, 
                      paste0("res.", res, "_diff_age_by_cluster.csv")))
  update_record(
    record_path, 
    paste0(
      "markers/res.", res, "_conserved_markers.csv \t", 
      "table of marker genes for ", res, 
      " clustering resolution, differential age by cluster"
    )
  )
}
