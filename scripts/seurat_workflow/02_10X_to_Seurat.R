# 10X to Seurat
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
update_record(record_path, "\nRun script 02_10X_to_Seurat.R")

if ("data_dir" %in% names(config)) {
  data_dir <- here("data", data_dir)
} else {
  data_dir <- here("data")
}

# read in data ----

data_list <- map(
  config$sample_names, 
  ~Read10X(file.path(data_dir,..1, "outs/filtered_feature_bc_matrix"))
)
names(data_list) <- config$sample_names
saveRDS(data_list, file.path(obj_dir, "data_list.rds"))
update_record(
  record_path, 
  "objects/data_list.rds \tList of 10X data objects from filtered_feature_bc_matrix"
)

# create seurat objects ----

seurat_list <- imap(
  data_list, 
  ~CreateSeuratObject(count = ..1[["Gene Expression"]], 
                      project = ..2, 
                      min.cells = 3)
)
saveRDS(seurat_list, file.path(obj_dir, "seurat_list_raw.rds"))
update_record(
  record_path, 
  "objects/seurat_list_raw.rds \tList of raw Seurat objects from filtered_feature_bc_matrix"
)

# QC RNA feature counts ----

add_mito <- function(seurat) {
  seurat[["percent.mt"]] <- PercentageFeatureSet(seurat, pattern = "^mt-")
  return(seurat)
}
seurat_list <- map(seurat_list, add_mito)

imap(
  seurat_list, 
  ~..1 %>% 
    set_idents("orig.ident") %>% 
    VlnPlot(features = c("nCount_RNA", "nFeature_RNA", "percent.mt")) %>% 
    ggsave(file.path(output_dir, paste0("qc_RNAcount_vln_", ..2, ".png")), 
           plot = .)
)
update_record(
  record_path, 
  "qc_RNAcount_vln_{sample}.png \tQuality control plots showing nCount_RNA, nFeature_RNA, percent.mt"
)

seurat_list <- map(
  seurat_list, 
  ~..1 %>% 
    set_assay("RNA") %>% 
    subset(subset = nFeature_RNA > 200 & 
             nFeature_RNA < 6000 & 
             percent.mt < 10) %>% 
    NormalizeData() %>% 
    FindVariableFeatures() %>% 
    ScaleData() %>% 
    RunPCA(reduction.name = "rnaPCA", 
           reduction.key = "rnaPC_")
)

saveRDS(seurat_list, file.path(obj_dir, "seurat_list.rds"))
update_record(
  record_path, 
  "objects/seurat_list.rds\tList of Seurat objects, filtered, normalized, and scaled by RNA"
)
