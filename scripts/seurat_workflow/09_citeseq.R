# citeseq analysis
# Amy Huang
# 11/24/2020

# libraries ----
library(tidyverse)
library(here)
library(Seurat)
library(myfunctions)
source(here("data/config/config.R"))

# directory paths and begin record ----

output_dir <- here("results", output_dir_name)
if (!dir.exists(output_dir)) stop("output_dir does not exist; stop. ")
obj_dir <- file.path(output_dir, "objects")
if (!dir.exists(obj_dir)) stop("obj_dir does not exist; stop. ")

record_file <- file.path(output_dir, "record.txt")
if (!file.exists(record_file)) stop("record.txt file does not exist; stop. ")
record_text <- c("\nRun script 08_citeseq.R")

citeseq_dir <- file.path(output_dir, "citeseq")
mkdir_if(citeseq_dir)
record_text <- c(
  record_text, 
  "Make citeseq_dir for citeseq analysis output"
)

# scale and normalize antibodies ----

seurat_obj <- seurat_obj %>% 
  set_assay("Ab") %>% 
  NormalizeData(assay = "Ab", normalization.method = "CLR") %>% 
  ScaleData(assay = "Ab")

saveRDS(seurat_obj, file.path(obj_dir, "seurat_obj.rds"))
record_text <- c(
  record_text, 
  "Overwriting seurat_obj.rds with citeseq normalization"
)

# Plot all antibodies ----

ab_names <- seurat_obj@assays[["Ab"]]@counts@Dimnames[[1]]

map(
  ab_names, 
  ~FeaturePlot(seurat_obj, features = ..1, 
               min.cutoff = "q10", max.cutoff = "q90") %>% 
    ggsave(file.path(citeseq_dir, paste0("umap_", ..1, ".png")), 
           plot = .)
)

record_text <- c(
  record_text, 
  "Plot citeseq antibodies projected on umap in citeseq/"
)

# antibody heatmap ---- 

Idents(seurat_obj) <- "SCT_snn_res.0.7"
DoHeatmap(seurat_obj, features = ab_names, assay = "Ab", angle = 90) + 
  NoLegend() + 
  ggtitle("SCT_snn_res.0.7")
ggsave(file.path(citeseq_dir, "heatmap.png"))

# add to records ----

write(record_text, 
      file = record_file, 
      append = TRUE)
