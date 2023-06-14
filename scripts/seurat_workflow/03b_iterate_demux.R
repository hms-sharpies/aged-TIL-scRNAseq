# Demux
# Amy Huang
# 11/24/2020

# libraries ----
library(tidyverse)
library(here)
library(Seurat)
library(myfunctions)
library(ggrepel)
library(cowplot)
source(here("data/config/config.R"))
source(here("scripts/03_scRNAseq_workflow_reorganize_output/update_record.R"))

# directory paths and begin record ----

output_dir <- here("results", output_dir_name)
if (!dir.exists(output_dir)) stop("output_dir does not exist; stop. ")
obj_dir <- file.path(output_dir, "objects")
if (!dir.exists(obj_dir)) stop("obj_dir does not exist; stop. ")

record_path <- file.path(output_dir, "record.txt")
update_record(record_path, "\nRun script 03b_iterate_demux.R")

# read in data ----

seurat_list <- readRDS(file.path(obj_dir, "seurat_list.rds"))

# subset negatives ----

neg_dir <- file.path(output_dir, "hto_neg_dir")
mkdir_if(neg_dir)

seurat_hto_negs <- seurat_list %>% 
  map(
    ~..1 %>% 
      set_assay("HTO") %>% 
      set_idents("hash.ID") %>% 
      subset(idents = "Negative")
  )

# demux negatives ----

seurat_hto_negs <- map(seurat_hto_negs, HTODemux)

# neg demux viz ----

hto_idents <- c("orig.ident", "hash.ID", 
                "HTO_maxID", "HTO_classification.global")

plot_count_idents <- function(seurat, output_dir, idents) {
  seurat_metadata <- as_tibble(seurat@meta.data, rownames = "cell")
  for (ident in idents) {
    count_data <- seurat_metadata %>% 
      group_by(.data[[ident]]) %>% 
      summarise(count = n()) %>% 
      ggplot() + 
      aes(.data[[ident]], count, label = count) + 
      geom_col() + 
      geom_label() + 
      ggtitle(ident)
    ggsave(file.path(output_dir, 
                     sprintf("%s_counts_%s.png", 
                             seurat@project.name, 
                             ident)), 
           height = 7, width= 7)
  }
}
map(
  seurat_hto_negs, 
  ~plot_count_idents(..1, neg_dir, hto_idents)
)
update_record(record_path, "Plot HTO negative count idents")

plot_ridge_idents <- function(seurat, output_dir, idents) {
  for (ident in idents) {
    seurat %>% 
      set_assay("HTO") %>% 
      set_idents(ident) %>% 
      RidgePlot(features = rownames(seurat[["HTO"]]),
                ncol = 2) %>%
      ggsave(file.path(
        neg_dir, 
        sprintf("%s_ridge_%s.png", seurat@project.name, ident)
      ),
      plot = .)
  }
}
map(
  seurat_hto_negs, 
  ~plot_ridge_idents(..1, neg_dir, idents)
)
update_record(record_path, "Plot HTO negative ridge idents")

# merge iterative demux into original ----

seurat_orig_metadata <- map(seurat_list, 
                            ~as_tibble(..1@meta.data, rownames = "cell"))
seurat_neg_metadata <- map(seurat_hto_negs, 
                           ~as_tibble(..1@meta.data, rownames = "cell"))

hash.ID.iter_metadata <- map2(
  seurat_orig_metadata, 
  seurat_neg_metadata, 
  ~full_join(..1, ..2, by = "cell") %>% 
    mutate(hash.ID.iter = case_when(
      hash.ID.x == "Negative" ~ hash.ID.y, 
      TRUE ~ hash.ID.x
    )) %>% 
    select("cell", hash.ID.iter) %>% 
    column_to_rownames(var = "cell")
)

seurat_list <- map2(
  seurat_list, 
  hash.ID.iter_metadata, 
  ~AddMetaData(..1, ..2)
)

saveRDS(seurat_list, file.path(obj_dir, "seurat_list.rds"))
update_record(record_path, "Overwrite seurat_list.rds with neg iter")

# iterative demux viz ----

hto_iter_dir <- file.path(output_dir, "hto_iter")
mkdir_if(hto_iter_dir)

map(
  seurat_list, 
  ~plot_count_idents(..1, hto_iter_dir, hto_idents)
)
update_record(record_path, "Plot HTO iterative count idents")

map(
  seurat_list, 
  ~plot_ridge_idents(..1, hto_iter_dir, hto_idents)
)
update_record(record_path, "Plot HTO iterative ridge idents")

# extract singlets ----

seurat_singlets <- map(
  seurat_list, 
  ~..1 %>% 
    set_assay("HTO") %>% 
    set_idents("hash.ID.iter") %>% 
    subset(idents = unique(..1@meta.data$HTO_maxID)) %>% 
    set_assay("RNA")
)

saveRDS(seurat_singlets, file.path(
  obj_dir, 
  "seurat_singlets.rds"
))
update_record(record_path, "Created seurat_singlets.rds, after SCTransform")

