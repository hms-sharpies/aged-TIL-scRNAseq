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
obj_dir <- file.path(output_dir, "3_objects")
record_path <- here("results", config$output_dir_name, "record.txt")
update_record(record_path, "\nRun script 10_add_age_metadata.R")

# read in data ----

seurat_obj <- readRDS(file.path(obj_dir, "seurat_obj.rds"))

# add age metadata ----

metadata <- as_tibble(seurat_obj@meta.data, rownames = "cell")
age_metadata <- metadata %>% 
  select(cell, orig.ident) %>% 
  mutate(age = str_sub(orig.ident, end = 1)) %>% 
  select(cell, age) %>% 
  column_to_rownames("cell")

seurat_obj <- AddMetaData(seurat_obj, metadata = age_metadata)
saveRDS(seurat_obj, file.path(obj_dir, "seurat_obj.rds"))
update_record(
  record_path, 
  "3_objects/seurat_obj.rds \tOverwriting Seurat object with new object with age metadata"
)
