# VISION
# Amy Huang
# 11/10/2020

# Load libraries ----

library(tidyverse)
library(myfunctions)
library(here)
library(Seurat)
library(VISION)

# Set directories ----

prev_dir <- here("results/2021-01-04_workflow")
output_dir <- here("results/2021-01-05_pathway_analysis")
mkdir_if(output_dir)

# Load data ----

seurat_obj <- readRDS(file.path(prev_dir, "3_objects/seurat_obj.rds"))

# VISION ----

seurat_obj <- seurat_obj %>% 
  set_assay("RNA")

sigs <- c(
  here("data/signatures/c2.all.v7.2.symbols.gmt"), 
  here("data/signatures/c7.all.v7.2.symbols.gmt"), 
  here("data/signatures/h.all.v7.2.symbols.gmt"), 
  here("data/signatures/brian_exhaustion_sigs.gmt"), 
  here("data/signatures/LCMVGeneSets.gmt")
)

vision_obj <- Vision(seurat_obj, signatures = sigs, dimRed = "SCT_pca")
vision_obj <- analyze(vision_obj)
# viewResults(vision_obj)

saveRDS(vision_obj, file.path(output_dir, "vision_obj.rds"))
vision_obj <- readRDS(file.path(output_dir, "vision_obj.rds"))