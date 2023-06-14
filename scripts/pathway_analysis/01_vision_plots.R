# VISION
# Amy Huang
# 11/10/2020

# Load libraries ----

library(tidyverse)
library(myfunctions)
library(here)
library(Seurat)
library(VISION)
library(viridis)
library(cowplot)

# Set directories ----

prev_dir <- here("results/2020-11-10_rerun_workflow/")
output_dir <- here("results/2020-11-11_pathway_analysis")
if (!dir.exists(output_dir)) dir.create(output_dir)

# plot vision scores with seurat obj ----

exh_scores <- as_tibble(vision_exh@SigScores, rownames = "cell") %>% 
  column_to_rownames("cell")
seurat_obj_scored <- AddMetaData(seurat_obj, exh_scores)

# Viz signatures ----

FeaturePlot(seurat_obj_scored, feature = "Proliferating_sig_from_LCMV", 
            min.cutoff = "q05", max.cutoff = "q95") + 
  scale_color_viridis()
FeaturePlot(seurat_obj_scored, feature = "Terminally_Exh_sig_from_LCMV", 
            min.cutoff = "q05", max.cutoff = "q95") + 
  scale_color_viridis()
FeaturePlot(seurat_obj_scored, feature = "Progenitor_Exh_sig_from_LCMV", 
            min.cutoff = "q05", max.cutoff = "q95") + 
  scale_color_viridis()

exh_sigs <- seurat_obj_scored@meta.data %>% names() %>% 
  str_subset("_sig_")

exh_sig_plots <- map(
  exh_sigs, 
  ~(FeaturePlot(seurat_obj_scored, 
                feature = ..1, 
                min.cutoff = "q05", 
                max.cutoff = "q95") + 
      scale_color_viridis()) %>% 
    ggsave(file.path(output_dir, paste0(..1, ".png")), 
           plot = ., 
           height = 7, width = 8)
)