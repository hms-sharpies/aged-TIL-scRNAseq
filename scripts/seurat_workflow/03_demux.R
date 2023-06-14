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
source(here("scripts/seurat_workflow/parse_config.R"))
source(here("scripts/seurat_workflow/update_record.R"))

# directory paths and begin record ----

output_dir <- here("results", config$output_dir_name)
obj_dir <- file.path(output_dir, "objects")
record_path <- here("results", config$output_dir_name, "record.txt")
update_record(record_path, "\nRun script 03_demux.R")

output_dir <- here("results", config$output_dir_name)
demux_dir_name <- "0_demux"
demux_dir <- file.path(output_dir, demux_dir_name)
mkdir_if(demux_dir)
update_record(
  record_path, 
  paste0(
    demux_dir_name, 
    "\tDirectory created for demultiplexing results"
  )
)

# read in data ----

data_list <- readRDS(file.path(obj_dir, "data_list.rds"))
seurat_list <- readRDS(file.path(obj_dir, "seurat_list.rds"))

# Split HTOs and antibodies ----

split_htos_from_antibodies <- function(data, nhto, nabs) {
  oligo_names <- data[["Antibody Capture"]]@Dimnames[[1]]
  hto_names <- str_subset(oligo_names, "hashtag.")
  ab_names <- str_subset(oligo_names, "hashtag.", negate = TRUE)
  data[["HTO"]] <- data[["Antibody Capture"]][hto_names,]
  if (length(ab_names) == 0) {
    data[["Antibody Capture"]] <- NULL
  } else {
    data[["Antibody Capture"]] <- data[["Antibody Capture"]][ab_names,]
  }
  return(data)
}
data_list <- map(data_list, split_htos_from_antibodies)

add_antibody_assays <- function(seurat, data) {
  assay_names <- names(data)
  assay_names <- assay_names[-str_which(assay_names, "Gene Expression")]
  for (name in assay_names) {
    assay_data <- data[[name]][, colnames(x = seurat)]
    seurat[[name]] <- CreateAssayObject(assay_data)
    seurat <- NormalizeData(
      seurat, 
      assay = name, 
      normalization.method = "CLR"
    )
  }
  return(seurat)
}
seurat_list <- map2(seurat_list, data_list, add_antibody_assays)

# Selecting quantile for demultiplexing ----

get_hto_classification_for_quantile <- function(seurat, q) {
  seurat_demux <- HTODemux(seurat, assay = "HTO", positive.quantile = q)
  table(seurat_demux$HTO_classification.global) %>% 
    as.matrix() %>% t() %>% as_tibble()
}

get_quantiles_df <- function(seurat, quantiles) {
  names(quantiles) <- quantiles
  map_dfr(quantiles, 
          ~get_hto_classification_for_quantile(seurat, .), 
          .id = "quantile") %>% 
    pivot_longer(cols = c("Doublet", "Negative", "Singlet"), 
                 names_to = "HTO_classification.global", 
                 values_to = "count")
}   

plot_quantiles_df <- function(quantiles_df) {
  max_singlet <- quantiles_df %>% 
    filter(`HTO_classification.global` == "Singlet") %>% 
    arrange(desc(`count`)) %>% 
    head(1)
  ggplot(quantiles_df) + 
    aes(quantile, count, color = `HTO_classification.global`) + 
    geom_point() + 
    facet_wrap(~`HTO_classification.global`, scale="free_y") + 
    geom_label_repel(data = max_singlet, aes(label = quantile))
}

find_max_singlet_quantile <- function(
  seurat,         
  quantile_range = c(0.9, 0.99), 
  save_plot = FALSE, 
  plot_title = "Testing quantiles", 
  filepath = here()
) {
  quantiles <- seq(quantile_range[1], quantile_range[2], 0.01)
  quantiles_df <- get_quantiles_df(seurat, quantiles)
  if (save_plot == TRUE) {
    quantiles_plot <- plot_quantiles_df(quantiles_df) + 
      ggtitle(plot_title)
    ggsave(
      filepath, 
      plot = quantiles_plot, 
      width = 7, 
      height = 7
    )
  }
  quantiles_df %>% 
    filter(`HTO_classification.global` == "Singlet") %>% 
    arrange(desc(`count`)) %>% 
    head(1) 
}

max_singlet_quantile <- map2_dfr(
  seurat_list, 
  names(seurat_list), 
  ~find_max_singlet_quantile(
    ..1, 
    save_plot = TRUE, 
    plot_title = sprintf("Testing quantiles, sample %s", ..2), 
    filepath = file.path(
      demux_dir, 
      sprintf("plot_quantiles_%s.png", ..2)
    )
  )
) %>% 
  mutate(quantile = as.numeric(quantile))

update_record(
  record_path, 
  paste0(
    demux_dir_name, "/plot_quantiles_{sample}.png \t", 
    "Plot showing singlet maximizing quantile threshold"
  )
)

write_csv(max_singlet_quantile, 
          file.path(demux_dir, "max_singlet_quantile.csv"))
update_record(
  record_path, 
  paste0(
    demux_dir_name, "/max_singlet_quantile.csv \t", 
    "Table of singlet maximizing quantile threshold"
  )
)

# Demultiplexing ----

seurat_list <- map2(seurat_list, 
                    max_singlet_quantile$quantile, 
                    ~HTODemux(..1, assay = "HTO", positive.quantile = ..2))

# save seurat_list ----

saveRDS(seurat_list, file.path(obj_dir, "seurat_list.rds"))
# seurat_list <- readRDS(file.path(obj_dir, "seurat_list.rds"))
update_record(
  record_path, 
  "objects/seurat_list.rds \tWrite over Seurat object; quantile selected, demux completed"
)

# extract singlets ----

extract_singlets <- function(seurat) {
  DefaultAssay(seurat) <- "HTO"
  Idents(seurat) <- "HTO_classification.global"
  seurat_singlets <- subset(seurat, idents = "Singlet")
  DefaultAssay(seurat_singlets) <- "RNA"
  seurat_singlets 
}

seurat_singlets <- map(seurat_list, extract_singlets)

saveRDS(seurat_singlets, file.path(obj_dir, "seurat_singlets.rds"))
update_record(
  record_path, 
  "objects/seurat_singlets.rds \tSeurat object; subset singlets after HTODemux"
)

# HTO ID viz ----

# ridge plots

viz_idents <- c("orig.ident", "hash.ID", "HTO_maxID", "HTO_classification.global")
for (ident in viz_idents) {
  imap(
    seurat_list, 
    ~..1 %>% 
      set_assay("HTO") %>% 
      set_idents(ident) %>% 
      RidgePlot(features = rownames(..1[["HTO"]]),
                ncol = 2) %>%
      ggsave(file.path(demux_dir, sprintf("ridge_hto_%s_%s.png", ..2, ident)),
             plot = .)
  )
}
update_record(
  record_path,
  paste0(
    demux_dir_name, "/ridge_hto_{sample}_{ident}.png \tPlots, HTO idents, ridge"
  )
)

# heatmap

imap(
  seurat_list, 
  ~(HTOHeatmap(..1) + ggtitle(paste0("Sample ", ..2))) %>% 
      ggsave(file.path(demux_dir, paste0("heatmap_", ..2, ".png")), 
             plot = .)
)
# note that HTOHeatmap downsamples
update_record(
  record_path,
  paste0(
    demux_dir_name, "/heatmap_{sample}.png \tPlot heatmap, note downsampling"
  )
)

# feature scatter

feature_scatter_dir <- file.path(demux_dir, "feature_scatter")
mkdir_if(feature_scatter_dir)

update_record(
  record_path,
  paste0(
    demux_dir_name, "/feature_scatter \tDirectory for feature scatter of HTO pairs"
  )
)

all_hto_feature_scatter <- function(seurat, name, ident) {
  seurat <- seurat %>% 
    set_assay("HTO") %>% 
    set_idents(ident) 
  seurat_htos <- rownames(seurat@assays$HTO)
  cross(list(
      hto1 = seurat_htos, 
      hto2 = seurat_htos
    )) %>% map(
      ~FeatureScatter(seurat, 
                      feature1 = ..1[[1]], 
                      feature2 = ..1[[2]], 
                      pt.size = 0.5) %>% 
        ggsave(file.path(feature_scatter_dir, 
                         sprintf("%s_%s_%s_%s.png", name, ident, ..1[[1]], ..1[[2]])), 
               plot = .)
    )
}
imap(seurat_list, ~all_hto_feature_scatter(..1, ..2, "HTO_maxID"))
# imap(seurat_list, ~all_hto_feature_scatter(..1, ..2, "HTO_classification.global"))
      
update_record(
  record_path,
  paste0(
    demux_dir_name, "/feature_scatter/* \tPlot, feature scatter of HTO pairs"
  )
)

# save tables ----

map_dfr(seurat_list , ~..1@meta.data$HTO_maxID %>% table(), .id = "sample") %>% 
  write_csv(file.path(demux_dir, "hto_id.csv"))

map_dfr(seurat_list , ~..1@meta.data$HTO_classification.global %>% table(), 
        .id = "sample") %>% 
  write_csv(file.path(demux_dir, "hto_classification.csv"))

update_record(
  record_path,
  paste0(
    demux_dir_name, "/hto_id.csv \tTable with HTO ID cell counts"
  )
)
update_record(
  record_path,
  paste0(
    demux_dir_name, "/hto_classification.csv \tTable with HTO classification cell counts"
  )
)