#!/usr/bin/env bash

R CMD BATCH 01_setup.R 
R CMD BATCH 02_10X_to_Seurat.R
R CMD BATCH 03_demux.R
R CMD BATCH 04_integrate.R
R CMD BATCH 05_dimreduc.R
R CMD BATCH 06_clustering.R
R CMD BATCH 07_markers.R

