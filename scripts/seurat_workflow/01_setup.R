# set up results directory
# Amy Huang
# 11/23/2020

# libraries ----

library(tidyverse)
library(here)
library(myfunctions)
source(here("scripts/seurat_workflow/parse_config.R"))

# make directories ----

mkdir_if(here("results/"))

output_dir <- here("results", config$output_dir_name)
mkdir_if(output_dir)

obj_dir <- file.path(output_dir, "objects")
mkdir_if(obj_dir)

dir_create_text <- paste0(
  "Creating directories: \n",
  "results/", config$output_dir_name, "\tResults directory\n", 
  "results/", config$output_dir_name, "/objects/ \tObjects directory"
)

# create output records document ---- 

date_text <- paste0("Created on: ", Sys.Date(), "\n")
config_text <- map2_chr(names(config), 
                        as.character(config), ~paste(..1, "\t", ..2))
write(c(config$output_dir_name, date_text, 
        "Config settings:", config_text, "\n", 
        dir_create_text), 
      file = file.path(output_dir, "record.txt"))
