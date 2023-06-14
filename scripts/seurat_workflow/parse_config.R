library(rjson)
library(here)
library(purrr)

# config file

config <- fromJSON(file = here("data/config/config.json"))

# checks

if ("output_dir_name" %in% names(config)) {
  print(paste0("Output directory is: results/", config$output_dir_name))
} else {
  stop("config file needs var 'output_dir_name'")
}

if ("data_dir_name" %in% names(config)) {
  print(paste0("Data directory is: data/", config$data_dir_name))
} else {
  print("No data directory specified. Directory will be data/")
}

if ("sample_names" %in% names(config)) {
  print(paste0("Sample names are: ", reduce(config$sample_names, paste)))
} else {
  stop("config file needs var 'sample_names'")
}

if ("markers_filename" %in% names(config)) {
  print(paste0("Marker file name is: ", config$markers_filename))
} else {
  stop("config file needs var 'markers_filename'")
}

if ("normalization_method" %in% names(config)) {
  print(paste0("Normalization method is: ", config$normalization_method))
} else {
  stop("config file needs var 'normalization_method'")
}

if ("integration_method" %in% names(config)) {
  if (!(config$integration_method %in% c("integrate", "merge"))) {
    stop("integration method not valid; must be one of 'integrate' or 'merge'")
  } else {
    print(paste0("Integration method is: ", config$integration_method))
  }
} else {
  stop("config file needs var 'integration_method'")
}

if ("iteration" %in% names(config)) {
  print(paste0("Iteration number is: ", config$iteration))
} else {
  config$iteration <- 1
  print("Iteration number not specified; default to 1")
}

if ("louvain_res" %in% names(config)) {
  print(paste0("Chosen Louvain resolution is: ", config$louvain_res))
} else {
  config$iteration <- 0.5
  print("Louvain resolution not specified; default to 0.5")
}

# dependent parameters 

if (config$integration_method == "merge") {
  config$clustering_assay <- "SCT"
} else if(config$integration_method == "integrate") {
  config$clustering_assay <- "integrated"
} else {
  stop("Invalid integration method; cannot set clustering assay")
}
