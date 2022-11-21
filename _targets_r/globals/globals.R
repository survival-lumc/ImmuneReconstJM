# Use capsule/renv eventually docker (!!)
options(tidyverse.quiet = TRUE)

# Global objects
cell_lines <- list("cell_line" = paste0("CD", c(3, 4, 8)))

# Load packages and functions necessary for pipeline
source(here::here("packages.R"))
source(here::here("data-raw/prepare-raw-data.R"))
functions <- list.files(here("R"), full.names = TRUE)
invisible(lapply(functions, source))
rm("functions")

# Run whole pipeline in parallel
plan(callr)

# Other pipeline options
tar_option_set(error = "continue")
