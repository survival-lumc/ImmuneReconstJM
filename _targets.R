library(targets)
library(tarchetypes)

# Set target-specific options such as packages.
tar_option_set(imports = "ImmuneReconstJM")

# Source functions not in R/
source("data-raw/prepare-raw-data.R")

# Pipeline
list(
  tar_target(
    lymphocytes_raw,
    data.table::data.table(readRDS("data-raw/2021-04-09_v4/lymphocytes.rds"))
  ),
  tar_target(
    variables_raw,
    data.table::data.table(readRDS("data-raw/2021-04-09_v4/variables.rds"))
  ),
  tar_target(dat_merged, prepare_raw_data(lymphocytes_raw, variables_raw)),
  tarchetypes::tar_render("20210406_meeting", path = "analysis/2021-04-06_meeting.Rmd")
)


