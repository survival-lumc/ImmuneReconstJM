# Global pipeline set-up --------------------------------------------------

# See following for setting up final pipline:
# https://github.com/epiforecasts/evaluate-delta-for-forecasting/blob/4d8449fbedba71690ac6f1320a8438bdc10a4f44/_targets.R
# - Also use capsule for renv

# Workhorse packages
library("targets")
library("tarchetypes")
library("future")
library("future.batchtools")
library("future.callr")

# All packages used by the projects - this is not good for renv
project_pkgs <- c(
  "data.table",
  "JM",
  "JMbayes2",
  "ggplot2",
  "mstate",
  "nlme",
  "kableExtra",
  "ggrepel",
  "future"
)

tar_option_set(packages = project_pkgs, error = "continue")
# Uncomment if running scripts interactively:
# sapply(project_pkgs, function(pkg) require(pkg, character.only = TRUE)); rm(project_pkgs)

# Everything except specific target is sequential
plan(sequential) # plan(callr) instead if local parallel


# Miscellaneous objects ---------------------------------------------------


# Reference values for the cells
cell_reference_values <- cbind.data.frame(
  "cell_type" = c("CD3", "CD4", "CD8", "NK", "CD19"),
  "lower_limit" = c(860, 560, 260, 40, 60),
  "upper_limit" = c(2490, 1490, 990, 390, 1000)
)


# Analysis pipeline -------------------------------------------------------


# Source support functions
source("data-raw/prepare-raw-data.R")
source("R/modelling-helpers.R")

# Pipeline (parts are in separate files):
targets_list <- list(

  # Part 1: Data preparation
  tar_target(
    lymphocytes_raw,
    data.table::data.table(readRDS("data-raw/2021-11-17_v8/lymphocytes.rds"))
  ),
  tar_target(
    variables_raw,
    data.table::data.table(readRDS("data-raw/2021-11-17_v8/variables.rds"))
  ),
  tar_target(dat_merged, prepare_raw_data(lymphocytes_raw, variables_raw)),
  tar_target(reference_values, data.table(cell_reference_values)),

  # Prepare datasets for analysis (pre-DLI ones)
  tar_target(
    NMA_preDLI_datasets,
    get_preDLI_datasets(
      # <CHECK THIS IS RIGHT?> 171?
      dat_merged[TCD2 %in% c(#"NMA UD: ALT", - these are weird protocol
                             "NMA RD: ALT", "UD: ALT + ATG")],
      admin_cens = 6
    )
  )
  #tarchetypes::tar_render(analysis_summary, path = "analysis/2020-09_analysis-summary.Rmd")
)

# Source cohort-specific targets
source("R/NMA-preDLI-models.R")
#...

targets_list <- c(targets_list, NMA_targets)
targets_list
