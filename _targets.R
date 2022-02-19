# Global pipeline set-up --------------------------------------------------

# See following for setting up final pipline:
# https://github.com/epiforecasts/evaluate-delta-for-forecasting/blob/4d8449fbedba71690ac6f1320a8438bdc10a4f44/_targets.R
# - Also use capsule for renv

# Workhorse packages
library("targets")
library("tarchetypes")

# All packages used by the projects - this is not good for renv
project_pkgs <- c(
  "data.table",
  "JM",
  "JMbayes2",
  "ggplot2",
  "mstate",
  "nlme",
  "kableExtra"
)

tar_option_set(packages = project_pkgs, error = "continue")
# Uncomment if running scripts interactively:
# sapply(project_pkgs, function(pkg) require(pkg, character.only = TRUE)); rm(project_pkgs)


# Miscellaneous objects ---------------------------------------------------


#...


# Analysis pipeline -------------------------------------------------------


# Source support functions
source("data-raw/prepare-raw-data.R")
source("R/modelling-helpers.R")
source("R/plotting-helpers.R")

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
  tar_target(
    reference_values,
    cbind.data.frame(
      "cell_type" = c("CD3", "CD4", "CD8", "NK", "CD19"),
      "lower_limit" = c(860, 560, 260, 40, 60),
      "upper_limit" = c(2490, 1490, 990, 390, 1000)
    )
  ),

  # Prepare datasets for analysis (pre-DLI ones)
  tar_target(
    NMA_preDLI_datasets,
    get_preDLI_datasets(
      dat_merged[TCD2 %in% c("NMA RD: ALT", "UD: ALT + ATG")], # NMA UD: ALT" are weird protocol
      admin_cens = 6
    )
  )

  # Post-DLI datasets here - or put in other file
  # tar_target(
  #   NMA_postDLI_datasets,
  #   get_postDLI_datasets(
  #     dat_merged = dat_merged[TCD2 %in% c("NMA RD: ALT", "UD: ALT + ATG")],
  #     admin_cens_dli = 12
  #   )
  # )

  #tarchetypes::tar_render(analysis_summary, path = "analysis/2020-09_analysis-summary.Rmd")
  #tarchetypes::tar_render([and rmd with raw data visualisations, also after data prep..
  #.. interactive with plotly?])
  # Or with shiny??
)

# Source cohort-specific targets
source("R/NMA-preDLI-models.R")
source("R/NMA-postDLI-models.R")

#...

targets_list <- c(targets_list, preDLI_targets, postDLI_targets)
targets_list
