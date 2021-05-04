# Pipeline workhorse packages
library(targets)
library(tarchetypes)
library(future)
library(future.callr)

# See https://books.ropensci.org/targets/hpc.html#future for parallelizing
plan(callr)

# Set target-specific options such as packages.
tar_option_set(
  packages = c(
    "here",
    "data.table",
    "magrittr",
    "JM",
    "JMbayes2",
    "ggplot2",
    "mstate",
    "nlme",
    "kableExtra"
  )
)

# Source functions
source("data-raw/prepare-raw-data.R")
source("R/individual-cell-models.R")

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
  #tar_target(cell_lines, paste0(c("CD8", "CD4", "NK"), "_abs_log")),
  tar_target(
    CD8_model_noInter,
    fit_indiv_JM(
      long_outcome = "CD8_abs_log",
      datasets = prepare_JM_data(dat_merged, "CD8_abs_log", admin_cens_time = 7),
      fform = ~ trans2 + trans3 + trans4 - 1
    )
  ),
  tar_target(
    CD8_model_Inter,
    fit_indiv_JM(
      long_outcome = "CD8_abs_log",
      datasets = prepare_JM_data(dat_merged, "CD8_abs_log", admin_cens_time = 7),
      fform = ~ trans2 + ATG.2:trans2 + trans3 + ATG.3:trans3 + trans4 - 1
    )
  ),
  tar_target(
    CD4_model_noInter,
    fit_indiv_JM(
      long_outcome = "CD4_abs_log",
      prepare_JM_data(dat_merged, "CD4_abs_log", admin_cens_time = 7),
      fform = ~ trans2 + trans3 + trans4 - 1
    )
  ),
  tar_target(
    CD4_model_Inter,
    fit_indiv_JM(
      long_outcome = "CD4_abs_log",
      prepare_JM_data(dat_merged, "CD4_abs_log", admin_cens_time = 7),
      fform = ~ trans2 + ATG.2:trans2 + trans3 + ATG.3:trans3 + trans4 - 1
    )
  ),
  tar_target(
    NK_model_noInter,
    fit_indiv_JM(
      long_outcome = "NK_abs_log",
      prepare_JM_data(dat_merged, "NK_abs_log", admin_cens_time = 7),
      fform = ~ trans2 + trans3 + trans4 - 1
    )
  ),
  tar_target(
    NK_model_Inter,
    fit_indiv_JM(
      long_outcome = "NK_abs_log",
      prepare_JM_data(dat_merged, "NK_abs_log", admin_cens_time = 7),
      fform = ~ trans2 + ATG.2:trans2 + trans3 + ATG.3:trans3 + trans4 - 1
    )
  ),
  tarchetypes::tar_render(analysis_summary, path = "analysis/analysis-summary.Rmd")
)


