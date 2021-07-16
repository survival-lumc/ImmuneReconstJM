# Global pipeline set-up --------------------------------------------------


# Workhorse packages
library(targets)
library(tarchetypes)
library(future)
library(future.callr)

# See https://books.ropensci.org/targets/hpc.html#future for parallelizing
plan(callr)

# All packages used by the projects
project_pkgs <- c(
  "checkmate",
  "data.table",
  "JM",
  "JMbayes2",
  "ggplot2",
  "mstate",
  "nlme",
  "kableExtra",
  "ggrepel"
)

tar_option_set(packages = project_pkgs)
# Uncomment if running scripts interactively:
# sapply(project_pkgs, function(pkg) require(pkg, character.only = TRUE))


# Analysis pipeline -------------------------------------------------------


# Source functions
source("data-raw/prepare-raw-data.R")
source("R/submodel-wrappers.R")
source("R/joint-model-wrappers.R")



# Pipeline
list(
  tar_target(
    lymphocytes_raw,
    data.table::data.table(readRDS("data-raw/2021-05-31_v6/lymphocytes.rds"))
  ),
  tar_target(
    variables_raw,
    data.table::data.table(readRDS("data-raw/2021-05-31_v6/variables.rds"))
  ),
  tar_target(dat_merged, prepare_raw_data(lymphocytes_raw, variables_raw)),
  tar_target(
    reference_values,
    data.table::data.table(
      cbind.data.frame(
        "cell_type" = c("CD3", "CD4", "CD8", "NK", "CD19"),
        "lower_limit" = c(860, 560, 260, 40, 60),
        "upper_limit" = c(2490, 1490, 990, 390, 1000)
      )
    )
  ),
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

  # Start of sub-pipeline
  tar_target(
    datasets,
    get_datasets(dat_merged)
  ),
  tar_target(
    dli_msdata,
    prepare_dli_msdata(datasets$wide)
  ),
  tar_target(
    long_submodels,
    run_longitudinal_submodels(datasets$long)
  ),
  tar_target(
    cox_submodel_all_dli,
    run_cox_submodel(
      dli_msdata = dli_msdata,
      form = Surv(Tstart, Tstop, status) ~
        hirisk.1 + DLI.1 + ATG.1 + # Relapse submodel
        DLI.2 + ATG.2 + # GVHD submodel
        DLI.3 + ATG.3 + # NRF_other submodel
        strata(trans)
    )
  ),
  tar_target(
    cox_submodel_no_dli, # Maybe also one with just DLI on gvhd trans
    run_cox_submodel(
      dli_msdata = dli_msdata,
      form = Surv(Tstart, Tstop, status) ~
        hirisk.1 + ATG.1 + # Relapse submodel
        ATG.2 + # GVHD submodel
        ATG.3 + # NRF_other submodel
        strata(trans)
    )
  ),
  tar_target(
    CD4_all_dli,
    run_JM(
      long_obj = long_submodels$CD4_abs_log,
      surv_obj = cox_submodel_all_dli,
      fform = ~ trans1 + trans1:DLI.1 + trans2 + trans2:DLI.2 +
        trans3 + trans3:DLI.3 - 1
    )
  ),
  tar_target(
    CD4_all_dli_avfform,
    run_JM(
      long_obj = long_submodels$CD4_abs_log,
      surv_obj = cox_submodel_all_dli,
      fform = ~ trans1 + trans2 + trans3 - 1
    )
  ),
  tar_target(
    CD8_all_dli,
    run_JM(
      long_obj = long_submodels$CD8_abs_log,
      surv_obj = cox_submodel_all_dli,
      fform = ~ trans1 + trans1:DLI.1 + trans2 + trans2:DLI.2 +
        trans3 + trans3:DLI.3 - 1
    )
  ),
  tar_target(
    CD8_all_dli_avfform,
    run_JM(
      long_obj = long_submodels$CD8_abs_log,
      surv_obj = cox_submodel_all_dli,
      fform = ~ trans1 + trans2 + trans3 - 1
    )
  )#,
  # Rest of cell-line hereafter, also multivar JMbayes2
  #tarchetypes::tar_render(analysis_summary, path = "analysis/analysis-summary.Rmd")
)


