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
rm(project_pkgs)


# Miscellaneous objects ---------------------------------------------------


# Reference values for the cells
cell_reference_values <- cbind.data.frame(
  "cell_type" = c("CD3", "CD4", "CD8", "NK", "CD19"),
  "lower_limit" = c(860, 560, 260, 40, 60),
  "upper_limit" = c(2490, 1490, 990, 390, 1000)
)


# Analysis pipeline -------------------------------------------------------


# Source functions
source("data-raw/prepare-raw-data.R")
source("R/submodel-wrappers.R")
source("R/joint-model-wrappers.R")


# Pipeline:
list(

  # Part 1: Data preparation
  tar_target(
    lymphocytes_raw,
    data.table::data.table(readRDS("data-raw/2021-05-31_v6/lymphocytes.rds"))
  ),
  tar_target(
    variables_raw,
    data.table::data.table(readRDS("data-raw/2021-05-31_v6/variables.rds"))
  ),
  tar_target(dat_merged, prepare_raw_data(lymphocytes_raw, variables_raw)),
  tar_target(reference_values, data.table(cell_reference_values)),
  tar_target(datasets, get_datasets(dat_merged)),
  tar_target(dli_msdata, prepare_dli_msdata(datasets$wide)),

  # Part 2: Prepare submodels
  tar_target(
    long_submodels,
    run_longitudinal_submodels(
      datasets$long,
      which_cells = c("CD4_abs_log", "CD8_abs_log"),
      df_splines = 4,
      ranef_structure = "diagonal"
    )
  ),
  tar_target(
    cox_all_dli,
    run_cox_submodel(
      dli_msdata = dli_msdata,
      form = Surv(Tstart, Tstop, status) ~
        hirisk.1 + DLI.1 + ATG.1 + # Relapse submodel
        DLI.2 + ATG.2 + # GVHD submodel
        DLI.3 + ATG.3 + # NRF_other submodel
        strata(trans)
    )
  ),

  # Part 3a: Run univariate joint models with both packages
  tar_target(
    JM_CD4_allDLI_nointer,
    run_jointModel(
      long_obj = long_submodels$CD4_abs_log,
      surv_obj = cox_all_dli,
      fform = ~ trans1 + trans2 + trans3 - 1
    )
  ),
  tar_target(
    JM_CD4_allDLI_inter, # this we keep for illustration
    run_jointModel(
      long_obj = long_submodels$CD4_abs_log,
      surv_obj = cox_all_dli,
      fform = ~ trans1 + trans1:DLI.1 + trans2 + trans2:DLI.2 + trans3 + trans3:DLI.3 - 1
    )
  ),
  tar_target(
    JMbayes2_CD4_allDLI_nointer, # name is non case sensitive!
    jm(
      Surv_object = cox_all_dli,
      Mixed_objects = long_submodels$CD4_abs_log,
      time_var = "intSCT2_5",
      functional_forms = list(
        "CD4_abs_log" = ~ value(CD4_abs_log):(trans1 + trans2 + trans3 - 1)
      ),
      data_Surv = dli_msdata,
      n_iter = 15000L,
      n_burnin = 5000L#,
      # priors = list(
      #   # Local ridge priors for each alpha ~ N(0, 1/2) - change 3
      #   Tau_alphas = lapply(seq_len(3), function(alpha) matrix(data = 1)) # change
      # )
    )
  ),

  # Part 3b: Run bivariate jm (now alphas of CD8 and CD4 are in same model)
  tar_target(
    multivar_allDLI_nointer, # name is non case sensitive!
    jm(
      Surv_object = cox_all_dli,
      Mixed_objects = long_submodels,
      time_var = "intSCT2_5",
      functional_forms = list(
        "CD4_abs_log" = ~ value(CD4_abs_log):(trans1 + trans2 + trans3 - 1),
        "CD8_abs_log" = ~ value(CD8_abs_log):(trans1 + trans2 + trans3 - 1)
      ),
      data_Surv = dli_msdata, # try with this
      n_iter = 15000L,
      n_burnin = 5000L
    )
  ),

  tar_target(
    multivar_allDLI_slopeGVHD, # name is non case sensitive!
    jm(
      Surv_object = cox_all_dli,
      Mixed_objects = long_submodels,
      time_var = "intSCT2_5",
      functional_forms = list(
        "CD4_abs_log" = ~ value(CD4_abs_log):(trans1 + trans2 + trans3 - 1) +
          slope(CD4_abs_log):trans2,
        "CD8_abs_log" = ~ value(CD8_abs_log):(trans1 + trans2 + trans3 - 1) +
          slope(CD8_abs_log):trans2
      ),
      data_Surv = dli_msdata, # try with this
      n_iter = 15000L,
      n_burnin = 5000L
    )
  ),

  # Rest of cell-line hereafter, also multivar JMbayes2
  tarchetypes::tar_render(analysis_summary, path = "analysis/2020-09_analysis-summary.Rmd")
)


