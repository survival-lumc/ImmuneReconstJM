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
#rm(project_pkgs)


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
    data.table::data.table(readRDS("data-raw/2021-09-22_v7/lymphocytes.rds"))
  ),
  tar_target(
    variables_raw,
    data.table::data.table(readRDS("data-raw/2021-09-22_v7/variables.rds"))
  ),
  tar_target(dat_merged, prepare_raw_data(lymphocytes_raw, variables_raw)),
  tar_target(reference_values, data.table(cell_reference_values)),

  # For now just non-myeloablative
  tar_target(
    datasets,
    get_datasets(
      dat_merged[TCD2 %in% c("NMA RD: ALT", "UD: ALT + ATG")], # only NMA cohort now
      admin_cens = 24
    )
  ),
  #tar_target(datasets, get_datasets(dat_merged, admin_cens = 24)),
  tar_target(dli_msdata, prepare_dli_msdata(datasets$wide)),

  # Part 2: Prepare submodels
  tar_target(
    long_submodels,
    run_longitudinal_submodels(
      datasets$long,
      which_cells = c("CD4_abs_log", "CD8_abs_log", "CD3_abs_log"),
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
        DLI.3 + ATG.3 + #NRF_other submodel
        strata(trans)
    )
  ),

  # Part 3a: Run univariate joint models with both packages
  tar_target(
    JM_CD4_allDLI,
    run_jointModel(
      long_obj = long_submodels$CD4_abs_log,
      surv_obj = cox_all_dli,
      fform = ~ trans1 + trans2 + trans3 - 1,
      control = list("iter.EM" = 500)
    )
  ),
  tar_target(
    JM_CD8_allDLI,
    run_jointModel(
      long_obj = long_submodels$CD8_abs_log,
      surv_obj = cox_all_dli,
      fform = ~ trans1 + trans2 + trans3 - 1,
      control = list("iter.EM" = 500)
    )
  ),
  tar_target(
    JM_CD3_allDLI,
    run_jointModel(
      long_obj = long_submodels$CD3_abs_log,
      surv_obj = cox_all_dli,
      fform = ~ trans1 + trans2 + trans3 - 1,
      control = list("iter.EM" = 500)
    )
  )#,
  # Do NK and CD19 too?

  # Try bayesian ones??
  # tar_target(
  #   JM_CD4_allDLI_bayes,
  #   jm(
  #     Surv_object = cox_all_dli,
  #     Mixed_objects = list(long_submodels$CD4_abs_log),
  #     time_var = "intSCT2_5",
  #     functional_forms = ~ value(CD4_abs_log):(trans1 + trans2 + trans3) - 1,
  #     data_Surv = dli_msdata,
  #     control = list("n_burnin" = 5000, "n_iter" = 15000)
  #   )
  # )

  #tarchetypes::tar_render(analysis_summary, path = "analysis/2020-09_analysis-summary.Rmd")
)


