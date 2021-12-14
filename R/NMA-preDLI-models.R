# Prepare Slurm (figure out how to get this in another file)
login <- future::tweak(future::remote, workers = "direct_clust")
slurm_run <- future::tweak(
  strategy = future.batchtools::batchtools_slurm,
  template = "~/future_slurm.tmpl",
  resources = list(
    walltime = 59, # in minutes
    ntasks = 1,
    ncpus = 1,
    memory = '4G',
    partition = 'short'
  )
)

slurm_three_cores <- future::tweak(
  strategy = future.batchtools::batchtools_slurm,
  template = "~/future_slurm.tmpl",
  resources = list(
    walltime = 180, # in minutes
    ntasks = 1,
    ncpus = 3,
    memory = '4G'#,
    #partition = 'short'
  )
)

# List of targets NMA
NMA_targets <- list(
  tar_target(
    NMA_preDLI_CD3_long,
    run_preDLI_longitudinal(
      cell_line = "CD3_abs_log",
      form_fixed = "ns(intSCT2_5, 3) * hirisk * ATG + CMV_PD", # add three-way interaction later
      form_random = list(IDAA = pdDiag(~ ns(intSCT2_5, 3))),
      dat = NMA_preDLI_datasets$long
    )
  ),
  tar_target(
    NMA_preDLI_cox,
    run_preDLI_cox(
      form = Surv(time, status) ~ hirisk.1 + ATG.2 + hirisk.2 + ATG.3 + strata(trans),
      dat_wide = NMA_preDLI_datasets$wide
    )
  ),
  tar_target(
    NMA_preDLI_CD3_jointModel_both,
    command = value(
      future({
        run_jointModel(
          long_obj = NMA_preDLI_CD3_long,
          surv_obj = NMA_preDLI_cox,
          timeVar = "intSCT2_5",
          parameterization = "both",
          iter.EM = 750,
          interFact = list(
            "value" = ~ strata(trans) - 1,
            "slope" = ~ strata(trans) - 1
          ),
          derivForm = list(
            fixed = ~ 0 + dns(intSCT2_5, 3),
            random = ~ 0 + dns(intSCT2_5, 3),
            indFixed = c(2:4),
            indRandom = c(2:4)
          )
        )
      })
    ),
    deployment = "worker",
    resources = tar_resources(future = tar_resources_future(plan = list(login, slurm_run)))
  ),
  # Try jmbayes2
  tar_target(
    NMA_preDLI_survdat,
    command = {
      dat_wide <- NMA_preDLI_datasets$wide
      tmat <- trans.comprisk(K = 3, names = c("gvhd", "relapse", "other_nrf"))
      covs <- c("CMV_PD", "hirisk", "ATG", "earlylow_DLI",
                "HLAmismatch_GvH", "relation", "SCT_May2010")

      msdat <- msprep(
        time = c(NA, rep("endpoint7", 3)),
        status = with(
          dat_wide, cbind(
            NA,
            1 * (endpoint7_s == "gvhd"),
            1 * (endpoint7_s == "relapse"),
            1 * (endpoint7_s == "other_nrf")
          )
        ),
        data = data.frame(dat_wide),
        trans = tmat,
        keep = covs,
        id = "IDAA"
      )

      expand.covs(msdat, covs, append = TRUE, longnames = FALSE)
    }
  ),
  tar_target(
    NMA_preDLI_CD3_jmbayes2_both,
    #value(future({
    #  c(future::availableCores(), parallel::detectCores()) # there does seem to be 3 cores available
    #})),
    command = value(
      future({
        jm(
          Surv_object = NMA_preDLI_cox,
          Mixed_objects = list(NMA_preDLI_CD3_long),
          time_var = "intSCT2_5",
          functional_forms = ~ value(CD3_abs_log) + slope(CD3_abs_log) +
            (value(CD3_abs_log) + slope(CD3_abs_log)):(strata(trans) - 1),
          data_Surv = NMA_preDLI_survdat, # Add later as intermediate data..
          control = list("n_burnin" = 5000, "n_iter" = 15000, "cores" = 3)
        )
      })
    ),
    deployment = "worker", # set deployment = main? parallel on main machine?
    resources = tar_resources(future = tar_resources_future(plan = list(
      login, slurm_three_cores)))#, future.callr::callr)))
  )
)

# See https://github.com/ropensci/targets/discussions/359

# Optional if all cell lines share same submodel
# tar_map(
#   list(cells = c("CD3_abs_log", "CD4_abs_log", "CD8_abs_log")),
#   tar_target(
#     NMA_preDLI_long,
#     run_preDLI_longitudinal(
#       cell_line = cells,
#       form_fixed = "ns(intSCT2_5, 3) * hirisk + ATG + CMV_PD",
#       form_random = list(IDAA = pdDiag(~ ns(intSCT2_5, 3))),
#       dat = NMA_preDLI_datasets$long
#     )
#   ),
#   unlist = FALSE
# ),
