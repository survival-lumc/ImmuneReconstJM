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

# List of targets NMA
preDLI_targets <- list(
  tar_target(
    preDLI_CD3_long,
    run_preDLI_longitudinal(
      cell_line = "CD3_abs_log",
      form_fixed = "ns(intSCT2_5, 3) * hirisk * ATG + CMV_PD", # add three-way interaction later
      form_random = list(IDAA = pdDiag(~ ns(intSCT2_5, 3))),
      dat = NMA_preDLI_datasets$long
    )
  ),
  tar_target(
    preDLI_CD4_long,
    run_preDLI_longitudinal(
      cell_line = "CD4_abs_log",
      form_fixed = "ns(intSCT2_5, 3) * hirisk * ATG + CMV_PD", # add three-way interaction later
      form_random = list(IDAA = pdDiag(~ ns(intSCT2_5, 3))),
      dat = NMA_preDLI_datasets$long
    )
  ),
  tar_target(
    preDLI_CD8_long,
    run_preDLI_longitudinal(
      cell_line = "CD8_abs_log",
      form_fixed = "ns(intSCT2_5, 3) * hirisk * ATG + CMV_PD", # add three-way interaction later
      form_random = list(IDAA = pdDiag(~ ns(intSCT2_5, 3))),
      dat = NMA_preDLI_datasets$long
    )
  ),
  tar_target(
    preDLI_cox,
    run_preDLI_cox(
      form = Surv(time, status) ~ ATG.1 + hirisk.1 + ATG.2 + hirisk.2 + ATG.3 + strata(trans),
      dat_wide = NMA_preDLI_datasets$wide
    )
  ),
  tar_target(
    preDLI_CD3_jointModel_both,
    command = value(
      future({
        run_jointModel(
          long_obj = preDLI_CD3_long,
          surv_obj = preDLI_cox,
          timeVar = "intSCT2_5",
          parameterization = "both",
          iter.EM = 1000,
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
  tar_target(
    preDLI_CD4_jointModel_both,
    command = value(
      future({
        run_jointModel(
          long_obj = preDLI_CD4_long,
          surv_obj = preDLI_cox,
          timeVar = "intSCT2_5",
          parameterization = "both",
          iter.EM = 1000,
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
  tar_target(
    preDLI_CD8_jointModel_both,
    command = value(
      future({
        run_jointModel(
          long_obj = preDLI_CD8_long,
          surv_obj = preDLI_cox,
          timeVar = "intSCT2_5",
          parameterization = "both",
          iter.EM = 1000,
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
  )
)
