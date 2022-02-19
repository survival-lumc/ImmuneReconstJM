preDLI_targets <- list(

  # -- Longitudinal submodels (correlated random effects)
  tar_target(
    preDLI_CD3_long,
    run_preDLI_longitudinal(
      cell_line = "CD3_abs_log",
      form_fixed = "ns(intSCT2_5, 3) * hirisk * ATG + CMV_PD",
      form_random = ~ ns(intSCT2_5, 3) | IDAA,
      dat = NMA_preDLI_datasets$long
    )
  ),
  tar_target(
    preDLI_CD4_long,
    run_preDLI_longitudinal(
      cell_line = "CD4_abs_log",
      form_fixed = "ns(intSCT2_5, 3) * hirisk * ATG + CMV_PD",
      form_random = ~ ns(intSCT2_5, 3) | IDAA,
      dat = NMA_preDLI_datasets$long
    )
  ),
  tar_target(
    preDLI_CD8_long,
    run_preDLI_longitudinal(
      cell_line = "CD8_abs_log",
      form_fixed = "ns(intSCT2_5, 3) * hirisk * ATG + CMV_PD",
      form_random = ~ ns(intSCT2_5, 3) | IDAA,
      dat = NMA_preDLI_datasets$long
    )
  ),

  # -- Survival submodel
  tar_target(
    preDLI_cox,
    run_preDLI_cox(
      form = Surv(time, status) ~
        ATG.1 + hirisk.1 + # GVHD
        ATG.2 + hirisk.2 + # Relapse
        ATG.3 + strata(trans), # NRF other
      dat_wide = NMA_preDLI_datasets$wide
    )
  ),

  # -- Joint models using JM (frequentist): current value

  tar_target(
    preDLI_CD3_jointModel_value,
    jointModel(
      lmeObject = preDLI_CD3_long,
      survObject = preDLI_cox,
      CompRisk = TRUE,
      method = "spline-PH-aGH",
      timeVar = "intSCT2_5",
      parameterization = "value",
      iter.EM = 1000,
      interFact = list("value" = ~ strata(trans) - 1)
    )
  ),

  # -- Joint models using JM: current value + slope

  tar_target(
    preDLI_CD3_jointModel_both,
    update(
      preDLI_CD3_jointModel_value,
      parameterization = "both",
      interFact = list("value" = ~ strata(trans) - 1, "slope" = ~ strata(trans) - 1),
      derivForm = list(
        fixed = ~ 0 + dns(intSCT2_5, 3) +
          dns(intSCT2_5, 3):as.numeric(hirisk == "yes") +
          dns(intSCT2_5, 3):as.numeric(ATG == "ALT+ATG") +
          dns(intSCT2_5, 3):as.numeric(hirisk == "yes"):as.numeric(ATG == "ALT+ATG"),
        random = ~ 0 + dns(intSCT2_5, 3),
        indFixed = c(2:4, 8:13, 15:17),
        indRandom = c(2:4)
      )
    )
  )

  # -- Joint models with jmbayes2 (bayesian): current value
)
