
# Old targets -------------------------------------------------------------


preDLI_targets_old <- list(

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
      iter.EM = 2000,
      #equal.strata.knots = FALSE, # Different knots per strata
      #lng.in.kn = 4, # baseline haz is cubic spline (ord = 4), try only 4 knots instead of 5
      interFact = list("value" = ~ strata(trans) - 1)
    )
  ),

  tar_target(
    preDLI_CD4_jointModel_value,
    jointModel(
      lmeObject = preDLI_CD4_long,
      survObject = preDLI_cox,
      CompRisk = TRUE,
      method = "spline-PH-aGH",
      timeVar = "intSCT2_5",
      parameterization = "value",
      iter.EM = 2000,
      #equal.strata.knots = FALSE, # Different knots per strata
      #lng.in.kn = 4, # baseline haz is cubic spline (ord = 4), try only 4 knots instead of 5
      interFact = list("value" = ~ strata(trans) - 1)
    )
  ),

  tar_target(
    preDLI_CD8_jointModel_value,
    jointModel(
      lmeObject = preDLI_CD8_long,
      survObject = preDLI_cox,
      CompRisk = TRUE,
      method = "spline-PH-aGH",
      timeVar = "intSCT2_5",
      parameterization = "value",
      iter.EM = 2000,
      #equal.strata.knots = FALSE, # Different knots per strata
      #lng.in.kn = 4, # baseline haz is cubic spline (ord = 4), try only 4 knots instead of 5
      interFact = list("value" = ~ strata(trans) - 1)
    )
  ),

  # -- Joint models using JM: current value + slope

  tar_target(
    preDLI_CD3_jointModel_both,
    update(
      preDLI_CD3_jointModel_value,
      lmeObject = preDLI_CD3_long,
      survObject = preDLI_cox,
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
  ),

  tar_target(
    preDLI_CD4_jointModel_both,
    update(
      preDLI_CD4_jointModel_value,
      lmeObject = preDLI_CD4_long,
      survObject = preDLI_cox,
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
  ),

  tar_target(
    preDLI_CD8_jointModel_both,
    update(
      preDLI_CD8_jointModel_value,
      lmeObject = preDLI_CD8_long,
      survObject = preDLI_cox,
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

  # -- Sensitivity analysis with GVHD one week earlier.. (need to edit DGMs)

  # -- Joint models with jmbayes2 (bayesian): current value

  # -- To try: combine relapse and nrf
)




# New targets -------------------------------------------------------------


preDLI_targets <- list(

  # -- Pre-DLI Cox model (no tdep DLI)
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

  # -- Longitudinal models - correlated reffs, fixed intercept
  tar_target(
    preDLI_CD3_long_corr,
    run_preDLI_longitudinal(
      cell_line = "CD3_abs_log",
      form_fixed = "ns(intSCT2_7, 3) * hirisk * ATG + CMV_PD",
      form_random = ~ 0 + ns(intSCT2_7, 3) | IDAA,
      dat = NMA_preDLI_datasets$long
    )
  ),
  tar_target(
    preDLI_CD4_long_corr,
    run_preDLI_longitudinal(
      cell_line = "CD4_abs_log",
      form_fixed = "ns(intSCT2_7, 3) * hirisk * ATG + CMV_PD",
      form_random = ~ 0 + ns(intSCT2_7, 3) | IDAA,
      dat = NMA_preDLI_datasets$long
    )
  ),
  tar_target(
    preDLI_CD8_long_corr,
    run_preDLI_longitudinal(
      cell_line = "CD8_abs_log",
      form_fixed = "ns(intSCT2_7, 3) * hirisk * ATG + CMV_PD",
      form_random = ~ 0 + ns(intSCT2_7, 3) | IDAA,
      dat = NMA_preDLI_datasets$long
    )
  ),

  # -- Longitudinal models - uncorrelated reffs, fixed intercept
  tar_target(
    preDLI_CD3_long_indep,
    run_preDLI_longitudinal(
      cell_line = "CD3_abs_log",
      form_fixed = "ns(intSCT2_7, 3) * hirisk * ATG + CMV_PD",
      form_random = list("IDAA" = pdDiag(~ 0 + ns(intSCT2_7, 3))),
      dat = NMA_preDLI_datasets$long
    )
  ),
  tar_target(
    preDLI_CD4_long_indep,
    run_preDLI_longitudinal(
      cell_line = "CD4_abs_log",
      form_fixed = "ns(intSCT2_7, 3) * hirisk * ATG + CMV_PD",
      form_random = list("IDAA" = pdDiag(~ 0 + ns(intSCT2_7, 3))),
      dat = NMA_preDLI_datasets$long
    )
  ),
  tar_target(
    preDLI_CD8_long_indep,
    run_preDLI_longitudinal(
      cell_line = "CD8_abs_log",
      form_fixed = "ns(intSCT2_7, 3) * hirisk * ATG + CMV_PD",
      form_random = list("IDAA" = pdDiag(~ 0 + ns(intSCT2_7, 3))),
      dat = NMA_preDLI_datasets$long
    )
  ),

  # -- Joint models - correlated reffs
  tar_target(
    preDLI_CD3_JM_value_corr,
    jointModel(
      lmeObject = preDLI_CD3_long_corr,
      survObject = preDLI_cox,
      CompRisk = TRUE,
      method = "spline-PH-aGH",
      timeVar = "intSCT2_7",
      parameterization = "value",
      #iter.EM = 500,
      interFact = list("value" = ~ strata(trans) - 1)
    )
  ),
  tar_target(
    preDLI_CD4_JM_value_corr,
    jointModel(
      lmeObject = preDLI_CD4_long_corr,
      survObject = preDLI_cox,
      CompRisk = TRUE,
      method = "spline-PH-aGH",
      timeVar = "intSCT2_7",
      parameterization = "value",
      #iter.EM = 500,
      interFact = list("value" = ~ strata(trans) - 1)
    )
  ),
  tar_target(
    preDLI_CD8_JM_value_corr,
    jointModel(
      lmeObject = preDLI_CD8_long_corr,
      survObject = preDLI_cox,
      CompRisk = TRUE,
      method = "spline-PH-aGH",
      timeVar = "intSCT2_7",
      parameterization = "value",
      #iter.EM = 500,
      interFact = list("value" = ~ strata(trans) - 1)
    )
  ),


  # -- Joint model - indep effs
  tar_target(
    preDLI_CD4_JM_value_indep,
    jointModel(
      lmeObject = preDLI_CD4_long_indep,
      survObject = preDLI_cox,
      CompRisk = TRUE,
      method = "spline-PH-aGH",
      timeVar = "intSCT2_7",
      parameterization = "value",
      interFact = list("value" = ~ strata(trans) - 1)
    )
  ),
  tar_target(
    preDLI_CD4_JM_both_indep,
    update(
      preDLI_CD4_JM_value_indep,
      lmeObject = preDLI_CD4_long_indep,
      survObject = preDLI_cox,
      parameterization = "both",
      interFact = list("value" = ~ strata(trans) - 1, "slope" = ~ strata(trans) - 1),
      derivForm = list(
        fixed = ~ 0 + dns(intSCT2_7, 3) +
          dns(intSCT2_7, 3):as.numeric(hirisk == "yes") +
          dns(intSCT2_7, 3):as.numeric(ATG == "ALT+ATG") +
          dns(intSCT2_7, 3):as.numeric(hirisk == "yes"):as.numeric(ATG == "ALT+ATG"),
        random = ~ 0 + dns(intSCT2_7, 3),
        indFixed = c(2:4, 8:13, 15:17),
        indRandom = c(1:3)
      )
    )
  ),
  tar_target(
    preDLI_CD8_JM_value_indep,
    jointModel(
      lmeObject = preDLI_CD8_long_indep,
      survObject = preDLI_cox,
      CompRisk = TRUE,
      method = "spline-PH-aGH",
      timeVar = "intSCT2_7",
      parameterization = "value",
      #iter.EM = 500,
      #lng.in.kn = 4, # baseline haz is cubic spline (ord = 4), try only 4 knots instead of 5
      interFact = list("value" = ~ strata(trans) - 1)
    )
  ),
  tar_target(
    preDLI_CD8_JM_both_indep,
    update(
      preDLI_CD8_JM_value_indep,
      lmeObject = preDLI_CD8_long_indep,
      survObject = preDLI_cox,
      parameterization = "both",
      interFact = list("value" = ~ strata(trans) - 1, "slope" = ~ strata(trans) - 1),
      derivForm = list(
        fixed = ~ 0 + dns(intSCT2_7, 3) +
          dns(intSCT2_7, 3):as.numeric(hirisk == "yes") +
          dns(intSCT2_7, 3):as.numeric(ATG == "ALT+ATG") +
          dns(intSCT2_7, 3):as.numeric(hirisk == "yes"):as.numeric(ATG == "ALT+ATG"),
        random = ~ 0 + dns(intSCT2_7, 3),
        indFixed = c(2:4, 8:13, 15:17),
        indRandom = c(1:3)
      )
    )
  ),
  tar_target(
    preDLI_CD3_JM_value_indep,
    jointModel(
      lmeObject = preDLI_CD3_long_indep,
      survObject = preDLI_cox,
      CompRisk = TRUE,
      method = "spline-PH-aGH",
      timeVar = "intSCT2_7",
      parameterization = "value",
      #iter.EM = 500,
      interFact = list("value" = ~ strata(trans) - 1)
    )
  )
)

