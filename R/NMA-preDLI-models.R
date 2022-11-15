# preDLI models -----------------------------------------------------------


preDLI_targets <- list(

  # -- Pre-DLI Cox model
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
  tar_map(
    values = cell_lines,
    tar_target(
      preDLI_long_corr,
      run_preDLI_longitudinal(
        cell_line = paste0(cell_line, "_abs_log"),
        form_fixed = "ns(intSCT2_7, 3) * hirisk * ATG + CMV_PD",
        form_random = ~ 0 + ns(intSCT2_7, 3) | IDAA,
        dat = NMA_preDLI_datasets$long
      )
    )
  ),

  # Not many events for we lower default flexibility of baseline hazard (from 5 to 3 internal knots)
  tar_target(preDLI_basehaz_knots, 3L),

  # -- Joint models, correlated reffs
  tar_target(
    preDLI_JM_value_corr_CD3,
    jointModel(
      lmeObject = preDLI_long_corr_CD3,
      survObject = preDLI_cox,
      CompRisk = TRUE,
      method = "spline-PH-aGH",
      timeVar = "intSCT2_7",
      parameterization = "value",
      interFact = list("value" = ~ strata(trans) - 1),
      iter.qN = 1000,
      lng.in.kn = preDLI_basehaz_knots,
      numeriDeriv = "cd",
      eps.Hes = 1e-04
    )
  ),
  tar_target(
    preDLI_JM_value_corr_CD4,
    jointModel(
      lmeObject = preDLI_long_corr_CD4,
      survObject = preDLI_cox,
      CompRisk = TRUE,
      method = "spline-PH-aGH",
      timeVar = "intSCT2_7",
      parameterization = "value",
      interFact = list("value" = ~ strata(trans) - 1),
      iter.qN = 1000,
      lng.in.kn = preDLI_basehaz_knots,
      numeriDeriv = "cd",
      eps.Hes = 1e-04
    )
  ),
  tar_target(
    preDLI_JM_value_corr_CD8,
    jointModel(
      lmeObject = preDLI_long_corr_CD8,
      survObject = preDLI_cox,
      CompRisk = TRUE,
      method = "spline-PH-aGH",
      timeVar = "intSCT2_7",
      parameterization = "value",
      interFact = list("value" = ~ strata(trans) - 1),
      iter.qN = 1000,
      lng.in.kn = preDLI_basehaz_knots,
      numeriDeriv = "cd",
      eps.Hes = 1e-04,
      verbose = TRUE
    )
  ),

  # -- Both current value + slope, corr reffs

  # Store the deriv form (since same for all mods)
  tar_target(
    derivForm_preDLI,
    list(
      fixed = ~ 0 + dns(intSCT2_7, 3) +
        dns(intSCT2_7, 3):as.numeric(hirisk == "yes") +
        dns(intSCT2_7, 3):as.numeric(ATG == "ALT+ATG") +
        dns(intSCT2_7, 3):as.numeric(hirisk == "yes"):as.numeric(ATG == "ALT+ATG"),
      random = ~ 0 + dns(intSCT2_7, 3),
      indFixed = c(2:4, 8:13, 15:17),
      indRandom = c(1:3)
    )
  ),
  tar_target(
    preDLI_JM_both_corr_CD3,
    update(
      preDLI_JM_value_corr_CD3,
      lmeObject = preDLI_long_corr_CD3,
      survObject = preDLI_cox,
      parameterization = "both",
      interFact = list("value" = ~ strata(trans) - 1, "slope" = ~ strata(trans) - 1),
      derivForm = derivForm_preDLI,
      lng.in.kn = preDLI_basehaz_knots
    )
  ),
  # This one gets stuck local minima.. update iter.EM to 250
  tar_target(
    preDLI_JM_both_corr_CD4,
    update(
      preDLI_JM_value_corr_CD4,
      lmeObject = preDLI_long_corr_CD4,
      survObject = preDLI_cox,
      parameterization = "both",
      interFact = list("value" = ~ strata(trans) - 1, "slope" = ~ strata(trans) - 1),
      derivForm = derivForm_preDLI,
      lng.in.kn = preDLI_basehaz_knots,
      iter.EM = 250,
      verbose = TRUE
    )
  ),
  tar_target(
    preDLI_JM_both_corr_CD8,
    update(
      preDLI_JM_value_corr_CD8,
      lmeObject = preDLI_long_corr_CD8,
      survObject = preDLI_cox,
      parameterization = "both",
      interFact = list("value" = ~ strata(trans) - 1, "slope" = ~ strata(trans) - 1),
      derivForm = derivForm_preDLI,
      lng.in.kn = preDLI_basehaz_knots
    )
  ),

  # -- Try interaction analyses of current value with UD(+ATG)

  tar_target(
    preDLI_JM_value_corr_CD4_inter,
    update(
      preDLI_JM_value_corr_CD4,
      lmeObject = preDLI_long_corr_CD4,
      survObject = tweak_preDLI_modmat(preDLI_cox),
      parameterization = "value",
      lng.in.kn = preDLI_basehaz_knots,
      interFact = list(
        "value" = ~ trans1 + trans1:ATG.1 +
          trans2 + trans2:ATG.2 +
          # omit trans3 interaction since not interesting, could also do for slope models
          # i.e. include slope only relapse and gvhd
          trans3 - 1
      ),
      verbose = TRUE
    )
  ),
  tar_target(
    preDLI_JM_value_corr_CD8_inter,
    update(
      preDLI_JM_value_corr_CD8,
      lmeObject = preDLI_long_corr_CD8,
      survObject = tweak_preDLI_modmat(preDLI_cox),
      parameterization = "value",
      lng.in.kn = preDLI_basehaz_knots,
      interFact = list(
        "value" = ~ trans1 + trans1:ATG.1 +
          trans2 + trans2:ATG.2 + trans3 - 1
      ),
      verbose = TRUE
    )
  ),
  tar_target(
    preDLI_JM_value_corr_CD3_inter,
    update(
      preDLI_JM_value_corr_CD3,
      lmeObject = preDLI_long_corr_CD3,
      survObject = tweak_preDLI_modmat(preDLI_cox),
      parameterization = "value",
      lng.in.kn = preDLI_basehaz_knots,
      interFact = list(
        "value" = ~ trans1 + trans1:ATG.1 +
          trans2 + trans2:ATG.2 + trans3 - 1
      ),
      iter.EM = 250,
      verbose = TRUE
    )
  )
)

