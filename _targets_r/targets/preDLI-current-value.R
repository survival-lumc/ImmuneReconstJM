list(
  # -- Pre-DLI Cox model
  tar_target(
    preDLI_cox,
    run_preDLI_cox(
      form = Surv(time, status) ~
        ATG.1 + hirisk.1 + # GVHD
        ATG.2 + hirisk.2 + # Relapse
        ATG.3 + # NRF other
        strata(trans),
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
      iter.qN = 1000L,
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
      iter.qN = 1000L,
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
      iter.qN = 1000L,
      lng.in.kn = preDLI_basehaz_knots,
      numeriDeriv = "cd",
      eps.Hes = 1e-04,
      verbose = TRUE
    )
  )
)
