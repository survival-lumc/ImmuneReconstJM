list(
  tar_target(
    preDLI_long_corr_NK,
    run_preDLI_longitudinal(
      cell_line = "NK_abs_log",
      form_fixed = "ns(intSCT2_7, 3) * hirisk * ATG + CMV_PD",
      form_random = ~ 0 + ns(intSCT2_7, 3) | IDAA,
      dat = NMA_preDLI_datasets$long
    )
  ),
  tar_target(
    preDLI_JM_value_corr_NK,
    jointModel(
      lmeObject = preDLI_long_corr_NK,
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
    preDLI_tdc_dataset,
    make_tdc_dataset(
      dat_wide = NMA_preDLI_datasets$wide,
      dat_long = NMA_preDLI_datasets$long
    )
  )
)
