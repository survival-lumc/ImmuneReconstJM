list(
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
      iter.EM = 250L,
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
  )
)
