# preDLI models -----------------------------------------------------------


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

  # -- Longitudinal models - indep reffs, fixed intercept
  tar_map(
    values = cell_lines,
    tar_target(
      preDLI_long_indep,
      run_preDLI_longitudinal(
        cell_line = paste0(cell_line, "_abs_log"),
        form_fixed = "ns(intSCT2_7, 3) * hirisk * ATG + CMV_PD",
        form_random = list("IDAA" = pdDiag(~ 0 + ns(intSCT2_7, 3))),
        dat = NMA_preDLI_datasets$long
      )
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
      interFact = list("value" = ~ strata(trans) - 1)
    )
  ),
  # tar_target(
  #   preDLI_JM_value_corr_CD4,
  #   jointModel(
  #     lmeObject = preDLI_long_corr_CD4,
  #     survObject = preDLI_cox,
  #     CompRisk = TRUE,
  #     method = "spline-PH-aGH",
  #     timeVar = "intSCT2_7",
  #     parameterization = "value",
  #     interFact = list("value" = ~ strata(trans) - 1),
  #     iter.EM = 250
  #   )
  # ),
  tar_target(
    preDLI_JM_value_corr_CD8,
    jointModel(
      lmeObject = preDLI_long_corr_CD8,
      survObject = preDLI_cox,
      CompRisk = TRUE,
      method = "spline-PH-aGH",
      timeVar = "intSCT2_7",
      parameterization = "value",
      interFact = list("value" = ~ strata(trans) - 1)
    )
  ),

  # -- Joint models, independent reffs
  tar_target(
    preDLI_JM_value_indep_CD3,
    jointModel(
      lmeObject = preDLI_long_indep_CD3,
      survObject = preDLI_cox,
      CompRisk = TRUE,
      method = "spline-PH-aGH",
      timeVar = "intSCT2_7",
      parameterization = "value",
      interFact = list("value" = ~ strata(trans) - 1)
    )
  ),
  tar_target(
    preDLI_JM_value_indep_CD4,
    jointModel(
      lmeObject = preDLI_long_indep_CD4,
      survObject = preDLI_cox,
      CompRisk = TRUE,
      method = "spline-PH-aGH",
      timeVar = "intSCT2_7",
      parameterization = "value",
      interFact = list("value" = ~ strata(trans) - 1)
    )
  ),
  tar_target(
    preDLI_JM_value_indep_CD8,
    jointModel(
      lmeObject = preDLI_long_indep_CD8,
      survObject = preDLI_cox,
      CompRisk = TRUE,
      method = "spline-PH-aGH",
      timeVar = "intSCT2_7",
      parameterization = "value",
      interFact = list("value" = ~ strata(trans) - 1)
    )
  ),

  # -- Joint models with value + slope

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

  # -- Both, independent reffs
  tar_target(
    preDLI_JM_both_indep_CD3,
    update(
      preDLI_JM_value_indep_CD3,
      lmeObject = preDLI_long_indep_CD3,
      survObject = preDLI_cox,
      parameterization = "both",
      interFact = list("value" = ~ strata(trans) - 1, "slope" = ~ strata(trans) - 1),
      derivForm = derivForm_preDLI,
      iter.EM = 1000
    )
  ),
  tar_target(
    preDLI_JM_both_indep_CD4,
    update(
      preDLI_JM_value_indep_CD4,
      lmeObject = preDLI_long_indep_CD4,
      survObject = preDLI_cox,
      parameterization = "both",
      interFact = list("value" = ~ strata(trans) - 1, "slope" = ~ strata(trans) - 1),
      derivForm = derivForm_preDLI,
      iter.EM = 1000
    )
  ),
  tar_target(
    preDLI_JM_both_indep_CD8,
    update(
      preDLI_JM_value_indep_CD8,
      lmeObject = preDLI_long_indep_CD8,
      survObject = preDLI_cox,
      parameterization = "both",
      interFact = list("value" = ~ strata(trans) - 1, "slope" = ~ strata(trans) - 1),
      derivForm = derivForm_preDLI,
      iter.EM = 1000
    )
  )#,

  # -- Both, corr reffs
  # tar_target(
  #   preDLI_JM_both_corr_CD3,
  #   update(
  #     preDLI_JM_value_corr_CD3,
  #     lmeObject = preDLI_long_corr_CD3,
  #     survObject = preDLI_cox,
  #     parameterization = "both",
  #     interFact = list("value" = ~ strata(trans) - 1, "slope" = ~ strata(trans) - 1),
  #     derivForm = derivForm_preDLI,
  #     iter.EM = 500
  #   )
  # ),
  # tar_target(
  #   preDLI_JM_both_corr_CD4,
  #   update(
  #     preDLI_JM_value_corr_CD4,
  #     lmeObject = preDLI_long_corr_CD4,
  #     survObject = preDLI_cox,
  #     parameterization = "both",
  #     interFact = list("value" = ~ strata(trans) - 1, "slope" = ~ strata(trans) - 1),
  #     derivForm = derivForm_preDLI,
  #     iter.EM = 500
  #   )
  # ),
  # tar_target(
  #   preDLI_JM_both_corr_CD8,
  #   update(
  #     preDLI_JM_value_corr_CD8,
  #     lmeObject = preDLI_long_corr_CD8,
  #     survObject = preDLI_cox,
  #     parameterization = "both",
  #     interFact = list("value" = ~ strata(trans) - 1, "slope" = ~ strata(trans) - 1),
  #     derivForm = derivForm_preDLI,
  #     iter.EM = 500
  #   )
  # )

)



# GVHD sensitivity analysis -----------------------------------------------

# new targets list
