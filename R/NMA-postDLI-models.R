# Targets after DLI - update: in actual earlylow subset

# To-do:
# - correlated long sub-models
# - no more need for predicting from previous model
# - use same specification as pre-DLI (see new proposal)

postDLI_targets <- list(

  tar_map(
    values = cell_lines,
    tar_target(
      postDLI_long_corr,
      run_preDLI_longitudinal(
        cell_line = paste0(cell_line, "_abs_log"),
        form_fixed = "ns(intDLI1, 2) * ATG + CMV_PD", # remove cmvpd
        form_random = ~ ns(intDLI1, 2) | IDAA,
        dat = NMA_postDLI_datasets$long
      )
    )
  ),

  # Post-DLI cox: single model needed
  tar_target(
    postDLI_cox,
    run_postDLI_cox(
      Surv(time, status) ~ ATG.1 + strata(trans), # ATG.2 wanted but not enough events
      dat_wide = NMA_postDLI_datasets$wide
    )
  ),
  tar_target(postDLI_basehaz_knots, 2L), # 3L used pre-DLI, consider 1-2 here
  tar_target(
    postDLI_JM_corr_CD3,
    jointModel(
      lmeObject = postDLI_long_corr_CD3,
      survObject = postDLI_cox,
      CompRisk = TRUE,
      parameterization = "value",
      interFact = list("value" = ~ strata(trans) - 1),
      method = "spline-PH-aGH",
      timeVar = "intDLI1",
      lng.in.kn = postDLI_basehaz_knots,
      iter.EM = 1000,
      numeriDeriv = "cd",
      eps.Hes = 1e-04,
      verbose = TRUE
    )
  ),
  tar_target(
    postDLI_JM_corr_CD4,
    jointModel(
      lmeObject = postDLI_long_corr_CD4,
      survObject = postDLI_cox,
      CompRisk = TRUE,
      parameterization = "value",
      interFact = list("value" = ~ strata(trans) - 1),
      method = "spline-PH-aGH",
      timeVar = "intDLI1",
      lng.in.kn = postDLI_basehaz_knots,
      iter.EM = 1000,
      numeriDeriv = "cd",
      eps.Hes = 1e-04,
      verbose = TRUE
    )
  ),
  tar_target(
    postDLI_JM_corr_CD8,
    jointModel(
      lmeObject = postDLI_long_corr_CD8,
      survObject = postDLI_cox,
      CompRisk = TRUE,
      parameterization = "value",
      interFact = list("value" = ~ strata(trans) - 1),
      method = "spline-PH-aGH",
      timeVar = "intDLI1",
      lng.in.kn = postDLI_basehaz_knots,
      iter.EM = 1000,
      numeriDeriv = "cd",
      eps.Hes = 1e-04,
      verbose = TRUE
    )
  )
)


# GVHD sensitivity analysis
