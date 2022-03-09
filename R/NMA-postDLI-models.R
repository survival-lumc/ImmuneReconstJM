# Targets after DLI - update: in ITT subset

# To-do:
# - correlated long sub-models
# - no more need for predicting from previous model
# - use same specification as pre-DLI (see new proposal)

postDLI_targets <- list(
  #tar_target(postDLI_CD3_long, run_preDLI_longitudinal(...)),
  # For CD4 and CD8 targets..
  tar_target(
    postDLI_CD4_long_corr,
    run_preDLI_longitudinal(
      cell_line = "CD4_abs_log",
      form_fixed = "ns(intDLI1, 2) * ATG + CMV_PD", # remove cmvpd
      form_random = ~ ns(intDLI1, 2) | IDAA,
      #form_random = list(IDAA = pdDiag(~ 0 + ns(intDLI1, 2))),
      dat = NMA_postDLI_datasets$long
    )
  ),
  tar_target(
    postDLI_CD3_long_corr,
    run_preDLI_longitudinal(
      cell_line = "CD3_abs_log",
      form_fixed = "ns(intDLI1, 2) * ATG + CMV_PD", # remove cmvpd
      form_random = ~ ns(intDLI1, 2) | IDAA,
      #form_random = list(IDAA = pdDiag(~ 0 + ns(intDLI1, 2))),
      dat = NMA_postDLI_datasets$long
    )
  ),
  tar_target(
    postDLI_CD8_long_corr,
    run_preDLI_longitudinal(
      cell_line = "CD8_abs_log",
      form_fixed = "ns(intDLI1, 2) * ATG + CMV_PD", # remove cmvpd
      form_random = ~ ns(intDLI1, 2) | IDAA,
      #form_random = list(IDAA = pdDiag(~ 0 + ns(intDLI1, 2))),
      dat = NMA_postDLI_datasets$long
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
  tar_target(
    postDLI_CD4_JM_value_corr,
    jointModel(
      lmeObject = postDLI_CD4_long_corr,
      survObject = postDLI_cox,
      CompRisk = TRUE,
      parameterization = "value",
      interFact = list("value" = ~ strata(trans) - 1),
      method = "spline-PH-aGH",
      timeVar = "intDLI1",
      lng.in.kn = 2,
      iter.EM = 500#,
      #verbose = TRUE
    )
  )

  # Joint models with only current val (try CD3 with slope)
  # .. later with jmbayes2 tests

  # Also run tests with gvhd a week earlier..
)
