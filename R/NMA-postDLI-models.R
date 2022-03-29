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
  tar_map(
    values = cell_lines,
    tar_target(
      postDLI_long_indep,
      run_preDLI_longitudinal(
        cell_line = paste0(cell_line, "_abs_log"),
        form_fixed = "ns(intDLI1, 2) * ATG + CMV_PD", # remove cmvpd
        form_random = list(IDAA = pdDiag(~ ns(intDLI1, 2))),
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
  tar_target(postDLI_basehaz_knots, 2L), # 2L is probs best
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
      iter.EM = 500,
      iter.qN = 500,
      verbose = TRUE
    )
  )#,
  # tar_target(
  #   postDLI_CD8_JM_value_corr,
  #   jointModel(
  #     lmeObject = postDLI_CD8_long_corr,
  #     survObject = postDLI_cox,
  #     CompRisk = TRUE,
  #     parameterization = "value",
  #     interFact = list("value" = ~ strata(trans) - 1),
  #     method = "spline-PH-aGH",
  #     timeVar = "intDLI1",
  #     lng.in.kn = postDLI_basehaz_knots,
  #     iter.EM = 500#,
  #     #verbose = TRUE
  #   )
  # )

  # Joint models with only current val
  # .. later with jmbayes2 tests

  # Also run tests with gvhd a week earlier..
)


# GVHD sensitivity analysis
