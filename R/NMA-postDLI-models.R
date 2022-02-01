postDLI_targets <- list(
  tar_target(
    postDLI_CD3_long,
    run_preDLI_longitudinal(
      cell_line = "CD3_abs_log",
      form_fixed = "ns(intSCT2_5_reset, 4) * earlylow_DLI + ATG", # add three-way? interaction later
      form_random = list(IDAA = pdDiag(~ ns(intSCT2_5_reset, 4))),
      dat = NMA_postDLI_datasets_CD3$long
    )
  ),
  tar_target(
    postDLI_CD4_long,
    run_preDLI_longitudinal(
      cell_line = "CD4_abs_log",
      form_fixed = "ns(intSCT2_5_reset, 3) * earlylow_DLI + ATG", # add three-way? interaction later
      form_random = list(IDAA = pdDiag(~ ns(intSCT2_5_reset, 3))),
      dat = NMA_postDLI_datasets_CD4$long
    )
  ),
  tar_target(
    postDLI_CD8_long,
    run_preDLI_longitudinal(
      cell_line = "CD8_abs_log",
      form_fixed = "ns(intSCT2_5_reset, 3) * earlylow_DLI + ATG",
      form_random = list(IDAA = pdDiag(~ ns(intSCT2_5_reset, 3))),
      dat = NMA_postDLI_datasets_CD8$long
    )
  ),
  tar_target(
    postDLI_CD3_cox,
    run_postDLI_cox(
      form = Surv(time, status) ~
        earlylow_DLI.1 + preds_subj.1  + # GVHD
        hirisk.2 + preds_subj.2 + # REL and NRF
        strata(trans),
      dat_wide = NMA_postDLI_datasets_CD3$wide
    )
  ),
  tar_target(
    postDLI_CD4_cox,
    run_postDLI_cox(
      form = Surv(time, status) ~
        earlylow_DLI.1 + preds_subj.1  + # GVHD
        hirisk.2 + preds_subj.2 + # REL and NRF
        strata(trans),
      dat_wide = NMA_postDLI_datasets_CD4$wide
    )
  ),
  tar_target(
    postDLI_CD8_cox,
    run_postDLI_cox(
      form = Surv(time, status) ~
        earlylow_DLI.1 + preds_subj.1  + # GVHD
        hirisk.2 + preds_subj.2 + # REL and NRF
        strata(trans),
      dat_wide = NMA_postDLI_datasets_CD8$wide
    )
  ),
  # Test first JM
  tar_target(
    postDLI_CD3_jointModel_both,
    jointModel(
      lmeObject = postDLI_CD3_long,
      survObject = postDLI_CD3_cox,
      CompRisk = TRUE,
      timeVar = "intSCT2_5_reset",
      method = "spline-PH-aGH",
      # derivForm = list(
      #   fixed = ~ 0 + dns(intSCT2_5_reset, 3) + dns(intSCT2_5_reset, 3):as.numeric(earlylow_DLI == "yes"),
      #   random = ~ 0 + dns(intSCT2_5_reset, 3),
      #   indFixed = c(2:4, 6:8),
      #   indRandom = c(2:4)
      # ),
      interFact = list(
        "value" = ~ strata(trans) - 1,
        "slope" = ~ strata(trans) - 1
      ),
      parameterization = "value",
      control = list("iter.EM" = 1000)
    )
  ),
  tar_target(
    postDLI_CD4_jointModel_both,
    jointModel(
      lmeObject = postDLI_CD4_long,
      survObject = postDLI_CD4_cox,
      CompRisk = TRUE,
      timeVar = "intSCT2_5_reset",
      method = "spline-PH-aGH",
      derivForm = list(
        fixed = ~ 0 + dns(intSCT2_5_reset, 3) + dns(intSCT2_5_reset, 3):as.numeric(earlylow_DLI == "yes"),
        random = ~ 0 + dns(intSCT2_5_reset, 3),
        indFixed = c(2:4, 6:8),
        indRandom = c(2:4)
      ),
      interFact = list(
        "value" = ~ strata(trans) - 1,
        "slope" = ~ strata(trans) - 1
      ),
      parameterization = "both",
      control = list("iter.EM" = 1000)
    )
  ),
  tar_target(
    postDLI_CD8_jointModel_both,
    jointModel(
      lmeObject = postDLI_CD8_long,
      survObject = postDLI_CD8_cox,
      CompRisk = TRUE,
      timeVar = "intSCT2_5_reset",
      method = "spline-PH-aGH",
      derivForm = list(
        fixed = ~ 0 + dns(intSCT2_5_reset, 3) + dns(intSCT2_5_reset, 3):as.numeric(earlylow_DLI == "yes"),
        random = ~ 0 + dns(intSCT2_5_reset, 3),
        indFixed = c(2:4, 6:8),
        indRandom = c(2:4)
      ),
      interFact = list(
        "value" = ~ strata(trans) - 1,
        "slope" = ~ strata(trans) - 1
      ),
      parameterization = "both",
      control = list("iter.EM" = 1000)
    )
  ),

  # For CD3 test things properly
  tar_target(
    postDLI_CD3_long_corr,
    run_preDLI_longitudinal(
      cell_line = "CD3_abs_log",
      form_fixed = "ns(intSCT2_5_reset, 3) * earlylow_DLI + ATG", # add three-way? interaction later
      form_random = ~ ns(intSCT2_5_reset, 3) | IDAA,
      dat = NMA_postDLI_datasets_CD3_corr$long
    )
  ),
  tar_target(
    postDLI_CD3_cox_corr,
    run_postDLI_cox(
      form = Surv(time, status) ~
        earlylow_DLI.1 + preds_subj.1  + # GVHD
        hirisk.2 + preds_subj.2 + # REL and NRF
        strata(trans),
      dat_wide = NMA_postDLI_datasets_CD3_corr$wide
    )
  ),
  tar_target(
    postDLI_CD3_jointModel_corr_both,
    jointModel(
      lmeObject = postDLI_CD3_long_corr,
      survObject = postDLI_CD3_cox_corr,
      CompRisk = TRUE,
      timeVar = "intSCT2_5_reset",
      method = "spline-PH-aGH",
      derivForm = list(
        fixed = ~ 0 + dns(intSCT2_5_reset, 3) + dns(intSCT2_5_reset, 3):as.numeric(earlylow_DLI == "yes"),
        random = ~ 0 + dns(intSCT2_5_reset, 3),
        indFixed = c(2:4, 6:8),
        indRandom = c(2:4)
      ),
      interFact = list(
        "value" = ~ strata(trans) - 1,
        "slope" = ~ strata(trans) - 1
      ),
      parameterization = "both",
      control = list("iter.EM" = 1000)
    )
  ),
  tar_target(
    postDLI_CD3_jointModel_corr_value,
    jointModel(
      lmeObject = postDLI_CD3_long_corr,
      survObject = postDLI_CD3_cox_corr,
      CompRisk = TRUE,
      timeVar = "intSCT2_5_reset",
      method = "spline-PH-aGH",
      interFact = list(
        "value" = ~ strata(trans) - 1
      ),
      parameterization = "value",
      control = list("iter.EM" = 1000)
    )
  )
)
