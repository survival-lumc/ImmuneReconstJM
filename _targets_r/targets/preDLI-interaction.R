list(
  tar_target(
    preDLI_JM_value_corr_CD4_inter,
    update(
      preDLI_JM_value_corr_CD4,
      lmeObject = preDLI_long_corr_CD4,
      survObject = tweak_preDLI_modmat(preDLI_cox),
      parameterization = "value",
      lng.in.kn = preDLI_basehaz_knots,
      interFact = list(
        "value" = ~ trans1 + trans1:ATG.1 + trans2 + trans2:ATG.2 +
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
        "value" = ~ trans1 + trans1:ATG.1 + trans2 + trans2:ATG.2 + trans3 - 1
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
        "value" = ~ trans1 + trans1:ATG.1 + trans2 + trans2:ATG.2 + trans3 - 1
      ),
      iter.EM = 250L,
      verbose = TRUE
    )
  )
)
