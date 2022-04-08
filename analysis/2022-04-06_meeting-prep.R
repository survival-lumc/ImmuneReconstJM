tar_load(
  c(
    NMA_preDLI_datasets,
    NMA_postDLI_datasets,
    preDLI_cox,
    preDLI_long_corr_CD4,

    # Pre-DLI interaction models
    preDLI_JM_value_corr_CD3_inter,
    preDLI_JM_value_corr_CD4_inter,
    preDLI_JM_value_corr_CD8_inter
  )
)



# First look at pre-DLI contrasts -----------------------------------------


UD_preDLI_long <- NMA_preDLI_datasets$long[ATG == "UD(+ATG)"]
UD_preDLI_wide <- NMA_preDLI_datasets$wide[ATG == "UD(+ATG)"]
RD_preDLI_long <- NMA_preDLI_datasets$long[ATG == "UD"]
RD_preDLI_wide <- NMA_preDLI_datasets$wide[ATG == "UD"]

# All UD hirisk patients get gvhd - a cox model for the subset simplifies to
# just hirisk at predictor for relapse (but not events probably with current value + interaction)
table(UD_preDLI_wide$endpoint7_s)
table(UD_preDLI_wide$endpoint7_s, UD_preDLI_wide$hirisk)

# Too few events generally in RD
table(RD_preDLI_wide$endpoint7_s)
table(RD_preDLI_wide$endpoint7_s, RD_preDLI_wide$hirisk)


# Post-DLI contrasts ------------------------------------------------------


UD_postDLI_long <- NMA_postDLI_datasets$long[ATG == "UD(+ATG)"]
UD_postDLI_wide <- NMA_postDLI_datasets$wide[ATG == "UD(+ATG)"]
RD_postDLI_long <- NMA_postDLI_datasets$long[ATG == "UD"]
RD_postDLI_wide <- NMA_postDLI_datasets$wide[ATG == "UD"]

# Cox model before only had UD as predictor, so now it would only have current value
# This analysis is not possible post DLI.
table(UD_postDLI_wide$sec_endpoint2_s)
table(RD_postDLI_wide$sec_endpoint2_s)


# Look at pre-DLI interaction models --------------------------------------


# Perhaps more efficient to use interaction parametrisation (just for current value)

# Look at CD4 first - interesting result for relapse, but loose power for gvhd
# For relapse: effect current value only for RD, not for UD
# GVHD: (despite low power) - effects same magnitude in both
round(summary(tar_read(preDLI_JM_value_corr_CD4))$`CoefTable-Event`[1:8, ], 3)
round(summary(tar_read(preDLI_JM_value_corr_CD4_inter))$`CoefTable-Event`[1:11, ], 3)

# CD8
round(summary(tar_read(preDLI_JM_value_corr_CD8))$`CoefTable-Event`[1:8, ], 3)
round(summary(tar_read(preDLI_JM_value_corr_CD8_inter))$`CoefTable-Event`[1:11, ], 3)

# CD3
round(summary(tar_read(preDLI_JM_value_corr_CD3))$`CoefTable-Event`[1:8, ], 3)
round(summary(tar_read(preDLI_JM_value_corr_CD3_inter))$`CoefTable-Event`[1:11, ], 3)

# To do: distribution of current values -----------------------------------


# ...
