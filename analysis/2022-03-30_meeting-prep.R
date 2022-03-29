# Load objects ------------------------------------------------------------


mod_names <- c(
  paste0("preDLI_JM_value_corr_CD", c(3, 4, 8)),
  paste0("preDLI_JM_both_corr_CD", c(3, 4, 8)),
  paste0("postDLI_JM_corr_CD", c(3, 4, 8))
)
names(mod_names) <- mod_names
mods <- lapply(mod_names, tar_read_raw)

# Other objects
tar_load(c(NMA_preDLI_datasets, NMA_postDLI_datasets, reference_values))

# General theme
theme_set(theme_bw(base_size = 14))

# Load datasets
dat_long <- NMA_preDLI_datasets$long
dat_wide <- NMA_preDLI_datasets$wide
dat_long_postDLI <- NMA_postDLI_datasets$long
dat_wide_postDLI <- NMA_postDLI_datasets$long

# "Newdata" for predictions
newdat_preDLI <- expand.grid(
  "ATG" = levels(dat_long$ATG),
  "CMV_PD" = levels(dat_long$CMV_PD),
  "hirisk" = levels(dat_long$hirisk),
  "intSCT2_7" = seq(0, 6, by = 0.02)
)

newdat_postDLI <- expand.grid(
  "ATG" = levels(dat_long_postDLI$ATG),
  "CMV_PD" = levels(dat_long_postDLI$CMV_PD),
  "intDLI1" = seq(0, 3, by = 0.01)
)


# Helper functions --------------------------------------------------------



# Summary pre-DLI models --------------------------------------------------


# Marginals
mods_preDLI_value <- list(
  "CD3" = mods$preDLI_JM_value_corr_CD3,
  "CD4" = mods$preDLI_JM_value_corr_CD4,
  "CD8" = mods$preDLI_JM_value_corr_CD8
)

marg_preds_preDLI <- lapply(mods_preDLI_value, function(mod) {
  predict(
    mod,
    newdata = newdat_preDLI,
    type = "Marginal",
    idVar = "IDAA",
    returnData = TRUE,
    interval = "confidence"
  )
})


# Good plot with one
rbindlist(marg_preds_preDLI, idcol = "cell_line") |>
  ggplot(aes(x = intSCT2_7, y = pred, group = interaction(ATG, hirisk))) +
  geom_ribbon(aes(ymin = low, ymax = upp), fill = "gray", alpha = 0.5, col = NA) +
  geom_line(aes(linetype = hirisk, col = ATG), size = 1.5, alpha = 0.75) +
  # Add repel labels?
  facet_grid(facets = cell_line ~ CMV_PD) +
  labs(
    x = "Time since alloHCT (months)",
    y = "Cell counts",
    linetype = "ITT",
    col = "Donor type"
  ) +
  scale_y_continuous(
    breaks = log(c(5, 25, 100, 500, 1500)),
    labels = c(5, 25, 100, 500, 1500)
  ) +
  coord_cartesian(xlim = c(0, 6)) + # add a common ylim for all
  scale_color_manual(
    labels = c("Related", "Unrelated (+ATG)"),
    values = c("brown", "darkblue")
  ) +
  scale_linetype_manual(
    labels = c("No", "Yes (high risk)"),
    values = c("solid", "dotdash")
  ) +
  theme(legend.position = "bottom")



# Summary of post DLI models ----------------------------------------------

table(dat_wide_postDLI$sec_endpoint2_s)
table(dat_wide_postDLI$sec_endpoint2_s, dat_wide_postDLI$ATG)


mods_postDLI_value <- list(
  "CD3" = mods$postDLI_JM_corr_CD3,
  "CD4" = mods$postDLI_JM_corr_CD4,
  "CD8" = mods$postDLI_JM_corr_CD8
)

marg_preds_postDLI <- lapply(mods_postDLI_value, function(mod) {
  predict(
    mod,
    newdata = newdat_postDLI,
    type = "Marginal",
    idVar = "IDAA",
    returnData = TRUE,
    interval = "confidence"
  )
})


# Good plot with one
rbindlist(marg_preds_postDLI, idcol = "cell_line") |>
  ggplot(aes(x = intDLI1, y = pred, group = ATG)) +
  geom_ribbon(aes(ymin = low, ymax = upp), fill = "gray", alpha = 0.5, col = NA) +
  geom_line(aes(col = ATG), size = 1.5, alpha = 0.75) +
  # Add repel labels?
  facet_grid(facets = cell_line ~ CMV_PD) +
  labs(
    x = "Time since early/low DLI (months)",
    y = "Cell counts",
    linetype = "ITT",
    col = "Donor type"
  ) +
  scale_y_continuous(
    breaks = log(c(5, 25, 100, 500, 1500)),
    labels = c(5, 25, 100, 500, 1500)
  ) +
  coord_cartesian(xlim = c(0, 3)) + # add a common ylim for all
  scale_color_manual(
    labels = c("Related", "Unrelated (+ATG)"),
    values = c("brown", "darkblue")
  ) +
  theme(legend.position = "bottom")
