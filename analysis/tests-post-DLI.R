# Post-DLI, raw plots -----------------------------------------------------

tar_load(c(NMA_postDLI_datasets, postDLI_JM_corr_CD3))

dat_long <- NMA_postDLI_datasets$long
dat_wide <- NMA_postDLI_datasets$wide
#ids_oneobs <- dat_long[, .N, by = IDAA][N == 1][["IDAA"]]


dat_long |>
  melt.data.table(
    measure.vars = patterns("*_log$"),
    variable.name = "cell_line",
    value.name = "count"
  ) |>
  ggplot(aes(intDLI1, count, col = cell_line)) +
  geom_point(aes(shape = cell_line), size = 4, alpha = 0.75) +
  facet_wrap(sec_endpoint2_s ~ IDAA) +
  scale_y_continuous(
    breaks = log(c(5, 25, 100, 500, 1500)),
    labels = c(5, 25, 100, 500, 1500)
  ) +
  theme_bw(base_size = 14) +
  labs(y = "Count", x = "Time since early/low DLI (months)")


dat_long |>
  melt.data.table(
    measure.vars = patterns("*_log$"),
    variable.name = "cell_line",
    value.name = "count"
  ) |>
  ggplot(aes(intDLI1, count, col = ATG, group = IDAA)) +
  geom_line(size = 1) +
  scale_y_continuous(
    breaks = log(c(5, 25, 100, 500, 1500)),
    labels = c(5, 25, 100, 500, 1500)
  ) +
  facet_grid(sec_endpoint2_s ~ cell_line) +
  theme_bw(base_size = 14) +
  labs(y = "Count", x = "Time since early/low DLI (months)")



# Try random effects plot -------------------------------------------------


GGally::ggpairs(
  data = as.data.frame(ranef(postDLI_JM_corr_CD3))
) +
  theme_bw(base_size = 14)

# Post-DLI CD4 ------------------------------------------------------------

# Note postDLI cox omits ATG.2 due to event numbers
# Also postDLI model has less baseline hazard knots

summary(postDLI_JM_corr_CD3)

# Predict marginal
newdat_postDLI <- expand.grid(
  "ATG" = levels(dat_long$ATG),
  "CMV_PD" = levels(dat_long$CMV_PD),
  "intDLI1" = seq(0, 3, by = 0.01)
)

predict(
  postDLI_JM_corr_CD3,
  newdata = newdat_postDLI,
  interval = "confidence",
  type = "Marginal",
  returnData = TRUE
) |>
  ggplot(aes(intDLI1, pred, col = ATG, group = ATG)) +
  geom_ribbon(
    aes(ymin = low, ymax = upp),
    fill = "gray",
    alpha = 0.5,
    col = NA
  ) +
  geom_line(size = 1.25) +
  scale_y_continuous(
    breaks = log(c(5, 25, 100, 500, 1500)),
    labels = c(5, 25, 100, 500, 1500)
  ) +
  coord_cartesian(ylim = c(0, 8)) +
  facet_wrap(~ CMV_PD) +
  theme_bw(base_size = 14) +
  labs(y = "Count", x = "Time since early/low DLI (months)")


# Try individual predictions
dat_long[, "fitted_ind" := fitted(
  postDLI_JM_corr_CD3,
  process = "Longitudinal",
  type = "Subject"
)]


ggplot(dat_long, aes(intDLI1, CD3_abs_log)) +
  geom_point(
    size = 3.5, pch = 21,
    alpha = 0.8,
    col = "#359fda",
    fill = colorspace::lighten("#359fda", amount = 0.3)
  ) +
  geom_line(aes(y = fitted_ind), size = 1.25) +
  geom_point(aes(y = fitted_ind), size = 1.5) +
  facet_wrap(~ IDAA) +
  scale_y_continuous(
    breaks = log(c(5, 25, 100, 500, 1500)),
    labels = c(5, 25, 100, 500, 1500)
  ) +
  coord_cartesian(ylim = c(0, 8)) +
  theme_bw(base_size = 14) +
  labs(y = "Count", x = "Time since early/low DLI (months)")



endpoints_df <- dat_long[order(intDLI1), .SD[.N], by = IDAA]
dat_long[, .(.N, "endp" = unique(sec_endpoint2_s)), by = IDAA]

ggplot(dat_long, aes(intDLI1, fitted_JM, group = IDAA)) +
  geom_line(alpha = 0.75, size = 1, aes(col = IDAA)) +
  facet_grid(sec_endpoint2_s ~ ATG) +
  geom_point(data = endpoints_df, aes(sec_endpoint2, CD4_abs_log, shape = sec_endpoint2_s)) +
  theme(legend.position = "none")

