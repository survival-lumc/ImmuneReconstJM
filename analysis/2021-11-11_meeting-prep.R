# Focus on NMA (extended cohort)

# Load object
tar_load(
  c(
    JM_CD3_allDLI,
    JM_CD4_allDLI,
    JM_CD8_allDLI,
    cox_all_dli,
    long_submodels,
    datasets,
    dat_merged,
    reference_values
  )
)

# Datasets
dat_wide <- datasets$wide
dat_long <- datasets$long

dat_long[, cr_ind := factor(
  x = ifelse(endpoint6_s != "cell_interv", as.character(endpoint6_s), "cens"),
  levels = c("cens", "REL", "NRF_gvhd", "NRF_other")
)]

# Data to predict
newdat <- expand.grid(
  "ATG" = levels(dat_wide$ATG),
  "CMV_PD" = levels(dat_wide$CMV_PD),
  "intSCT2_5" = seq(0.1, 24, by = 0.1)
)

# ggplot labels
theme_set(theme_bw(base_size = 14))



# ID samples --------------------------------------------------------------

set.seed(2018)
id_samps <- sample(dat_long$IDAA, size = 24, replace = FALSE)

ids_abs <- dat_long[IDAA %in% id_samps] |>
  ggplot(aes(intSCT2_5, exp(CD3_abs_log))) +
  geom_point(
    size = 3, pch = 21,
    alpha = 0.8,
    col = "#359fda",
    fill = colorspace::lighten("#359fda", amount = 0.3)
  ) +
  #geom_line(aes(group = cell_line)) +
  facet_wrap(~ IDAA, ncol = 6) +
  # scale_y_continuous(
  #   breaks = log(c(1, 5, 25, 100, 500, 2000)),
  #   labels = c(1, 5, 25, 100, 500, 2000)
  # ) +
  labs(x = "Time since HSCT (months)", y = "CD3 cell counts")

ids_log_noscale <- dat_long[IDAA %in% id_samps] |>
  ggplot(aes(intSCT2_5, CD3_abs_log)) +
  geom_point(
    size = 3, pch = 21,
    alpha = 0.8,
    col = "#359fda",
    fill = colorspace::lighten("#359fda", amount = 0.3)
  ) +
  #geom_line(aes(group = cell_line)) +
  facet_wrap(~ IDAA, ncol = 6) +
  # scale_y_continuous(
  #   breaks = log(c(1, 5, 25, 100, 500, 2000)),
  #   labels = c(1, 5, 25, 100, 500, 2000)
  # ) +
  labs(x = "Time since HSCT (months)", y = "Log CD3 cell counts")

ids_log <- dat_long[IDAA %in% id_samps] |>
  ggplot(aes(intSCT2_5, CD3_abs_log)) +
  geom_point(
    size = 3, pch = 21,
    alpha = 0.8,
    col = "#359fda",
    fill = colorspace::lighten("#359fda", amount = 0.3)
  ) +
  #geom_line(aes(group = cell_line)) +
  facet_wrap(~ IDAA, ncol = 6) +
  scale_y_continuous(
    breaks = log(c(1, 5, 25, 100, 500, 2000)),
    labels = c(1, 5, 25, 100, 500, 2000)
  ) +
  labs(x = "Time since HSCT (months)", y = "CD3 cell counts")


ids_log_step <- dat_long[IDAA %in% id_samps] |>
  ggplot(aes(intSCT2_5, CD3_abs_log)) +
  geom_point(
    size = 3, pch = 21,
    alpha = 0.8,
    col = "#359fda",
    fill = colorspace::lighten("#359fda", amount = 0.3)
  ) +
  geom_step(size = 1) +
  facet_wrap(~ IDAA, ncol = 6) +
  scale_y_continuous(
    breaks = log(c(1, 5, 25, 100, 500, 2000)),
    labels = c(1, 5, 25, 100, 500, 2000)
  ) +
  labs(x = "Time since HSCT (months)", y = "CD3 cell counts")

ids_flex_ls <- lapply(c(1, 4), function(df) {
  dat_long[IDAA %in% id_samps] |>
    ggplot(aes(intSCT2_5, CD3_abs_log)) +
    geom_point(
      size = 3, pch = 21,
      alpha = 0.8,
      col = "#359fda",
      fill = colorspace::lighten("#359fda", amount = 0.3)
    ) +
    geom_smooth(
      method = "lm",
      se = FALSE,
      formula = y ~ ns(x, df = df),
      col = "black",
      size = 1
    ) +
    facet_wrap(~ IDAA, ncol = 6) +
    scale_y_continuous(
      breaks = log(c(1, 5, 25, 100, 500, 2000)),
      labels = c(1, 5, 25, 100, 500, 2000)
    ) +
    labs(x = "Time since HSCT (months)", y = "CD3 cell counts")
})


# Export all
ggsave(
  filename = "C:/Users/efbonneville/Downloads/immume-rec-talk/ids_abs.png",
  plot = ids_abs,
  dpi = 300,
  width = 16,
  height = 10,
  units = "in",
  scale = 0.7
)
ggsave(
  filename = "C:/Users/efbonneville/Downloads/immume-rec-talk/ids_log_noscale.png",
  plot = ids_log_noscale,
  dpi = 300,
  width = 16,
  height = 10,
  units = "in",
  scale = 0.7
)
ggsave(
  filename = "C:/Users/efbonneville/Downloads/immume-rec-talk/ids_log.png",
  plot = ids_log,
  dpi = 300,
  width = 16,
  height = 10,
  units = "in",
  scale = 0.7
)
ggsave(
  filename = "C:/Users/efbonneville/Downloads/immume-rec-talk/ids_log_step.png",
  plot = ids_log_step,
  dpi = 300,
  width = 16,
  height = 10,
  units = "in",
  scale = 0.7
)
ggsave(
  filename = "C:/Users/efbonneville/Downloads/immume-rec-talk/ids_log_lin.png",
  plot = ids_flex_ls[[1]],
  dpi = 300,
  width = 16,
  height = 10,
  units = "in",
  scale = 0.7
)
ggsave(
  filename = "C:/Users/efbonneville/Downloads/immume-rec-talk/ids_log_flex.png",
  plot = ids_flex_ls[[2]],
  dpi = 300,
  width = 16,
  height = 10,
  units = "in",
  scale = 0.7
)


# Average CD3 -------------------------------------------------------------


avg_cd3 <- dat_long |>
  ggplot(aes(intSCT2_5, CD3_abs_log)) +
  geom_line(aes(group = IDAA), col = "gray") +
  geom_smooth(
    method = "lm",
    se = FALSE,
    formula = y ~ ns(x, df = 4),
    col = "black",
    size = 1
  ) +
  scale_y_continuous(
    breaks = log(c(1, 5, 25, 100, 500, 2000)),
    labels = c(1, 5, 25, 100, 500, 2000)
  ) +
  labs(x = "Time since HSCT (months)", y = "CD3 cell counts")

ggsave(
  filename = "C:/Users/efbonneville/Downloads/immume-rec-talk/avg_cd3.png",
  plot = avg_cd3,
  dpi = 300,
  width = 10,
  height = 10,
  units = "in",
  scale = 0.6
)


# CD3 by endpoint ---------------------------------------------------------


CD3_per_endpoint <- dat_long |>
  ggplot(aes(intSCT2_5, CD3_abs_log)) +
  geom_line(aes(group = IDAA), alpha = 0.7, col = "gray") +
  geom_smooth(
    method = "lm",
    se = FALSE,
    aes(group = ATG, col = ATG),
    size = 1.25,
    formula = y ~ ns(x, df = 4)
  ) +
  facet_wrap(
    . ~ cr_ind,
    ncol = 4,
    labeller = as_labeller(
      x = c(
        "cens" = "Event-free",
        "REL" = "Relapse",
        "NRF_gvhd" = "GVHD",
        "NRF_other" = "Other NRF failures"
      )
    )
  ) +
  scale_y_continuous(
    breaks = log(c(1, 5, 25, 100, 500, 2000)),
    labels = c(1, 5, 25, 100, 500, 2000)
  ) +
  labs(x = "Time since HSCT (months)", y = "CD3 cell counts") +
  theme(legend.position = "top")

ggsave(
  filename = "C:/Users/efbonneville/Downloads/immume-rec-talk/CD3_per_endpoint.png",
  plot = CD3_per_endpoint,
  dpi = 300,
  width = 12,
  height = 8,
  units = "in",
  scale = 0.6
)

# Descriptives ------------------------------------------------------------


NMA_melted <- melt.data.table(
  data = dat_long,
  id.vars = c("IDAA", "ATG", "CMV_PD", "intSCT2_5", "cr_ind"),
  measure.vars = paste0(c("CD4", "CD8", "CD19", "CD3", "NK"), "_abs_log"),
  variable.name = "cell_line",
  value.name = "cell_count"
)


# Overall (by cell line)
NMA_melted |>
  ggplot(aes(intSCT2_5, cell_count)) +
  geom_line(aes(group = IDAA), alpha = 0.5, col = "lightgray") +
  geom_smooth(
    method = "loess",
    se = FALSE,
    aes(group = ATG, col = ATG),
    size = 1.25,
    formula = y ~ x
  ) +
  facet_grid(. ~ cell_line) +
  scale_y_continuous(
    breaks = log(c(1, 5, 25, 100, 500, 2000)),
    labels = c(1, 5, 25, 100, 500, 2000)
  ) +
  labs(x = "Time since HSCT (months)", y = "Cell counts")

# Per endpoint
nma_allcells <- NMA_melted |>
  ggplot(aes(intSCT2_5, cell_count)) +
  geom_line(aes(group = IDAA), alpha = 0.8, col = "lightgray") +
  geom_smooth(
    method = "lm",
    se = FALSE,
    aes(group = ATG, col = ATG),
    size = 1.25,
    formula = y ~ ns(x, df = 4)
  ) +
  #facet_grid(endpoint6_s ~ cell_line) +
  facet_grid(
    cell_line ~ cr_ind,
    labeller = labeller(
      cr_ind = c(
        "cens" = "Event-free",
        "REL" = "Relapse",
        "NRF_gvhd" = "GVHD",
        "NRF_other" = "Other NRF failures"
      ),
      cell_line = c(
        "CD3_abs_log" = "CD3",
        "CD4_abs_log" = "CD4",
        "CD8_abs_log" = "CD8",
        "CD19_abs_log" = "CD19",
        "NK_abs_log" = "NK"
      )
    )
  ) +
  scale_y_continuous(
    breaks = log(c(1, 5, 25, 100, 500, 2000)),
    labels = c(1, 5, 25, 100, 500, 2000)
  ) +
  labs(x = "Time since HSCT (months)", y = "Cell counts")

ggsave(
  filename = "C:/Users/efbonneville/Downloads/immume-rec-talk/nma_allcells.png",
  plot = nma_allcells,
  dpi = 300,
  width = 14,
  height = 8,
  units = "in",
  scale = 0.8
)

# Mixed models on their own -----------------------------------------------


cbind.data.frame(
  newdat,
  "preds" = predict(
    long_submodels$CD8_abs_log,
    newdata = newdat,
    level = 0L
  )
) |>
  ggplot(aes(intSCT2_5, preds, col = ATG)) +
  geom_line(size = 1.25) +
  facet_wrap(~ CMV_PD)



# CD3 model ---------------------------------------------------------------


CD3_fit <- predict(
  JM_CD3_allDLI,
  newdat,
  type = "Marginal",
  interval = "confidence",
  returnData = TRUE
) |>
  ggplot(aes(intSCT2_5, pred, col = ATG, group = ATG)) +
  geom_ribbon(aes(ymin = low, ymax = upp), col = NA, fill = "lightgray", alpha = 0.7) +
  geom_line(size = 1.25) +
  facet_wrap(~ CMV_PD) +
  scale_y_continuous(
    breaks = log(c(1, 5, 25, 100, 500, 2000)),
    labels = c(1, 5, 25, 100, 500, 2000)
  ) +
  labs(x = "Time since HSCT (months)", y = "CD3 cell counts") +
  theme_minimal()

ggsave(
  filename = "C:/Users/efbonneville/Downloads/immume-rec-talk/CD3_fit.png",
  plot = CD3_fit,
  dpi = 300,
  width = 12,
  height = 8,
  units = "in",
  scale = 0.6
)


exp(confint(JM_CD3_allDLI))


# CD4 model ---------------------------------------------------------------


predict(
  JM_CD4_allDLI,
  newdat,
  type = "Marginal",
  interval = "confidence",
  returnData = TRUE
) |>
  ggplot(aes(intSCT2_5, pred, col = ATG, group = ATG)) +
  geom_ribbon(aes(ymin = low, ymax = upp), col = NA, fill = "lightgray", alpha = 0.7) +
  geom_line(size = 1.25) +
  facet_wrap(~ CMV_PD) +
  scale_y_continuous(
    breaks = log(c(1, 5, 25, 100, 500, 2000)),
    labels = c(1, 5, 25, 100, 500, 2000)
  ) +
  labs(x = "Time since HSCT (months)", y = "CD4 cell counts") +
  theme_minimal()



# CD8 model ---------------------------------------------------------------



predict(
  JM_CD8_allDLI,
  newdat,
  type = "Marginal",
  interval = "confidence",
  returnData = TRUE
) |>
  ggplot(aes(intSCT2_5, pred, col = ATG, group = ATG)) +
  geom_ribbon(aes(ymin = low, ymax = upp), col = NA, fill = "lightgray", alpha = 0.7) +
  geom_line(size = 1.25) +
  facet_wrap(~ CMV_PD) +
  scale_y_continuous(
    breaks = log(c(1, 5, 25, 100, 500, 2000)),
    labels = c(1, 5, 25, 100, 500, 2000)
  ) +
  labs(x = "Time since HSCT (months)", y = "CD4 cell counts") +
  theme_minimal()
