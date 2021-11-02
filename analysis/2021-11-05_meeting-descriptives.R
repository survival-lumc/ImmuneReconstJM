# (Looking at descriptives due to poor model convergence)

theme_set(theme_bw(base_size = 14))

# Load datasets
tar_load(c("datasets"))
dat_wide <- datasets$wide
dat_long <- datasets$long

# Event numbers (cell_interv = modified t-cell product)
gtsummary::tbl_summary(dat_wide[, c("TCD2", "endpoint6_s")], by = "TCD2")

# Number of measurements per person
gtsummary::tbl_summary(dat_wide[, c("TCD2", "n_measurements")], by = "TCD2")
dat_wide |>
  ggplot(aes(n_measurements)) +
  geom_histogram(bins = 20, col = "black", fill = "lightblue") +
  facet_wrap(~ TCD2)



# Closer look at NMA cohort -----------------------------------------------


NMA_melted <- melt.data.table(
  data = dat_long[grepl(TCD, pattern = "^NMA")],
  id.vars = c("IDAA", "ATG", "CMV_PD", "intSCT2_5", "endpoint6_s"),
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
NMA_melted |>
  ggplot(aes(intSCT2_5, cell_count)) +
  geom_line(aes(group = IDAA), alpha = 0.85, col = "gray") +
  geom_smooth(
    method = "loess",
    se = FALSE,
    aes(group = ATG, col = ATG),
    size = 1.25,
    formula = y ~ x
  ) +
  facet_grid(endpoint6_s ~ cell_line) +
  scale_y_continuous(
    breaks = log(c(1, 5, 25, 100, 500, 2000)),
    labels = c(1, 5, 25, 100, 500, 2000)
  ) +
  labs(x = "Time since HSCT (months)", y = "Cell counts")


# For a few individuals plot all of their cell trajectories
set.seed(2022)
NMA_melted[IDAA %in% sample(IDAA, size = 24, replace = FALSE)] |>
  ggplot(aes(intSCT2_5, cell_count, col = cell_line)) +
  geom_point() +
  geom_line(aes(group = cell_line)) +
  facet_wrap(~ IDAA, ncol = 6) +
  scale_y_continuous(
    breaks = log(c(1, 5, 25, 100, 500, 2000)),
    labels = c(1, 5, 25, 100, 500, 2000)
  ) +
  labs(x = "Time since HSCT (months)", y = "Cell count")


# Closer look at cell lines
cell_lines <- paste0(c("CD3", "CD4", "CD8", "CD19", "NK"), "_abs_log")
names(cell_lines) <- cell_lines

grid_plots <- lapply(cell_lines, function(cell) {

  cell <- rlang::sym(cell)
  dat_long |>
    ggplot(aes(intSCT2_5, !!cell, col = IDAA)) +
    geom_line(aes(group = IDAA), alpha = 0.9)+#, col = ATG)) +
    geom_point(aes(group = IDAA))+#, col = ATG)) +
    facet_grid(endpoint6_s ~ CMV_PD * ATG) + #  ~ IDAA
    coord_cartesian(xlim = c(0, 24)) +
    theme_bw() +
    theme(legend.position = "none") +
    labs(x = "Time since HSCT (months)", y = cell)
})

grid_plots$CD4_abs_log



# All cohorts compared ----------------------------------------------------


melt.data.table(
  data = dat_long,
  id.vars = c("IDAA", "ATG", "CMV_PD", "intSCT2_5", "endpoint6_s", "TCD2"),
  measure.vars = paste0(c("CD4", "CD8", "CD19", "CD3", "NK"), "_abs_log"),
  variable.name = "cell_line",
  value.name = "cell_count"
) |>
  ggplot(aes(intSCT2_5 , cell_count)) +
  geom_line(aes(group = IDAA), alpha = 0.8, col = "gray") +
  geom_smooth(
    method = "loess",
    formula = y ~ x,
    size = 1.25,
    se = FALSE
  ) +
  facet_grid(cell_line ~ TCD2) +
  scale_y_continuous(
    breaks = log(c(1, 5, 25, 100, 500, 2000)),
    labels = c(1, 5, 25, 100, 500, 2000)
  ) +
  labs(x = "Time since HSCT (months)", y = "Cell count")


# Closer look at "clean" MA cohort ----------------------------------------



MA_RD_melted <- melt.data.table(
  data = dat_long[grepl(TCD, pattern = "^MA RD")],
  id.vars = c("IDAA", "ATG", "CMV_PD", "intSCT2_5", "endpoint6_s"),
  measure.vars = paste0(c("CD4", "CD8", "CD19", "CD3", "NK"), "_abs_log"),
  variable.name = "cell_line",
  value.name = "cell_count"
)


# Overall (by cell line)
MA_RD_melted |>
  ggplot(aes(intSCT2_5, cell_count)) +
  geom_line(aes(group = IDAA), alpha = 0.5, col = "lightgray") +
  geom_smooth(
    method = "loess",
    se = FALSE,
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
MA_RD_melted |>
  ggplot(aes(intSCT2_5, cell_count)) +
  geom_line(aes(group = IDAA), alpha = 0.85, col = "gray") +
  geom_smooth(
    method = "loess",
    se = FALSE,
    #aes(group = ATG, col = ATG),
    size = 1.25,
    formula = y ~ x
  ) +
  facet_grid(endpoint6_s ~ cell_line) +
  scale_y_continuous(
    breaks = log(c(1, 5, 25, 100, 500, 2000)),
    labels = c(1, 5, 25, 100, 500, 2000)
  ) +
  labs(x = "Time since HSCT (months)", y = "Cell counts")
