# Presentation figure 2022/02/02


# Source targets/models and functions -------------------------------------


# Focus on CD3 - since CD4 struggled with correlated random effects
tar_load(
  c(
    dat_merged,
    NMA_preDLI_datasets,
    NMA_postDLI_datasets_CD3_corr,
    preDLI_CD3_jointModel_corr,
    postDLI_CD3_jointModel_corr_both,
    postDLI_CD3_jointModel_corr_value

    # Remaining possible models to load:
    #preDLI_CD4_jointModel_both,
    #preDLI_CD8_jointModel_both,
    #NMA_postDLI_datasets_CD4,
    #NMA_postDLI_datasets_CD8,
    #postDLI_CD4_jointModel_both,
    #postDLI_CD8_jointModel_both
  )
)

# Any helper functions needed
source("data-raw/prepare-raw-data.R")
source("R/modelling-helpers.R")
source("R/plotting-helpers.R")

theme_set(theme_bw(base_size = 14))


# Subject
fitted_sub <- fitted(
  preDLI_CD3_jointModel_corr,
  process = "Longitudinal",
  type = "Marginal"
)
resid_sub <- residuals(
  preDLI_CD3_jointModel_corr,
  process = "Longitudinal",
  type = "Subject"
)
plot(fitted_sub, resid_sub)
abline(h = 0, lty = 3, col = "grey", lwd = 2)
panel.smooth(fitted_sub, resid_sub, lwd = 2)

# Marginal

plot(preDLI_CD3_jointModel_corr)

# pre-DLI: intro ----------------------------------------------------------

# Read and add label
dat_long <- NMA_preDLI_datasets$long
dat_wide <- NMA_preDLI_datasets$wide
table(dat_wide$endpoint7_s)

ggplot(dat_wide, aes(n_measurements)) +
  geom_histogram(bins = 13, col = "black", fill = "lightblue") +
  labs(x = "# Repeated measurements per patient") +
  theme_minimal(base_size = 14)

dat_long[, "endpoint_lab" := fcase(
  endpoint7_s == "cens", "Censored",
  endpoint7_s == "gvhd", "GVHD",
  endpoint7_s == "relapse", "Relapse",
  endpoint7_s == "other_nrf", "Other NRF"
)]
dat_long[, endpoint_lab := factor(
  endpoint_lab, levels = c("Censored", "GVHD", "Relapse", "Other NRF")
)]

# Added fitted vals and slope
dat_long[, ':=' (
  curr_val = fitted(preDLI_CD3_jointModel_corr, type = "Subject"),
  slope = fitted_slopes_long(preDLI_CD3_jointModel_corr)
)]

# For plotting: get tangent intercept, and then get coordinates
dat_long[, int_tang := get_int_tangent(x = intSCT2_5, y = curr_val, slope = slope)]
delta_tan <- 0.15 # Width of tangent arrow
dat_long[, ':=' (
  start_x = intSCT2_5 - delta_tan,
  start_y = int_tang + slope * (intSCT2_5 - delta_tan),
  end_x = intSCT2_5 + delta_tan,
  end_y = int_tang + slope * (intSCT2_5 + delta_tan)
)]

# Make a data with last measurement
dat_long_last <- dat_long[, .SD[.N], by = "IDAA"]

# Sample 16 patients
set.seed(649846561)
IDAA_subs <- sample(levels(dat_long$IDAA), replace = FALSE, size = 16)


# pre-DLI: raw + indiv fits -----------------------------------------------


# First raw-plot
p_raw_indiv <- ggplot(dat_long[IDAA %in% IDAA_subs], aes(intSCT2_5, CD3_abs_log)) +
  #ggplot(dat_long, aes(intSCT2_5, CD3_abs_log)) +
  geom_point(
    size = 3.5, pch = 21,
    alpha = 0.8,
    col = "#359fda",
    fill = colorspace::lighten("#359fda", amount = 0.3)
  ) +
  geom_vline(aes(xintercept = endpoint7), linetype = "dashed") +
  facet_wrap(~ IDAA) +
  labs(x = "Time since alloHCT (months)", y = "CD3 cell counts") +
  scale_y_continuous(
    breaks = log(c(5, 25, 100, 500, 1500)),
    labels = c(5, 25, 100, 500, 1500)
  ) +

  # Label endpoint
  geom_label(
    data = dat_long_last[IDAA %in% IDAA_subs],
    aes(x = endpoint7 + 0.25, y = log(2.5), label = endpoint_lab),
    #size = 5,
    hjust = 0,
    lineheight = .8,
    inherit.aes = FALSE,
    label.size = NA
  ) +
  geom_curve(
    data = dat_long_last[IDAA %in% IDAA_subs],
    #data = data.frame(x = 2, y = 29, xend = 2.5, yend = 20),
    mapping = aes(x = endpoint7 + 0.25, y = log(2.5), xend = endpoint7 + 0.025, yend = log(5)),
    colour = "black",
    size = 0.75,
    curvature = 0.2,
    arrow = arrow(length = unit(0.05, "npc"), type = "open"),
    inherit.aes = FALSE
  ) +
  coord_cartesian(xlim = c(0, 7.5))

p_raw_indiv

# Add individual fits
p_raw_indiv +
  geom_line(aes(y = curr_val), size = 1, col = "darkblue")


# Try with all cell_lines (maybe on original scale)
melt(
  data = dat_long[IDAA %in% IDAA_subs],
  id.vars = c("IDAA", "endpoint7", "intSCT2_5"),
  measure.vars = patterns("*_abs_log$"),
  value.name = "count",
  variable.name = "cell_line"
) |>
  ggplot(aes(intSCT2_5, exp(count)^(1/2))) +
  geom_point(
    aes(fill = cell_line, shape = cell_line, col = cell_line),
    size = 3.5,
    #pch = 21,
    alpha = 0.8#,
    #col = "#359fda",
    #fill = colorspace::lighten("#359fda", amount = 0.3)
  ) +
  geom_line(aes(linetype = cell_line, col = cell_line)) +
  geom_vline(aes(xintercept = endpoint7), linetype = "dashed") +
  facet_wrap(~ IDAA) +
  labs(x = "Time since alloHCT (months)", y = "CD3 cell counts") +
  #scale_y_continuous(
  #  breaks = log(c(5, 25, 100, 500, 1500)),
  #  labels = c(5, 25, 100, 500, 1500)
  #) +

  # Label endpoint
  geom_label(
    data = dat_long_last[IDAA %in% IDAA_subs],
    aes(x = endpoint7 + 0.25, y = log(2.5), label = endpoint_lab),
    #size = 5,
    hjust = 0,
    lineheight = .8,
    inherit.aes = FALSE,
    label.size = NA
  ) +
  geom_curve(
    data = dat_long_last[IDAA %in% IDAA_subs],
    #data = data.frame(x = 2, y = 29, xend = 2.5, yend = 20),
    mapping = aes(x = endpoint7 + 0.25, y = log(2.5), xend = endpoint7 + 0.025, yend = log(5)),
    colour = "black",
    size = 0.75,
    curvature = 0.2,
    arrow = arrow(length = unit(0.05, "npc"), type = "open"),
    inherit.aes = FALSE
  ) +
  coord_cartesian(xlim = c(0, 7.5))


# pre-DLI: marginal fit ---------------------------------------------------


newdat_preDLI <- expand.grid(
  "ATG" = levels(dat_long$ATG),
  "CMV_PD" = levels(dat_long$CMV_PD),
  "hirisk" = levels(dat_long$hirisk),
  "intSCT2_5" = seq(0, 6, by = 0.05)
)

dat_preds <- predict(
  preDLI_CD3_jointModel_corr,
  newdata = newdat_preDLI,
  type = "Marginal",
  idVar = "IDAA",
  returnData = TRUE,
  interval = "confidence"
)

# This plot already answers first q
ggplot(
  data = dat_preds,
  aes(
    x = intSCT2_5,
    y = pred,
    group = interaction(ATG, hirisk),
    col = ATG
  )
) +
  geom_ribbon(
    aes(ymin = low, ymax = upp),
    fill = "gray",
    alpha = 0.5,
    col = NA
  ) +
  geom_line(
    aes(linetype = hirisk),
    size = 1.5,
    alpha = 0.75
  ) +
  # add repel labels
  facet_wrap(facets = ~ CMV_PD) +
  labs(
    x = "Time since alloHCT (months)",
    y = "CD3 cell counts",
    linetype = "ITT",
    col = "Donor type"
  ) +
  scale_y_continuous(
    breaks = log(c(5, 25, 100, 500, 1500)),
    labels = c(5, 25, 100, 500, 1500)
  ) +
  xlim(c(0, 6)) +
  scale_color_manual(
    labels = c("Related", "Unrelated (+ATG)"),
    values = c("brown", "darkblue")
  ) +
  scale_linetype_manual(
    labels = c("No", "Yes (high risk)"),
    values = c("solid", "dotdash")
  )

# Pre-DLI: current value + slope indiv plots ------------------------------

ggplot(dat_long[IDAA %in% IDAA_subs], aes(intSCT2_5, CD3_abs_log)) +
  geom_point(
    size = 3.5, pch = 21,
    alpha = 0.3,
    col = "#359fda",
    fill = colorspace::lighten("#359fda", amount = 0.3)
  ) +
  geom_vline(aes(xintercept = endpoint7), linetype = "dashed") +
  facet_wrap(~ IDAA) +
  labs(x = "Time since alloHCT (months)", y = "CD3 cell counts") +
  scale_y_continuous(
    breaks = log(c(5, 25, 100, 500, 1500)),
    labels = c(5, 25, 100, 500, 1500)
  ) +

  # Label endpoint
  geom_label(
    data = dat_long_last[IDAA %in% IDAA_subs],
    aes(x = endpoint7 + 0.25, y = log(2.5), label = endpoint_lab),
    #size = 5,
    hjust = 0,
    lineheight = .8,
    inherit.aes = FALSE,
    label.size = NA
  ) +
  geom_curve(
    data = dat_long_last[IDAA %in% IDAA_subs],
    #data = data.frame(x = 2, y = 29, xend = 2.5, yend = 20),
    mapping = aes(x = endpoint7 + 0.25, y = log(2.5), xend = endpoint7 + 0.025, yend = log(5)),
    colour = "black",
    size = 0.75,
    curvature = 0.2,
    arrow = arrow(length = unit(0.05, "npc"), type = "open"),
    inherit.aes = FALSE
  ) +
  coord_cartesian(xlim = c(0, 7.5)) +
  geom_line(
    aes(y = curr_val), size = 0.75, col = "darkblue", linetype = "dashed",
    alpha = 0.5
  ) +
  geom_point(
    aes(y = curr_val), size = 1.5, col = "darkblue"
  )

# Individual fit + slopes
ggplot(dat_long[IDAA %in% IDAA_subs], aes(intSCT2_5, CD3_abs_log)) +
  geom_point(
    size = 3.5,
    pch = 21,
    alpha = 0.3,
    col = "#359fda",
    fill = colorspace::lighten("#359fda", amount = 0.3)
  ) +
  geom_vline(aes(xintercept = endpoint7), linetype = "dashed") +
  facet_wrap(~ IDAA) +
  labs(x = "Time since alloHCT (months)", y = "CD3 cell counts") +
  scale_y_continuous(
    breaks = log(c(5, 25, 100, 500, 1500)),
    labels = c(5, 25, 100, 500, 1500)
  ) +

  # Label endpoint
  geom_label(
    data = dat_long_last[IDAA %in% IDAA_subs],
    aes(x = endpoint7 + 0.25, y = log(2.5), label = endpoint_lab),
    #size = 5,
    hjust = 0,
    lineheight = .8,
    inherit.aes = FALSE,
    label.size = NA
  ) +
  geom_curve(
    data = dat_long_last[IDAA %in% IDAA_subs],
    #data = data.frame(x = 2, y = 29, xend = 2.5, yend = 20),
    mapping = aes(x = endpoint7 + 0.25, y = log(2.5), xend = endpoint7 + 0.025, yend = log(5)),
    colour = "black",
    size = 0.75,
    curvature = 0.2,
    arrow = arrow(length = unit(0.05, "npc"), type = "open"),
    inherit.aes = FALSE
  ) +
  coord_cartesian(xlim = c(0, 7.5)) +
  geom_line(
    aes(y = curr_val), size = 0.75, col = "darkblue", linetype = "dashed",
    alpha = 0.5
  ) +
  geom_segment(
    aes(x = start_x, y = start_y, xend = end_x, yend = end_y),
    size = 1,
    arrow = arrow(length = unit(0.04, "npc"), type = "open")
  )


# Table summary -----------------------------------------------------------

table_mod_summary(preDLI_CD3_jointModel_corr)


# Slope coef pt. 1 - direction --------------------------------------------

p_all_traj <- dat_long |>
  ggplot(aes(intSCT2_5, curr_val, group = IDAA)) +
  geom_line(aes(col = IDAA), alpha = 0.75, size = 1) +
  labs(x = "Time since alloHCT (months)", y = "CD3 cell counts") +
  scale_y_continuous(
    breaks = log(c(5, 25, 100, 500, 1500)),
    labels = c(5, 25, 100, 500, 1500)
  ) +
  facet_wrap(endpoint_lab ~ .) + # ATG
  theme(legend.position = "none")


p_all_traj

p_all_traj +
  geom_segment(
    data = dat_long_last,
    aes(x = start_x, y = start_y, xend = end_x, yend = end_y),
    size = 1,
    arrow = arrow(length = unit(0.02, "npc"), type = "closed")
  )

# Summarise current values and slopes at last moment
dat_long_last |>
  ggplot(aes(x = endpoint_lab, y = slope, col = endpoint_lab)) +
  #geom_violin() +
  geom_boxplot(width = .2, outlier.shape = NA) +
  geom_point(position = position_dodge2(width = 0.1), alpha = 0.5) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  labs(y = "Slope at measurement prior to endpoint \n (CD3 increase p/month)",
       x = "Endpoint") +
  scale_y_continuous(
    breaks = log(c(5, 25, 100, 500, 1500)),
    labels = c(5, 25, 100, 500, 1500),
    sec.axis = sec_axis(trans = ~ ., name = "(Increase in log CD3 p/month)")
  ) +
  theme(legend.position = "none")

# Same plot for current value
dat_long_last |>
  ggplot(aes(x = endpoint_lab, y = curr_val, col = endpoint_lab)) +
  #geom_violin() +
  geom_boxplot(width = .2, outlier.shape = NA) +
  geom_point(position = position_dodge2(width = 0.1), alpha = 0.5) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  labs(y = "'True' (current val.) CD3 at last \n measurement prior to endpoint") +
  scale_y_continuous(
    breaks = log(c(5, 25, 100, 500, 1500)),
    labels = c(5, 25, 100, 500, 1500),
    sec.axis = sec_axis(trans = ~ ., name = "(on log CD3 scale)")
  ) +
  theme(legend.position = "none")

# Check averages
avg_curr_vals <- dat_long_last[, .(avg_curr_val = mean(curr_val)), by = "endpoint7_s"]
vec_curr_vals <- setNames(avg_curr_vals[["avg_curr_val"]], avg_curr_vals[["endpoint7_s"]])
vec_curr_vals[-1] - vec_curr_vals[1]


# Slope coef pt. 2 - magnitude --------------------------------------------


# Illustrate increase of one in single patient
IDAA_single <- "5875"

p_single_slope <- ggplot(dat_long[IDAA %in% IDAA_single], aes(intSCT2_5, CD3_abs_log)) +
  geom_point(
    size = 3.5,
    pch = 21,
    alpha = 0.3,
    col = "#359fda",
    fill = colorspace::lighten("#359fda", amount = 0.3)
  ) +
  geom_vline(aes(xintercept = endpoint7), linetype = "dashed") +
  labs(x = "Time since alloHCT (months)", y = "CD3 cell counts") +
  scale_y_continuous(
    breaks = log(c(5, 25, 100, 500, 1500)),
    labels = c(5, 25, 100, 500, 1500)
  ) +

  # Label endpoint
  geom_label(
    data = dat_long_last[IDAA %in% IDAA_single],
    aes(x = endpoint7 + 0.25, y = log(2.5), label = endpoint_lab),
    #size = 5,
    hjust = 0,
    lineheight = .8,
    inherit.aes = FALSE,
    label.size = NA
  ) +
  geom_curve(
    data = dat_long_last[IDAA %in% IDAA_single],
    #data = data.frame(x = 2, y = 29, xend = 2.5, yend = 20),
    mapping = aes(x = endpoint7 + 0.25, y = log(2.5), xend = endpoint7 + 0.025, yend = log(5)),
    colour = "black",
    size = 0.75,
    curvature = 0.2,
    arrow = arrow(length = unit(0.05, "npc"), type = "open"),
    inherit.aes = FALSE
  ) +
  coord_cartesian(xlim = c(0, 6), ylim = c(log(1), log(2000))) +
  geom_line(
    aes(y = curr_val), size = 0.75, col = "darkblue", linetype = "dashed",
    alpha = 0.5
  )

p_single_slope +
  geom_segment(
    data = dat_long_last[IDAA %in% IDAA_single],
    aes(x = start_x - 0.3, y = start_y, xend = end_x + 0.3, yend = end_y),
    size = 1,
    arrow = arrow(length = unit(0.02, "npc"), type = "closed")
  )

# Calculate what increase in slope looks like
slope_plus1 <- copy(dat_long_last[IDAA %in% IDAA_single])
slope_plus1[, int_tang := get_int_tangent(x = intSCT2_5, y = curr_val, slope = slope + 1)]
slope_plus1[, ':=' (
  start_x = intSCT2_5 - delta_tan,
  start_y = int_tang + (slope + 1) * (intSCT2_5 - delta_tan),
  end_x = intSCT2_5 + delta_tan,
  end_y = int_tang + (slope + 1) * (intSCT2_5 + delta_tan)
)]

p_single_slope +
  geom_segment(
    data = dat_long_last[IDAA %in% IDAA_single],
    alpha = 0.33,
    aes(x = start_x - 0.3, y = start_y, xend = end_x + 0.3, yend = end_y),
    size = 1,
    arrow = arrow(length = unit(0.02, "npc"), type = "closed")
  ) +
  geom_segment(
    data = slope_plus1,
    aes(x = start_x - 0.3, y = start_y, xend = end_x + 0.3, yend = end_y),
    size = 1,
    arrow = arrow(length = unit(0.02, "npc"), type = "closed")
  ) +
  annotate("text", label = "Increase in slope by 1 unit \n (increase log CD3 p/month)",
            x = 5, y = log(1000), hjust = 0, size = 5)


# Range of slope variable
ggplot(dat_long, aes(slope)) +
  geom_histogram(bins = 30, fill = "lightblue", col = "black") +
  facet_wrap(~ endpoint_lab)

ggplot(dat_long, aes(slope)) +
  geom_histogram(bins = 30, fill = "lightblue", col = "black") +
  scale_x_continuous(
    breaks = log(c(5, 25, 100, 500, 1500)),
    labels = c(5, 25, 100, 500, 1500)
  ) +
  facet_wrap(~ endpoint_lab)


# Slope over time
dat_long[IDAA %in% IDAA_single] |>
  ggplot(aes(intSCT2_5, slope)) +
  geom_line()


# Pretend baseline hazard is contant value 1, holding everything else fixed
# Maybe just show rizopolous plot p .103
dat_long[IDAA %in% IDAA_single] |>
  ggplot(aes(intSCT2_5, slope)) +
  geom_hline(yintercept = 1) +
  geom_line(
    aes(y = slope * preDLI_CD3_jointModel_corr$coefficients$Dalpha[2]) # relapse
  ) +
  labs(y = "Hazard")


# Post-DLI ----------------------------------------------------------------


dat_long_post <- NMA_postDLI_datasets_CD3_corr$long
dat_wide_post <- NMA_postDLI_datasets_CD3_corr$wide
table(dat_wide_post$sec_endpoint_s)
ggplot(dat_wide_post, aes(n_measurements)) +
  geom_histogram(bins = 13, col = "black", fill = "lightblue") +
  labs(x = "# Repeated measurements per patient") +
  theme_minimal(base_size = 14)

dat_long_post[, "endpoint_lab" := fcase(
  sec_endpoint_s == "cens", "Event-free",
  sec_endpoint_s == "gvhd", "GVHD",
  sec_endpoint_s == "rel_nrf", "Relapse or Other NRF"
)]

# Added fitted vals and slope
dat_long_post[, ':=' (
  curr_val = fitted(postDLI_CD3_jointModel_corr_both, type = "Subject"),
  slope = fitted_slopes_long(postDLI_CD3_jointModel_corr_both)
)]

# For plotting: get tangent intercept, and then get coordinates
dat_long_post[, int_tang := get_int_tangent(x = intSCT2_5_reset, y = curr_val, slope = slope)]
delta_tan <- 0.3 # Width of tangent arrow
dat_long_post[, ':=' (
  start_x = intSCT2_5_reset - delta_tan,
  start_y = int_tang + slope * (intSCT2_5_reset - delta_tan),
  end_x = intSCT2_5_reset + delta_tan,
  end_y = int_tang + slope * (intSCT2_5_reset + delta_tan)
)]

# Make a data with last measurement
dat_long_post_last <- dat_long_post[, .SD[.N], by = "IDAA"]


# Raw general
dat_long_post |>
  ggplot(aes(intSCT2_5, CD3_abs_log)) +
  geom_line(aes(col = IDAA), alpha = 0.75, size = 1) +
  scale_y_continuous(
    breaks = log(c(5, 25, 100, 500, 1500)),
    labels = c(5, 25, 100, 500, 1500)
  ) +
  coord_cartesian(ylim = c(0, 8), xlim = c(0, 24)) +
  labs(x = "Time since alloHCT (months)", y = "CD3 cell counts") +
  theme(legend.position = "none")

# Reset
dat_long_post |>
  ggplot(aes(intSCT2_5_reset, CD3_abs_log)) +
  geom_line(aes(col = IDAA), alpha = 0.75, size = 1) +
  scale_y_continuous(
    breaks = log(c(5, 25, 100, 500, 1500)),
    labels = c(5, 25, 100, 500, 1500)
  ) +
  coord_cartesian(ylim = c(0, 8), xlim = c(0, 12)) +
  geom_vline(xintercept = 12, linetype = "dashed") +
  labs(x = "Time since first DLI (months)", y = "CD3 cell counts") +
  theme(legend.position = "none") +
  annotate("text", label = "Admin. censoring at 12 months \n after first DLI",
           x = 9, y = log(10), hjust = 0.5, size = 5) +
  geom_segment(
    data = data.frame(start_x = 10.5, start_y = log(10), end_x = 12, end_y = log(25)),
    aes(x = start_x, y = start_y, xend = end_x, yend = end_y),
    size = 1,
    arrow = arrow(length = unit(0.02, "npc"), type = "closed")
  )


# Divide the trajectories by endpoint
dat_long_post |>
  ggplot(aes(intSCT2_5_reset, CD3_abs_log)) +
  geom_line(aes(col = IDAA), alpha = 0.75, size = 1) +
  scale_y_continuous(
    breaks = log(c(5, 25, 100, 500, 1500)),
    labels = c(5, 25, 100, 500, 1500)
  ) +
  coord_cartesian(ylim = c(0, 8), xlim = c(0, 12)) +
  labs(x = "Time since first DLI (months)", y = "CD3 cell counts") +
  theme(legend.position = "none") +
  facet_wrap(~ endpoint_lab)



# Raw post DLI: GVHD and Relapse/NRF --------------------------------------

p_post_gvhd <- ggplot(
  dat_long_post[sec_endpoint_s == "gvhd"],
  aes(intSCT2_5_reset, CD3_abs_log)
) +
  geom_point(
    size = 3.5,
    pch = 21,
    alpha = 0.8,
    col = "#359fda",
    fill = colorspace::lighten("#359fda", amount = 0.3)
  ) +
  geom_point(
    data = dat_wide_post[sec_endpoint_s == "gvhd"],
    aes(0, preds_subj),
    size = 3.5,
    pch = 21,
    alpha = 0.8,
    col = "black",
    fill = "darkblue"#colorspace::lighten("red", amount = 0.3)
  ) +
  labs(x = "Time since DLI (months)", y = "CD3 cell counts") +
  scale_y_continuous(
    breaks = log(c(5, 25, 100, 500, 1500)),
    labels = c(5, 25, 100, 500, 1500)
  ) +
  facet_wrap(~ IDAA) +
  coord_cartesian(ylim = c(0, 8), xlim = c(0, 12))

p_post_relnrf <- ggplot(
  dat_long_post[sec_endpoint_s == "rel_nrf"],
  aes(intSCT2_5_reset, CD3_abs_log)
) +
  geom_point(
    size = 3.5,
    pch = 21,
    alpha = 0.8,
    col = "#359fda",
    fill = colorspace::lighten("#359fda", amount = 0.3)
  ) +
  geom_point(
    data = dat_wide_post[sec_endpoint_s == "rel_nrf"],
    aes(0, preds_subj),
    size = 3.5,
    pch = 21,
    alpha = 0.8,
    col = "black",
    fill = "darkblue"#colorspace::lighten("red", amount = 0.3)
  ) +
  labs(x = "Time since DLI (months)", y = "CD3 cell counts") +
  scale_y_continuous(
    breaks = log(c(5, 25, 100, 500, 1500)),
    labels = c(5, 25, 100, 500, 1500)
  ) +
  facet_wrap(~ IDAA) +
  coord_cartesian(ylim = c(0, 8), xlim = c(0, 12))

p_post_gvhd
p_post_relnrf


# Include censored too
ggplot(
  dat_long_post[sec_endpoint_s == "cens"],
  aes(intSCT2_5_reset, CD3_abs_log)
) +
  geom_point(
    size = 3.5,
    pch = 21,
    alpha = 0.8,
    col = "#359fda",
    fill = colorspace::lighten("#359fda", amount = 0.3)
  ) +
  geom_point(
    data = dat_wide_post[sec_endpoint_s == "cens"],
    aes(0, preds_subj),
    size = 3.5,
    pch = 21,
    alpha = 0.8,
    col = "black",
    fill = "darkblue"#colorspace::lighten("red", amount = 0.3)
  ) +
  labs(x = "Time since DLI (months)", y = "CD3 cell counts") +
  scale_y_continuous(
    breaks = log(c(5, 25, 100, 500, 1500)),
    labels = c(5, 25, 100, 500, 1500)
  ) +
  facet_wrap(~ IDAA) +
  coord_cartesian(ylim = c(0, 8), xlim = c(0, 12)) +
  geom_line(aes(y = curr_val), size = 1, col = "darkblue")


# Fits --------------------------------------------------------------------


p_post_gvhd +
  geom_line(aes(y = curr_val), size = 1, col = "darkblue")

p_post_relnrf +
  geom_line(aes(y = curr_val), size = 1, col = "darkblue")

# What about the slopes
p_post_gvhd +
  geom_line(
    aes(y = curr_val), size = 0.75, col = "darkblue", linetype = "dashed",
    alpha = 0.5
  ) +
  geom_segment(
    aes(x = start_x, y = start_y, xend = end_x, yend = end_y),
    size = 1,
    arrow = arrow(length = unit(0.04, "npc"), type = "open")
  )

p_post_relnrf +
  geom_line(
    aes(y = curr_val), size = 0.75, col = "darkblue", linetype = "dashed",
    alpha = 0.5
  ) +
  geom_segment(
    aes(x = start_x, y = start_y, xend = end_x, yend = end_y),
    size = 1,
    arrow = arrow(length = unit(0.04, "npc"), type = "open")
  )



# Do the marginal fit -----------------------------------------------------


newdat_post <- expand.grid(
  "ATG" = levels(dat_long_post$ATG),
  "earlylow_DLI" = levels(dat_long_post$earlylow_DLI),
  "intSCT2_5_reset" = seq(0, 12, by = 0.05)
)

dat_preds_post <- predict(
  postDLI_CD3_jointModel_corr_both,
  newdata = newdat_post,
  type = "Marginal",
  idVar = "IDAA",
  returnData = TRUE,
  interval = "confidence"
)

# This plot already answers first q
ggplot(
  data = dat_preds_post,
  aes(
    x = intSCT2_5_reset,
    y = pred,
    group = interaction(ATG, earlylow_DLI),
    col = earlylow_DLI
  )
) +
  geom_ribbon(
    aes(ymin = low, ymax = upp),
    fill = "gray",
    alpha = 0.5,
    col = NA
  ) +
  geom_line(
    aes(linetype = ATG),
    size = 1.5,
    alpha = 0.75
  ) +
  scale_y_continuous(
    breaks = log(c(5, 25, 100, 500, 1500)),
    labels = c(5, 25, 100, 500, 1500)
  ) +
  coord_cartesian(ylim = c(0, 8), xlim = c(0, 12)) +
  labs(
    x = "Time since alloHCT (months)",
    y = "CD3 cell counts",
    col = "Actual DLI received",
    linetype = "Donor type"
  ) +
  scale_linetype_manual(
    labels = c("Related", "Unrelated (+ATG)"),
    values = c("solid", "dotdash")
  ) +
  scale_color_manual(
    labels = c("Early/low dose", "Late/high dose"),
    values = c("brown", "darkblue")
  )


# .. think about goodness of fit things?
plot(postDLI_CD3_jointModel_corr_both)



# Attempt at a table post DLI ---------------------------------------------


# Event summary extract
summ_mod <- summary(postDLI_CD3_jointModel_corr_both)
surv_tab <- summ_mod$`CoefTable-Event`

# Make table
summary_raw <- data.table::data.table(surv_tab, keep.rownames = TRUE)
data.table::setnames(summary_raw, "rn", "Coefficient")
dt <- summary_raw[grep("^bs", x = Coefficient, invert = TRUE)]

dt[, event_num := regmatches(
  x = Coefficient,
  m = regexpr(pattern = "[+1-9]$", text = Coefficient)
)]

dt[, event := factor(
  x = event_num,
  levels = seq_len(2),
  labels = c("GVHD", "Relapse or Other NRF")
)]

dt[, Coefficient := gsub(x = Coefficient, pattern = "\\:.*", replacement = "")]
dt[, Coefficient := gsub(x = Coefficient, pattern = "\\.[+1-9]$", replacement = "")]
dt[, Coefficient := factor(
  x = Coefficient,
  levels = c("preds_subj", "earlylow_DLI", "hirisk", "Assoct", "Assoct.s"),
  labels = c("log(CD3)* at DLI", "Early/low DLI",  "High risk (ITT)", "Current value", "Slope")
)]
data.table::setorder(dt, event, Coefficient)
dt_edit <- dt[, c("event", "Coefficient", "Value", "Std.Err", "p-value")]
data.table::setnames(dt_edit, new = c("Event", "Coefficient", "log(HR)", "SE", "p-value"))
dt_edit[, "HR [95% CI]" := paste0(
  as.character(round(exp(`log(HR)`), 3)), " [",
  as.character(round(exp(`log(HR)` - qnorm(0.975) * SE), 3)), ";",
  as.character(round(exp(`log(HR)` + qnorm(0.975) * SE), 3)), "]"
)]
dt_edit[, Event := NULL]
# Change order?
dt_edit <- dt_edit[, c("Coefficient", "HR [95% CI]", "log(HR)", "SE", "p-value")]
dt_edit[, c("log(HR)", "SE", "p-value") := lapply(.SD, function(x) round(x, 3)),
        .SDcols = c("log(HR)", "SE", "p-value")]
dt_edit

dt_edit %>%
  kableExtra::kbl(format = "html") %>%
  kable_paper("striped", full_width = F) %>%
  pack_rows("GVHD", 1, 4) %>%
  pack_rows("Relapse or Other NRF", 5, 8)


# Show vcov mat?
summary(postDLI_CD3_jointModel_corr_both)


# Same plots as before with slopes ----------------------------------------

# The fitted ones
p_all_traj_post <- dat_long_post |>
  ggplot(aes(intSCT2_5_reset, curr_val, group = IDAA)) +
  geom_line(aes(col = IDAA), alpha = 0.75, size = 1) +
  labs(x = "Time since alloHCT (months)", y = "CD3 cell counts") +
  scale_y_continuous(
    breaks = log(c(5, 25, 100, 500, 1500)),
    labels = c(5, 25, 100, 500, 1500)
  ) +
  facet_wrap(endpoint_lab ~ .) + # ATG
  theme(legend.position = "none") +
  coord_cartesian(ylim = c(0, 8), xlim = c(0, 12))

p_all_traj_post

p_all_traj_post +
  geom_segment(
    data = dat_long_post_last,
    aes(x = start_x, y = start_y, xend = end_x, yend = end_y),
    size = 1,
    arrow = arrow(length = unit(0.02, "npc"), type = "closed")
  )

# Mixed mod tries uncorrelated reffs --------------------------------------



dat_long[sec_endpoint_s == "rel_nrf"] |>
  ggplot(aes(intSCT2_5_reset, CD3_abs_log, group = IDAA)) +
  geom_point() +
  facet_wrap(~ IDAA)

mod <- run_preDLI_longitudinal(
  cell_line = "CD3_abs_log",
  form_fixed = "ns(intSCT2_5_reset, 3) * earlylow_DLI + ATG", #"ns(intSCT2_5_reset, 3)", #"ns(intSCT2_5_reset, 3) * earlylow_DLI + ATG",
  form_random = ~ ns(intSCT2_5_reset, 3) | IDAA,
  #form_random = list(IDAA = pdDiag(~ ns(intSCT2_5_reset, 3))),
  dat = dat_long
)

summary(mod)

dat_long[, "preds" := predict(mod)]

dat_long[sec_endpoint_s == "cens"] |>
  ggplot(aes(intSCT2_5_reset, CD3_abs_log, group = IDAA)) +
  geom_point() +
  #geom_point(aes(y = preds), col = "blue") +
  geom_line(aes(y = preds), col = "blue") +
  facet_wrap(~ IDAA) +
  coord_cartesian(ylim = c(0, 8))

dat_long[sec_endpoint_s == "gvhd"] |>
  ggplot(aes(intSCT2_5_reset, CD3_abs_log, group = IDAA)) +
  geom_point() +
  #geom_point(aes(y = preds), col = "blue") +
  geom_line(aes(y = preds), col = "blue") +
  facet_wrap(~ IDAA) +
  coord_cartesian(ylim = c(0, 8))

dat_long[sec_endpoint_s == "rel_nrf"] |>
  ggplot(aes(intSCT2_5_reset, CD3_abs_log, group = IDAA)) +
  geom_point() +
  #geom_point(aes(y = preds), col = "blue") +
  geom_line(aes(y = preds), col = "blue") +
  facet_wrap(~ IDAA) +
  coord_cartesian(ylim = c(0, 8))

# Marginals
newdat_jm <- expand.grid(
  "ATG" = levels(dat_long$ATG),
  #"CMV_PD" = levels(dat_long$CMV_PD),
  "earlylow_DLI" = levels(dat_long$earlylow_DLI),
  "intSCT2_5_reset" = seq(0, 12, by = 0.05)
)

newdat_jm$preds_marg <- predict(mod, newdata = newdat_jm, level = 0L)

newdat_jm |>
  ggplot(aes(intSCT2_5_reset, preds_marg)) +
  geom_line(aes(group = interaction(ATG, earlylow_DLI),
                col = interaction(ATG, earlylow_DLI),
                linetype = interaction(ATG, earlylow_DLI)))  +
  coord_cartesian(ylim = c(0, 8))



# Script-specific functions -----------------------------------------------


# Make this eventually into one big table/forest plot
table_mod_summary <- function(model) {

  # Event summary extract
  summ_mod <- summary(model)
  surv_tab <- summ_mod$`CoefTable-Event`

  # Make table
  summary_raw <- data.table::data.table(surv_tab, keep.rownames = TRUE)
  data.table::setnames(summary_raw, "rn", "Coefficient")
  dt <- summary_raw[grep("^bs", x = Coefficient, invert = TRUE)]

  dt[, event_num := regmatches(
    x = Coefficient,
    m = regexpr(pattern = "[+1-9]$", text = Coefficient)
  )]

  dt[, event := factor(
    x = event_num,
    levels = seq_len(3),
    labels = c("GVHD", "Relapse", "NRF: Other")
  )]

  dt[, Coefficient := gsub(x = Coefficient, pattern = "\\:.*", replacement = "")]
  dt[, Coefficient := gsub(x = Coefficient, pattern = "\\.[+1-9]$", replacement = "")]
  dt[, Coefficient := factor(
    x = Coefficient,
    levels = c("ATG", "hirisk", "Assoct", "Assoct.s"),
    labels = c("Unrel. donor (+ATG)", "High risk (ITT)", "Current value", "Slope")
  )]
  data.table::setorder(dt, event, Coefficient)
  dt_edit <- dt[, c("event", "Coefficient", "Value", "Std.Err", "p-value")]
  data.table::setnames(dt_edit, new = c("Event", "Coefficient", "log(HR)", "SE", "p-value"))
  dt_edit[, "HR [95% CI]" := paste0(
    as.character(round(exp(`log(HR)`), 3)), " [",
    as.character(round(exp(`log(HR)` - qnorm(0.975) * SE), 3)), ";",
    as.character(round(exp(`log(HR)` + qnorm(0.975) * SE), 3)), "]"
  )]
  dt_edit[, Event := NULL]
  # Change order?
  dt_edit <- dt_edit[, c("Coefficient", "HR [95% CI]", "log(HR)", "SE", "p-value")]
  dt_edit[, c("log(HR)", "SE", "p-value") := lapply(.SD, function(x) round(x, 3)),
          .SDcols = c("log(HR)", "SE", "p-value")]
  dt_edit

  dt_edit %>%
    kableExtra::kbl(format = "html") %>%
    kable_paper("striped", full_width = F) %>%
    pack_rows("GVHD", 1, 4) %>%
    pack_rows("Relapse", 5, 8) %>%
    pack_rows("NRF: Other", 9, 11)
}


