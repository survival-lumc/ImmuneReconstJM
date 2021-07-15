# Load all objects
tar_load(
  c(
    datasets,
    dli_msdata,
    long_submodels,
    cox_submodel_all_dli,
    cox_submodel_no_dli,
    CD4_all_dli,
    CD4_all_dli_avfform,
    CD8_all_dli,
    CD8_all_dli_avfform,
    reference_values
  )
)

# Rename wide data
dat_wide <- datasets$wide

# For longitudinal predictions
newdat <- expand.grid(
  "ATG" = levels(dat_wide$ATG),
  "VCMVPAT_pre" = levels(dat_wide$VCMVPAT_pre),
  "intSCT2_5" = seq(0.1, 24, by = 0.1)
)

# Cox model ---------------------------------------------------------------

summary(cox_submodel_all_dli)


# CD4 JM: comparing cox ---------------------------------------------------


# Without interaction with DLI in association
summary(CD4_all_dli_avfform)
CD4_all_dli_avfform$coefficients$gammas
cox_submodel_all_dli$coefficients
CD4_all_dli_avfform$coefficients$alpha

# With interaction - DLI.2 explodes? Check with sim data
summary(CD4_all_dli)
CD4_all_dli$coefficients$gammas
cox_submodel_all_dli$coefficients
CD4_all_dli$coefficients$alpha

confint(CD4_all_dli)


long_submodels$CD4_abs_log$coefficients$fixed
CD4_all_dli$coefficients$betas



# Concerns on event number pre-post DLI? ----------------------------------

table(dat_wide$uDLI_s, dat_wide$endpoint6_s)


# CD4 longitudinal --------------------------------------------------------


JM::predict.jointModel(
  CD4_all_dli,
  newdata = newdat,
  type = "Marginal",
  idVar = "IDAA",
  returnData = TRUE,
  interval = "confidence"
) %>%
  ggplot(aes(x = intSCT2_5, y = pred)) +
  geom_ribbon(
    aes(ymin = low, ymax = upp, group = ATG, col = ATG),
    fill = "gray",
    alpha = 0.5,
    col = NA
  ) +
  geom_line(aes(group = ATG, col = ATG), size = 1.5, alpha = 0.75) +
  #geom_text(aes(label = label), hjust = 0, na.rm = TRUE, fontface = "bold") +
  facet_wrap(facets = ~ VCMVPAT_pre) +
  labs(x = "Time since alloHCT (months)", y = "CD4 cell counts") +
  # scale_y_continuous(
  #   breaks = log(c(5, 25, 100, 500, 1500)),
  #   labels = c(5, 25, 100, 500, 1500)
  # ) +
  theme_minimal() +
  #theme(legend.position = "none") +
  #coord_cartesian(xlim = c(0, 14)) +
  scale_color_brewer(palette = "Paired") +
  geom_vline(xintercept = 6, linetype = "dotted")


JM::predict.jointModel(
  CD4_all_dli,
  newdata = newdat,
  type = "Marginal",
  idVar = "IDAA",
  returnData = TRUE,
  interval = "confidence"
) %>%
  ggplot(aes(x = intSCT2_5, y = pred)) +
  geom_hline(
    yintercept = log(c())
  ) +
  geom_ribbon(
    aes(ymin = low, ymax = upp, group = ATG, col = ATG),
    fill = "gray",
    alpha = 0.5,
    col = NA
  ) +
  geom_line(aes(group = ATG, col = ATG), size = 1.5, alpha = 0.75) +
  #geom_text(aes(label = label), hjust = 0, na.rm = TRUE, fontface = "bold") +
  facet_wrap(facets = ~ VCMVPAT_pre) +
  labs(x = "Time since alloHCT (months)", y = "CD4 cell counts") +
   scale_y_continuous(
     breaks = log(c(5, 25, 100, 500, 1500)),
     labels = c(5, 25, 100, 500, 1500)
   ) +
  theme_minimal() +
  #theme(legend.position = "none") +
  #coord_cartesian(xlim = c(0, 14)) +
  scale_color_brewer(palette = "Paired")


# CD8 longitudinal --------------------------------------------------------


JM::predict.jointModel(
  CD8_all_dli,
  newdata = newdat,
  type = "Marginal",
  idVar = "IDAA",
  returnData = TRUE,
  interval = "confidence"
) %>%
  ggplot(aes(x = intSCT2_5, y = pred)) +
  geom_hline(yintercept = log(c(260, 990)), linetype = "dashed") +
  geom_ribbon(
    aes(ymin = low, ymax = upp, group = ATG, col = ATG),
    fill = "gray",
    alpha = 0.5,
    col = NA
  ) +
  geom_line(aes(group = ATG, col = ATG), size = 1.5, alpha = 0.75) +
  #geom_text(aes(label = label), hjust = 0, na.rm = TRUE, fontface = "bold") +
  facet_wrap(facets = ~ VCMVPAT_pre) +
  labs(x = "Time since alloHCT (months)", y = "CD8 cell counts") +
  scale_y_continuous(
    breaks = log(c(5, 25, 100, 500, 1500)),
    labels = c(5, 25, 100, 500, 1500)
  ) +
  theme_minimal() +
  #theme(legend.position = "none") +
  #coord_cartesian(xlim = c(0, 14)) +
  scale_color_brewer(palette = "Paired")


# CD8 cox -----------------------------------------------------------------




# Without interaction with DLI in association
summary(CD8_all_dli_avfform)
summary(CD8_all_dli)

CD4_all_dli_avfform$coefficients$gammas
cox_submodel_all_dli$coefficients
CD4_all_dli_avfform$coefficients$alpha

# With interaction - DLI.2 explodes? Check with sim data
summary(CD4_all_dli)
CD4_all_dli$coefficients$gammas
cox_submodel_all_dli$coefficients
CD4_all_dli$coefficients$alpha



# Da plot -----------------------------------------------------------------



dat_long_melt <- melt.data.table(
  data = datasets$long,
  id.vars = c("IDAA", "intSCT2_5", "endpoint6_s", "VCMVPAT_pre", "ATG"),
  variable.name = "cell_type",
  value.name = "counts",
  measure.vars = paste0(c("CD4", "CD8"), "_abs_log")
)
dat_long_melt[, cell_type := factor(cell_type, labels = c("CD4", "CD8"))]

preds_df_CD4 <- JM::predict.jointModel(
  CD4_all_dli,
  newdata = newdat,
  type = "Marginal",
  idVar = "IDAA",
  returnData = TRUE,
  interval = "confidence"
)

preds_df_CD8 <- JM::predict.jointModel(
  CD8_all_dli,
  newdata = newdat,
  type = "Marginal",
  idVar = "IDAA",
  returnData = TRUE,
  interval = "confidence"
)

preds_df <- rbindlist(list(preds_df_CD4, preds_df_CD8), idcol = "cell_type")
preds_df[, cell_type := factor(cell_type, levels = c(1, 2), labels = c("CD4", "CD8"))]


# Plot 1: predictions and overlay of endpoints
preds_df %>%
  ggplot(aes(x = intSCT2_5, y = pred)) +
  geom_hline(
    data = reference_values[cell_type %in% c("CD4", "CD8")],
    aes(yintercept = log(lower_limit)),
    linetype = "dotted"
  ) +
  geom_hline(
    data = reference_values[cell_type %in% c("CD4", "CD8")],
    aes(yintercept = log(upper_limit)),
    linetype = "dotted"
  ) +
  geom_ribbon(
    aes(ymin = low, ymax = upp, group = ATG, col = ATG),
    fill = "gray",
    alpha = 0.5,
    col = NA
  ) +
  geom_line(aes(group = ATG, col = ATG), size = 1.5, alpha = 0.75) +
  labs(x = "Time since alloHCT (months)", y = "Cell counts") +
  scale_y_continuous(
    breaks = log(c(5, 25, 100, 500, 1500)),
    labels = c(5, 25, 100, 500, 1500)
  ) +
  theme_minimal() +
  scale_color_brewer(palette = "Paired") +
  ggnewscale::new_scale_colour() +
  geom_point(
    data = dat_long_melt[, .SD[.N], by = c("IDAA", "cell_type")],
    aes(shape = endpoint6_s, colour = endpoint6_s, y = counts)
    #aes(shape = ATG, colour = endpoint6_s, y = counts)
  ) +
  facet_grid(facets = cell_type ~ VCMVPAT_pre)


# Raw plots ---------------------------------------------------------------


set.seed(1986)
IDAA_samp <- dat_wide[, .(
  IDAA_samp = sample(IDAA, size = floor(.N / 3), replace = FALSE)
), by = "endpoint6_s"][["IDAA_samp"]]

# Try dividing three groups
IDAA_grps <- dat_wide[, .(
  IDAA = unique(IDAA),
  IDAA_grps = sample(seq_len(3), size = length(unique(IDAA)), replace = TRUE, rep(1/3, 3))
), by = "endpoint6_s"]


tar_load(dat_merged)
raw_long <- dat_merged %>%
  melt.data.table(
    id.vars = c("IDAA", "intSCT2_5", "endpoint6_s", "VCMVPAT_pre"),
    variable.name = "cell_type",
    value.name = "counts",
    measure.vars = paste0(c("CD4", "CD8"), "_abs_log")
  )
raw_long <- raw_long[intSCT2_5 <= 24]
raw_long[, ':=' (
  cell_type = factor(cell_type, labels = c("CD4", "CD8")),
  endpoint6_s = factor(
    endpoint6_s,
    levels = c(
      "censored",
      "7 days after cellular intervention",
      "relapse",
      "non-relapse failure: other",
      "non-relapse failure: GvHD"
    ),
    labels = c("cens", "cell_interv", "REL", "NRF_other", "NRF_gvhd")
  )
)]

grp <- 1
raw_long[, highlight_pats := ifelse(IDAA %in% IDAA_grps[IDAA_grps == grp]$IDAA, 1, 0.65)]

p <- raw_long[endpoint6_s != "cell_interv"] %>%
  ggplot(aes(intSCT2_5, counts)) +
  geom_hline(
    data = reference_values[cell_type %in% c("CD4", "CD8")],
    aes(yintercept = log(lower_limit)),
    linetype = "dotted"
  ) +
  geom_hline(
    data = reference_values[cell_type %in% c("CD4", "CD8")],
    aes(yintercept = log(upper_limit)),
    linetype = "dotted"
  ) +
  geom_line(aes(group = IDAA, col = IDAA, alpha = highlight_pats)) +
  geom_point(
    data = dat_long_melt[endpoint6_s != "cell_interv", .SD[.N], by = c("IDAA", "cell_type")][
      , highlight_pats := ifelse(IDAA %in% IDAA_grps[IDAA_grps == grp]$IDAA, 1, 0.6)
    ],
    aes(y = counts, alpha = highlight_pats),
    shape = 4
    #aes(shape = ATG, colour = endpoint6_s, y = counts)
  ) +
  scale_y_continuous(
    breaks = log(c(5, 25, 100, 500, 1500)),
    labels = c(5, 25, 100, 500, 1500)
  ) +
  facet_grid(facets = cell_type ~ endpoint6_s) +
  theme_minimal() +
  theme(legend.position = "none")

plotly::ggplotly(p)

# - try multivar model again

p <- raw_long[endpoint6_s != "cell_interv"] %>%
  ggplot(aes(intSCT2_5, counts)) +
  geom_hline(
    data = reference_values[cell_type %in% c("CD4", "CD8")],
    aes(yintercept = log(lower_limit)),
    linetype = "dotted"
  ) +
  geom_hline(
    data = reference_values[cell_type %in% c("CD4", "CD8")],
    aes(yintercept = log(upper_limit)),
    linetype = "dotted"
  ) +
  geom_line(aes(group = IDAA, col = IDAA)) +
  geom_point(
    data = dat_long_melt[endpoint6_s != "cell_interv", .SD[.N], by = c("IDAA", "cell_type")][
      , highlight_pats := ifelse(IDAA %in% IDAA_grps[IDAA_grps == grp]$IDAA, 1, 0.6)
    ],
    aes(y = counts, col = IDAA),
    shape = 4
    #aes(shape = ATG, colour = endpoint6_s, y = counts)
  ) +
  scale_y_continuous(
    breaks = log(c(5, 25, 100, 500, 1500)),
    labels = c(5, 25, 100, 500, 1500)
  ) +
  facet_grid(facets = cell_type ~ endpoint6_s) +
  theme_minimal() +
  theme(legend.position = "none")

p
library(plotly)
gg <- plotly::ggplotly(p, tooltip = "IDAA")
plotly::highlight(gg, on = "plotly_hover", dynamic = TRUE)

# .. then back to Fine-Gray coding..

d <- highlight_key(txhousing, ~city)
p <- ggplot(d, aes(date, median, group = city)) + geom_line() + facet_wrap(~ month)
gg <- ggplotly(p, tooltip = "city")
highlight(gg, dynamic = TRUE, color = "red")



# Full attempt ------------------------------------------------------------

