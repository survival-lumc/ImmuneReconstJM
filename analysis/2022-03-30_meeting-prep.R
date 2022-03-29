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
dat_wide_postDLI <- NMA_postDLI_datasets$wide

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


# Make table summary, ready to feed into kable
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

  return(dt_edit)
}


# Summary pre-DLI models --------------------------------------------------

table(dat_wide$endpoint7_s)


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
  coord_cartesian(xlim = c(0, 6), ylim = c(0, log(1500))) + # add a common ylim for all
  scale_color_manual(
    labels = c("Related", "Unrelated (+ATG)"),
    values = c("brown", "darkblue")
  ) +
  scale_linetype_manual(
    labels = c("No", "Yes (high risk)"),
    values = c("solid", "dotdash")
  ) +
  theme(legend.position = "bottom")


# -- Current value summaries
table_mod_summary(mods$preDLI_JM_value_corr_CD3) %>%
  kableExtra::kbl(format = "html") %>%
  kable_paper("striped", full_width = F) %>%
  pack_rows("GVHD", 1, 3) %>%
  pack_rows("Relapse", 4, 6) %>%
  pack_rows("NRF: Other", 7, 8)

table_mod_summary(mods$preDLI_JM_value_corr_CD4) %>%
  kableExtra::kbl(format = "html") %>%
  kable_paper("striped", full_width = F) %>%
  pack_rows("GVHD", 1, 3) %>%
  pack_rows("Relapse", 4, 6) %>%
  pack_rows("NRF: Other", 7, 8)

table_mod_summary(mods$preDLI_JM_value_corr_CD8) %>%
  kableExtra::kbl(format = "html") %>%
  kable_paper("striped", full_width = F) %>%
  pack_rows("GVHD", 1, 3) %>%
  pack_rows("Relapse", 4, 6) %>%
  pack_rows("NRF: Other", 7, 8)

# -- Current value + current slope summaries
table_mod_summary(mods$preDLI_JM_both_corr_CD3) %>%
  kableExtra::kbl(format = "html") %>%
  kable_paper("striped", full_width = F) %>%
  pack_rows("GVHD", 1, 4) %>%
  pack_rows("Relapse", 5, 8) %>%
  pack_rows("NRF: Other", 9, 11)

table_mod_summary(mods$preDLI_JM_both_corr_CD4) %>%
  kableExtra::kbl(format = "html") %>%
  kable_paper("striped", full_width = F) %>%
  pack_rows("GVHD", 1, 4) %>%
  pack_rows("Relapse", 5, 8) %>%
  pack_rows("NRF: Other", 9, 11)

table_mod_summary(mods$preDLI_JM_both_corr_CD8) %>%
  kableExtra::kbl(format = "html") %>%
  kable_paper("striped", full_width = F) %>%
  pack_rows("GVHD", 1, 4) %>%
  pack_rows("Relapse", 5, 8) %>%
  pack_rows("NRF: Other", 9, 11)


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
  ggplot(aes(x = intDLI1 + 3, y = pred, group = ATG)) +
  geom_ribbon(aes(ymin = low, ymax = upp), fill = "gray", alpha = 0.5, col = NA) +
  geom_line(aes(col = ATG), size = 1.5, alpha = 0.75) +
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
  coord_cartesian(xlim = c(0, 6), ylim = c(0, log(1500))) +
  scale_color_manual(
    labels = c("Related", "Unrelated (+ATG)"),
    values = c("brown", "darkblue")
  ) +
  theme(legend.position = "bottom") +
  geom_vline(xintercept = 3, linetype = "dotted")


# -- Model summaries
table_mod_summary(mods$postDLI_JM_corr_CD3) %>%
  kableExtra::kbl(format = "html") %>%
  kable_paper("striped", full_width = F) %>%
  pack_rows("GVHD", 1, 2) %>%
  pack_rows("Relapse & other NRF", 3, 3)

table_mod_summary(mods$postDLI_JM_corr_CD4) %>%
  kableExtra::kbl(format = "html") %>%
  kable_paper("striped", full_width = F) %>%
  pack_rows("GVHD", 1, 2) %>%
  pack_rows("Relapse & other NRF", 3, 3)

table_mod_summary(mods$postDLI_JM_corr_CD8) %>%
  kableExtra::kbl(format = "html") %>%
  kable_paper("striped", full_width = F) %>%
  pack_rows("GVHD", 1, 2) %>%
  pack_rows("Relapse & other NRF", 3, 3)



# Extra raw plots ---------------------------------------------------------


# Raw counts whole post DLI cohort
dat_long_postDLI |>
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


# Trajectories post DLI
dat_long_postDLI |>
  melt.data.table(
    measure.vars = patterns("*_log$"),
    variable.name = "cell_line",
    value.name = "count"
  ) |>
  ggplot(aes(intDLI1, count, linetype = ATG, col = ATG, group = IDAA)) +
  geom_line(size = 1) +
  scale_y_continuous(
    breaks = log(c(5, 25, 100, 500, 1500)),
    labels = c(5, 25, 100, 500, 1500)
  ) +
  facet_grid(sec_endpoint2_s ~ cell_line) +
  theme_bw(base_size = 14) +
  labs(y = "Count", x = "Time since early/low DLI (months)")

# Random effects plots
GGally::ggpairs(data = as.data.frame(ranef(mods$postDLI_JM_corr_CD8))) +
  theme_bw(base_size = 14)

# Plot individual fits together for each pat
postDLI_indiv <- cbind(
  dat_long_postDLI,
  cbind.data.frame(
    "CD3_subj" = fitted(mods$postDLI_JM_corr_CD3, process = "Longitudinal", type = "Subject"),
    "CD4_subj" = fitted(mods$postDLI_JM_corr_CD4, process = "Longitudinal", type = "Subject"),
    "CD8_subj" = fitted(mods$postDLI_JM_corr_CD8, process = "Longitudinal", type = "Subject")
  )
)

# Add last measurement before endpoint
dat_long_postDLI[, "endpoint_lab" := fcase(
  sec_endpoint2_s == "cens", "Censored",
  sec_endpoint2_s == "gvhd", "GVHD",
  sec_endpoint2_s == "rel_nrf", "Relapse or \n Other NRF"
)]
dat_long_postDLI[, endpoint_lab := factor(
  endpoint_lab, levels = c("Censored", "GVHD", "Relapse or \n Other NRF")
)]

endpoints_df <- dat_long_postDLI[order(intDLI1), .SD[.N], by = IDAA]
dat_long_postDLI[, .(.N, "endp" = unique(sec_endpoint2_s)), by = IDAA] # n measurments

melt.data.table(
  data = postDLI_indiv,
  measure.vars = paste0("CD", c(3, 4, 8), "_subj"),
  variable.name = "cell_line",
  value.name = "subj_pred"
) |>
  ggplot(aes(x = intDLI1, y = subj_pred, col = cell_line,)) +
  geom_line(aes(linetype = cell_line, group = cell_line), size = 1) +
  geom_point(size = 1.5) +
  geom_label(
    data = endpoints_df,
    aes(x = sec_endpoint2 + 0.25, y = log(3), label = endpoint_lab),
    #size = 5,
    hjust = 0,
    lineheight = .8,
    inherit.aes = FALSE,
    label.size = NA
  ) +
  geom_vline(aes(xintercept = sec_endpoint2), linetype = "dashed") +
  facet_wrap(~ IDAA) +
  scale_y_continuous(
    breaks = log(c(5, 25, 100, 500, 1500)),
    labels = c(5, 25, 100, 500, 1500)
  ) +
  scale_color_brewer(palette = "Dark2") +
  coord_cartesian(xlim = c(0, 5.5), ylim = c(0, log(1500))) +
  labs(y = "Count", x = "Time since early/low DLI (months)")
