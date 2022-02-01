tar_load(
  c(
    dat_merged,
    NMA_preDLI_datasets,
    NMA_postDLI_datasets_CD3,
    NMA_postDLI_datasets_CD4,
    NMA_postDLI_datasets_CD8,
    preDLI_CD3__jointModel_both,
    preDLI_CD4_jointModel_both,
    preDLI_CD3_jointModel_corr,
    preDLI_CD8_jointModel_both,
    postDLI_CD3_jointModel_both,
    postDLI_CD4_jointModel_both,
    postDLI_CD8_jointModel_both
  )
)

# Check

theme_set(theme_bw(base_size = 14))
dat_long <- NMA_preDLI_datasets$long
dat_wide <- NMA_preDLI_datasets$wide

newdat_jm <- expand.grid(
  "ATG" = levels(dat_long$ATG),
  "CMV_PD" = levels(dat_long$CMV_PD),
  "hirisk" = levels(dat_long$hirisk),
  "intSCT2_5" = seq(0, 6, by = 0.05)
)



# Correlated one ----------------------------------------------------------

plot_preDLI_average_trajectories(
  model = preDLI_CD3_jointModel_corr,
  newdat = newdat_jm,
  ylab = "CD3 cell counts"
)

table_mod_summary(preDLI_CD3_jointModel_corr)

dat_long[, preds_subj := fitted(
  preDLI_CD3_jointModel_corr,
  process = "Longitudinal",
  type = "Subject"
)]
# CD3 joint model ---------------------------------------------------------

plot_preDLI_average_trajectories(
  model = preDLI_CD3__jointModel_both,
  newdat = newdat_jm,
  ylab = "CD3 cell counts"
)

summary(preDLI_CD3__jointModel_both)

table_mod_summary(preDLI_CD3__jointModel_both)

# (For Hein/Liesbeth)
# - subject-specific predictions (within scope of data)
# - get last measure preDLI
dat_long[, preds_subj := fitted(
  preDLI_CD3__jointModel_both,
  process = "Longitudinal",
  type = "Subject"
)]

# Try cumhaz
data_id <- copy(preDLI_CD3__jointModel_both$data.id)
data_id[, "cumhaz" := fitted(
  preDLI_CD3__jointModel_both,
  process = "Event",
  type = "Marginal",
  #scale = "survival"
  scale = "cumulative-Hazard"
)]
data_id[, "cr_ind" := factor(
  x = rep(seq_len(3), .N/3),
  levels = seq_len(3),
  labels = c("GVHD", "REL", "Other")
)]
data_id[, "covar_group" := interaction(ATG, hirisk, CMV_PD)]
setorder(data_id, "cr_ind", "covar_group", "intSCT2_5")
data_id |>
  ggplot(aes(intSCT2_5, cumhaz)) +
  geom_step(aes(group = covar_group)) +
  facet_grid(covar_group ~ cr_ind)


preDLI_CD3__jointModel_both$data.id

# Sample third of patients
set.seed(20220119)
IDAA_subs <- sample(levels(dat_long$IDAA), replace = FALSE, size = 56) #size = 16)#size = 56)

ggplot(
  dat_long[IDAA %in% IDAA_subs], aes(intSCT2_5, CD3_abs_log)
) +
  geom_point() +
  geom_line(aes(y = preds_subj, group = IDAA)) +
  facet_wrap(~ IDAA) +
  theme(legend.position = "none")

# Now lets try slope - get an error because dataset has multiple lines
# (stratas of competing risks)
fitted(
  preDLI_CD3__jointModel_both,
  process = "Longitudinal",
  type = "Slope"
)

# See other file
dat_long[, subj_spec_slope := ff]

# Make function for slop
slope <- 6
get_int_tangent <- function(x, y, slope) {
  m <- slope * x
  int <- -x * slope + y
  return(int)
}

get_int_tangent(x = 1, y = 3, slope = 6)
dat_long[, int_tang := get_int_tangent(
  x = intSCT2_5, y = preds_subj, slope = subj_spec_slope
)]

delta_tan <- 0.10
dat_long[, ':=' (
  start_x = intSCT2_5 - delta_tan,
  start_y = int_tang + subj_spec_slope * (intSCT2_5 - delta_tan),
  end_x = intSCT2_5 + delta_tan,
  end_y = int_tang + subj_spec_slope * (intSCT2_5 + delta_tan)
)]



dat_long[IDAA %in% IDAA_subs] |>
  ggplot(aes(intSCT2_5, preds_subj, group = IDAA)) +
  geom_point(aes(y = CD3_abs_log)) +
  geom_line(size = 0.5, aes(col = )) +
  #geom_segment(aes(x = start_x, y = start_y, xend = end_x, yend = end_y), size = 1,
  #             arrow = arrow(length = unit(0.03, "npc"))) +
  facet_wrap(~ IDAA) +
  theme(legend.position = "none")
  #geom_abline(aes(intercept = int_tang, slope = subj_spec_slope, y = preds_subj), linetype = "dotted")



# try something else
df_last_slopes  <- dat_long[, .SD[.N], by = "IDAA"]

dat_long |>
  ggplot(aes(intSCT2_5, preds_subj, group = IDAA)) +
  geom_line(aes(col = IDAA), alpha = 0.75) +
  geom_segment(
    data = df_last_slopes,
    aes(x = start_x, y = start_y, xend = end_x, yend = end_y),
    size = 0.75,
    arrow = arrow(length = unit(0.03, "npc"))
  ) +
  #geom_point(data = df_last_slopes, col = "red", shape = 4) +
 #facet_wrap(~ endpoint7_s) +
  facet_grid(endpoint7_s ~ ATG) +
  theme(legend.position = "none")

# Try density of the slopes?
df_last_slopes |>
  ggplot(aes(x = endpoint7_s, y = subj_spec_slope)) +
  geom_violin() +
  geom_point(position = position_dodge2(width = 0.1))
  #geom_histogram(col = "black", fill = "lightblue") +
  #geom_vline(xintercept = c(0, 2), linetype = "dashed") +
  #facet_grid(~ , scales = "fixed") +
  theme(legend.position = "none")

# Try excalidraw for illustration slope??

# CD4 joint model ---------------------------------------------------------


plot_preDLI_average_trajectories(
  model = preDLI_CD4_jointModel_both,
  newdat = newdat_jm,
  ylab = "CD4 cell counts"
)

table_mod_summary(preDLI_CD4_jointModel_both)



# CD8 joint model ---------------------------------------------------------


plot_preDLI_average_trajectories(
  model = preDLI_CD8_jointModel_both,
  newdat = newdat_jm,
  ylab = "CD8 cell counts"
)

table_mod_summary(preDLI_CD8_jointModel_both)


# Post DLI model ----------------------------------------------------------


postDLI_CD3_jointModel_both |>  summary()
NMA_postDLI_datasets_CD3$long |>
  ggplot(aes(intSCT2_5_reset, CD3_abs_log, col = IDAA, group = IDAA)) +
  geom_line() +
  facet_grid(earlylow_DLI ~ sec_endpoint_s) +
  theme(legend.position = "none")

dat_long_CD3 <- NMA_postDLI_datasets_CD3$long
dat_long_CD3[, preds_subj := fitted(
  postDLI_CD3_jointModel_both,
  process = "Longitudinal",
  type = "Subject"
)]

# Fit single values
ggplot(
  dat_long_CD3, aes(intSCT2_5_reset, CD3_abs_log)
) +
  geom_point() +
  geom_line(aes(y = preds_subj, group = IDAA)) +
  facet_wrap(~ IDAA) +
  theme(legend.position = "none")

# Marginal
newdat_postDLI <- expand.grid(
  "ATG" = levels(dat_long_CD3$ATG),
  "earlylow_DLI" = levels(dat_long_CD3$earlylow_DLI),
  #"hirisk" = levels(dat_long_CD3$hirisk),
  "intSCT2_5_reset" = seq(0, 12, by = 0.1)
)


dat_preds_post <- predict(
  postDLI_CD3_jointModel_both,
  newdata = newdat_postDLI,
  type = "Marginal",
  idVar = "IDAA",
  returnData = TRUE,
  interval = "confidence"
)

ggplot(
  data = dat_preds_post,
  aes(
    x = intSCT2_5_reset,
    y = pred,
    group = interaction(earlylow_DLI, ATG),
    col = earlylow_DLI
  )
) +
  geom_ribbon(
    aes(ymin = low, ymax = upp),
    fill = "gray",
    alpha = 0.3,
    col = NA
  ) +
  geom_line(
    aes(linetype = ATG),
    size = 1.5,
    alpha = 0.75
  ) +
  # add repel labels
  #facet_wrap(facets = ~ ATG) +
  labs(x = "Time since alloHCT (months)", y = "CD3 counts") #+
  scale_y_continuous(
    breaks = log(c(5, 25, 100, 500, 1500)),
    labels = c(5, 25, 100, 500, 1500)
  ) +
  xlim(c(0, 6))

# Extra functions ---------------------------------------------------------


plot_preDLI_average_trajectories <- function(model,
                                             newdat,
                                             ylab) {

  dat_preds <- predict(
    model,
    newdata = newdat,
    type = "Marginal",
    idVar = "IDAA",
    returnData = TRUE,
    interval = "confidence"
  )

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
      alpha = 0.3,
      col = NA
    ) +
    geom_line(
      aes(linetype = hirisk),
      size = 1.5,
      alpha = 0.75
    ) +
    # add repel labels
    facet_wrap(facets = ~ CMV_PD) +
    labs(x = "Time since alloHCT (months)", y = ylab) +
    scale_y_continuous(
      breaks = log(c(5, 25, 100, 500, 1500)),
      labels = c(5, 25, 100, 500, 1500)
    ) +
    xlim(c(0, 6))
}

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
    labels = c("ATG", "High risk (ITT)", "Current value", "Slope")
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
