tar_load(
  c(
    dat_merged,
    NMA_preDLI_datasets,
    preDLI_CD3__jointModel_both,
    #preDLI_CD4_jointModel_both,
    #preDLI_CD8_jointModel_both
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

# NOTE ATG IS REVERSED!! (will fix)

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

# Sample third of patients
set.seed(20220119)
IDAA_subs <- sample(levels(dat_long$IDAA), replace = FALSE, size = 16)#size = 56)

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

dat_long[, ':=' (
  start_x = intSCT2_5 - 0.2,
  start_y = int_tang + subj_spec_slope * (intSCT2_5 - 0.25),
  end_x = intSCT2_5 + 0.2,
  end_y = int_tang + subj_spec_slope * (intSCT2_5 + 0.25)
)]



dat_long[IDAA %in% IDAA_subs] |>
  ggplot(aes(intSCT2_5, preds_subj, group = IDAA)) +
  geom_point(aes(y = CD3_abs_log)) +
  geom_line(size = 0.5) +
  geom_segment(aes(x = start_x, y = start_y, xend = end_x, yend = end_y), size = 1,
               arrow = arrow(length = unit(0.03, "npc"))) +
  facet_wrap(~ IDAA) +
  theme(legend.position = "none")
  #geom_abline(aes(intercept = int_tang, slope = subj_spec_slope, y = preds_subj), linetype = "dotted")


ggplot(
  dat_long[IDAA %in% IDAA_subs], aes(intSCT2_5, CD3_abs_log)
) +
  geom_point() +
  geom_line(aes(y = preds_subj)) +
  geom_line(aes(y = subj_spec_slope), col = "blue", linetype = "dashed") +
  facet_wrap(~ IDAA) +
  theme(legend.position = "none") +
  ylab("log CD3 counts/current slope") +
  geom_hline(yintercept = 0, linetype = "dashed")

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
