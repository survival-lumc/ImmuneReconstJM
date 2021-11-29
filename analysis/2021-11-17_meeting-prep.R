# Argumentation against current DLI model:



# Load and prep data  -----------------------------------------------------


# Load objects
tar_load(c("datasets", "dat_merged", "reference_values"))
theme_set(theme_bw(base_size = 14))

# Prepare wide and long datasets
dat_wide <- datasets$wide
dat_long <- datasets$long
dat_long_dli <- copy(dat_long)[uDLI_s == "uDLI"]
dat_wide_dli <- copy(dat_wide)[uDLI_s == "uDLI"]
table(dat_wide_dli$endpoint6_s)

# Define the competing endpoints
dat_long_dli[, cr_new := factor(
  fcase(
    endpoint6_s %in% c("REL", "NRF_other"), "REL_NRF",
    endpoint6_s == "NRF_gvhd", "GVHD",
    endpoint6_s == "cens", "event_free"
  ), levels = c("event_free", "GVHD", "REL_NRF")
)]
dat_wide_dli[, cr_new := factor(
  fcase(
    endpoint6_s %in% c("REL", "NRF_other"), "REL_NRF",
    endpoint6_s == "NRF_gvhd", "GVHD",
    endpoint6_s == "cens", "event_free"
  ), levels = c("event_free", "GVHD", "REL_NRF")
)]

# Make sure to clock-reset measurements!
dat_long_prepped <- dat_long_dli[uDLI < intSCT2_5]
dat_long_prepped[, intSCT2_5_reset := intSCT2_5 - uDLI]
dat_wide_prepped <- dat_wide_dli[IDAA %in% unique(dat_long_prepped$IDAA)]

# Example of patient with no measurements between DLI and GVHD
table(dat_wide_prepped$cr_new)
dat_long_dli[IDAA == "5875"]  |>  View()


# Prep submodels ----------------------------------------------------------


tmat <- trans.comprisk(K = 2, names = c("DLI", "GVHD", "REL_NRF"))
covs <- c("CMV_PD", "hirisk", "ATG", "DLI_type")

dat_wide_prepped[, ':=' (
  GVHD_ind = as.numeric(cr_new == "GVHD"),
  REL_NRF_ind = as.numeric(cr_new == "REL_NRF")
)]

msdat <- msprep(
  time = c(NA, "endpoint6", "endpoint6"),
  status = c(NA, "GVHD_ind", "REL_NRF_ind"),
  start = list(state = rep(1, nrow(dat_wide_prepped)), time = dat_wide_prepped$uDLI),
  data = data.frame(dat_wide_prepped),
  trans = tmat,
  keep = covs,
  id = "IDAA"
)

msdat_expand <- expand.covs(msdat, covs, append = TRUE, longnames = FALSE)

# Survival submodel (excluding pre-DLI info)
mod_comp <- coxph(
  Surv(time, status) ~
    DLI_type.1 + #CD3_predli.1  + # GVHD
    hirisk.2 + #CD3_predli.2 + # REL and NRF, hirisk is DLI_indication?
    strata(trans),
  cluster = IDAA,
  model = TRUE,
  x = TRUE,
  data = msdat_expand
)

summary(mod_comp)

# Pretend linear mixed model for now
# (cmv_pd not planned in analysis but apparently important)
lmeFit <- lme(
  fixed = CD3_abs_log ~ intSCT2_5_reset * CMV_PD, # should be DLI type and ATG
  #random = ~ intSCT2_5_reset | IDAA,
  random = list(IDAA = pdDiag(~ intSCT2_5_reset)),
  control = lmeControl(opt = "optim", msMaxIter = 200),
  data = data.frame(dat_long_prepped)
)

lmeFit <- lme(
  fixed = CD3_abs_log ~ ns(intSCT2_5_reset, 2) * CMV_PD, # should be DLI type and ATG
  random = ~ ns(intSCT2_5_reset, 2)| IDAA,
  #random = list(IDAA = pdDiag(~ intSCT2_5_reset)),
  control = lmeControl(opt = "optim", msMaxIter = 200),
  data = data.frame(dat_long_prepped)
)

summary(lmeFit)

# Plots -------------------------------------------------------------------


# Look at time since first DLI
ggplot(data = dat_long_prepped, aes(intSCT2_5_reset, CD3_abs_log)) +

  # Last measurement pre-DLI
  geom_point(
    data = dat_long_dli[
      intSCT2_5 <= uDLI & IDAA %in% unique(dat_long_prepped$IDAA),
      .SD[.N], by = IDAA
    ],
    aes(intSCT2_5 - uDLI, CD3_abs_log),
    col = "red"
  ) +
  geom_point() +
  facet_wrap(~ IDAA) +
  geom_vline(xintercept = 0, linetype = "dotted") +
  theme(legend.position = "none") +
  scale_y_continuous(
    breaks = log(c(1, 5, 25, 100, 500, 2000)),
    labels = c(1, 5, 25, 100, 500, 2000)
  ) +
  labs(x = "Time since first DLI (months)", y = "CD3 cell counts")

# Look at solely those with gvhd
ggplot(data = dat_long_prepped[cr_new == "GVHD"], aes(intSCT2_5_reset, CD3_abs_log)) +

  # Last measurement pre-DLI
  geom_point(
    data = dat_long_dli[
      intSCT2_5 <= uDLI & IDAA %in% unique(dat_long_prepped[cr_new == "GVHD"]$IDAA),
      .SD[.N], by = IDAA
    ],
    aes(intSCT2_5 - uDLI, CD3_abs_log),
    col = "red"
  ) +
  geom_point() +
  facet_wrap(~ IDAA) +
  geom_vline(xintercept = 0, linetype = "dotted") +
  theme(legend.position = "none") +
  scale_y_continuous(
    breaks = log(c(1, 5, 25, 100, 500, 2000)),
    labels = c(1, 5, 25, 100, 500, 2000)
  ) +
  labs(x = "Time since first DLI (months)", y = "CD3 cell counts")


# Plot individual predictions from mixed model
long_newdat <- expand.grid("IDAA" = dat_wide_prepped$IDAA, "intSCT2_5_reset" = seq(0, 12, by = 0.1))
df_preds <- merge(long_newdat, dat_wide_prepped)
df_preds <- data.table(cbind.data.frame(df_preds, "pred" = predict(lmeFit, newdata = df_preds)))


ggplot(
  data = df_preds,#[, .SD[.N > 1], by = IDAA],
  aes(intSCT2_5_reset, pred)
) +
  #geom_point(alpha = 0.8, col = "blue") +
  geom_point(data = dat_long_prepped,
             aes(y = CD3_abs_log), size = 1) +
  geom_line(aes(y = pred, group = IDAA), size = 0.5) +
  facet_wrap(~ IDAA) +
  theme(legend.position = "none") +
  scale_y_continuous(
    breaks = log(c(1, 5, 25, 100, 500, 2000)),
    labels = c(1, 5, 25, 100, 500, 2000)
  ) +
  labs(x = "Time since first DLI (months)", y = "CD3 cell counts")

# subset with GVHD
ggplot(
  data = df_preds[cr_new == "GVHD"],#[, .SD[.N > 1], by = IDAA],
  aes(intSCT2_5_reset, pred)
) +
  # Last measurement pre-DLI
  geom_point(
    data = dat_long_dli[
      intSCT2_5 <= uDLI & IDAA %in% unique(dat_long_prepped[cr_new == "GVHD"]$IDAA),
      .SD[.N], by = IDAA
    ],
    aes(intSCT2_5 - uDLI, CD3_abs_log),
    col = "red"
  ) +
  geom_point(
    data = dat_long_prepped[cr_new == "GVHD"],
    aes(y = CD3_abs_log), size = 1
  ) +
  geom_line(aes(y = pred, group = IDAA), size = 0.5) +
  geom_vline(xintercept = 0, linetype = "dotted") +
  facet_wrap(~ IDAA) +
  theme(legend.position = "none") +
  scale_y_continuous(
    breaks = log(c(1, 5, 25, 100, 500, 2000)),
    labels = c(1, 5, 25, 100, 500, 2000)
  ) +
  labs(x = "Time since first DLI (months)", y = "CD3 cell counts")



# Marginal trajectories ---------------------------------------------------


newdat_marg <- expand.grid(
  "CMV_PD" = levels(dat_long_prepped$CMV_PD),
  "intSCT2_5_reset" = seq(0, 12, by = 0.1)
)

df_preds_marg <- cbind(
  newdat_marg,
  "pred" = predict(lmeFit, newdata = newdat_marg, level = 0L)
)

ggplot(data = df_preds_marg, aes(intSCT2_5_reset, pred)) +
  geom_line(
    data = dat_long_prepped,
    aes(intSCT2_5_reset, CD3_abs_log, group = IDAA, col = IDAA), alpha = 0.8
  ) +
  geom_line(size = 1.25) +
  facet_grid(~ CMV_PD) + # cr_new here
  theme(legend.position = "none") +
  scale_y_continuous(
    breaks = log(c(1, 5, 25, 100, 500, 2000)),
    labels = c(1, 5, 25, 100, 500, 2000)
  ) +
  labs(x = "Time since first DLI (months)", y = "CD3 cell counts")

table(dat_wide_prepped$cr_new, dat_wide_prepped$CMV_PD)
