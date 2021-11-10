# Load objects
tar_load(c("datasets", "dat_merged", "reference_values"))
theme_set(theme_bw(base_size = 14))

dat_wide <- datasets$wide
dat_long <- datasets$long


# Raw plots event free NMA ------------------------------------------------


dat_evfree <- dat_long[endpoint6_s == "cens"] # also "7 days after cellular intervention" ?

dat_evfree |>
  ggplot(aes(intSCT2_5, CD3_abs_log)) +
  geom_point(
    size = 1.5,
    pch = 21,
    alpha = 0.8,
    col = "#359fda",
    fill = colorspace::lighten("#359fda", amount = 0)
  ) +
  facet_wrap(~ IDAA) +
  scale_y_continuous(
    breaks = log(c(1, 5, 25, 100, 500, 2000)),
    labels = c(1, 5, 25, 100, 500, 2000)
  ) +
  labs(x = "Time since HSCT (months)", y = "CD3 cell counts")


# Model in subset with DLI ------------------------------------------------


dat_long_dli <- dat_long[uDLI_s == "uDLI"]
dat_wide_dli <- dat_wide[uDLI_s == "uDLI"]

# This + below justifies composite REL and NRF endpoint
table(dat_wide_dli$endpoint6_s)

# First look at raw trajectories
dat_long_dli[uDLI < intSCT2_5] |>
  ggplot(aes(intSCT2_5, CD3_abs_log, col = ATG, group = IDAA, fill = ATG)) +
  geom_point(
    size = 2.5,
    pch = 21,
    alpha = 0.8
  ) +
  geom_line() +
  facet_grid(CMV_PD ~ endpoint6_s) +
  scale_y_continuous(
    breaks = log(c(1, 5, 25, 100, 500, 2000)),
    labels = c(1, 5, 25, 100, 500, 2000)
  ) +
  xlim(c(0, 18)) +
  labs(x = "Time since HSCT (months)", y = "CD3 cell counts")

# Needs a simpler longitudinal model, linear or spline with low df make be fine
dat_long_dli[uDLI < intSCT2_5] |>
  ggplot(aes(intSCT2_5, CD3_abs_log)) +
  geom_line(aes(group = IDAA, col = ATG), size = 1.25, alpha = 0.7) +
  facet_grid(ATG ~ CMV_PD) +
  scale_y_continuous(
    breaks = log(c(1, 5, 25, 100, 500, 2000)),
    labels = c(1, 5, 25, 100, 500, 2000)
  ) +
  xlim(c(0, 18)) +
  labs(x = "Time since HSCT (months)", y = "CD3 cell counts")


# Join REL and NRF into composite endpoint
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

table(dat_wide_dli$cr_new)

#... use only measurements after DLI, keep last measurement pre or
# at DLI as baseline covar (in mixed or survival submodel?)

# Number of measurements post DLI...
dat_long_dli[uDLI < intSCT2_5][, .(.N), by = IDAA][["N"]] |> table() # hist()
dat_long_dli[uDLI < intSCT2_5][, c(
  "IDAA", "CD4_abs_log", "cr_new", "uDLI", "intSCT2_5", "endpoint6"
)] |> View()

# Last measurement pre_dli
last_measurement <- dat_long_dli[uDLI > intSCT2_5, .SD[.N], by = IDAA]
last_measurement <- last_measurement[, c(
  "IDAA",
  paste0(c("CD4", "CD8", "CD19", "CD3", "NK"), "_abs_log"),
  "intSCT2_5",
  "cr_new",
  "endpoint6"
)]
data.table::setnames(
  x = last_measurement,
  old = paste0(c("CD4", "CD8", "CD19", "CD3", "NK"), "_abs_log"),
  new = paste0(c("CD4", "CD8", "CD19", "CD3", "NK"), "_last")
)
last_measurement

dat_wide_prepped <- merge(
  x = last_measurement,
  y = dat_wide_dli,
  by = c("IDAA", "cr_new", "endpoint6")
)

dat_long_prepped <- dat_long_dli[uDLI < intSCT2_5]
dat_wide_prepped <- dat_wide_prepped[IDAA %in% unique(dat_long_prepped$IDAA)]


# Joint model from DLI ----------------------------------------------------


tmat <- trans.comprisk(K = 2, names = c("DLI", "GVHD", "REL_NRF"))
tmat

covs <- c("CMV_PD", "hirisk", "ATG", "DLI_type", paste0(c("CD4", "CD8", "CD19", "CD3", "NK"), "_last"))
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
msdat |>  View()

msdat_expand <- mstate::expand.covs(
  msdat,
  covs,
  append = TRUE,
  longnames = FALSE
)

table(dat_wide_prepped$DLI_type, dat_wide_prepped$cr_new)
table(dat_wide_prepped$ATG, dat_wide_prepped$cr_new)
table(dat_wide_prepped$hirisk, dat_wide_prepped$cr_new)

# Model with DLI_type.2 complains
mod_comp <- coxph(
  Surv(Tstart, Tstop, status) ~
    DLI_type.1 + CD3_last.1  + # GVHD  # + ATG.1 + CD3_last.1 +
    CD3_last.2 + # REL and NRF # + ATG.2 + CD3_last.2 +
    strata(trans),
  cluster = IDAA,
  model = TRUE,
  x = TRUE,
  data = msdat_expand
)

summary(mod_comp)

# Model longitudinal, interaction ATG?
lmeFit <- lme(
  fixed = CD3_abs_log ~ ns(intSCT2_5, 2) + ATG,
  #random = ~ ns(intSCT2_5, 2) | IDAA,   #intSCT2_5
  random = list(IDAA = pdDiag(~ ns(intSCT2_5, 2))), #ns(intSCT2_5, 2)
  control = lmeControl(opt = "optim", msMaxIter = 200),
  data = data.frame(dat_long_prepped)
)

summary(lmeFit)
summary(lmeFit)$AIC

# Try joint model
jm_fit <- jointModel(
  lmeObject = lmeFit,
  survObject = mod_comp,
  CompRisk = TRUE,
  timeVar = "intSCT2_5",
  method = "spline-PH-aGH",
  interFact = list("value" = ~ strata(trans) - 1),
  control = list("iter.EM" = 500)#, "lng.in.kn" = 4)
)

summary(jm_fit)


newdat_jm <- expand.grid(
  "ATG" = levels(dat_long_prepped$ATG),
  "CMV_PD" = levels(dat_long_prepped$CMV_PD),
  "intSCT2_5" = seq(6, 24, by = 0.1)
)

predict(
  jm_fit,
  newdata = newdat_jm,
  type = "Marginal",
  idVar = "IDAA",
  returnData = TRUE,
  interval = "confidence"
) |>
  ggplot(aes(x = intSCT2_5, y = pred, group = ATG, col = ATG)) +
  geom_ribbon(
    aes(ymin = low, ymax = upp),
    fill = "gray",
    alpha = 0.3,
    col = NA
  ) +
  geom_line(size = 1.5, alpha = 0.75) +
  #geom_text(aes(label = label), hjust = 0, na.rm = TRUE, fontface = "bold") +
  facet_wrap(facets = ~ CMV_PD) +
  labs(x = "Time since alloHCT (months)", y = "CD3 cell counts") +
  scale_y_continuous(
    breaks = log(c(5, 25, 100, 500, 1500)),
    labels = c(5, 25, 100, 500, 1500)
  ) +
  xlim(c(0, 24)) +
  geom_vline(xintercept = 6, linetype = "dashed")



