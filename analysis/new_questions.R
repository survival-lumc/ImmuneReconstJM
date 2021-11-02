# Load object
tar_load(
  c(
    long_submodels,
    datasets,
    dat_merged,
    reference_values
  )
)

dat_wide <- datasets$wide
dat_long <- datasets$long

theme_set(theme_bw(base_size = 14))


# Try time-dep model  -----------------------------------------------------


dat_long_evfree <- dat_long[endpoint6_s == "cens"]
dat_wide_evfree <- dat_wide[endpoint6_s == "cens"]

dat_long_evfree[, ':=' (
  t_from_dli = round((intSCT2_5 >= uDLI) * (intSCT2_5 - uDLI), digits = 6),
  tdep_dli_ind = as.numeric(intSCT2_5 >= uDLI)
)]

lmeFit_tdep <- lme(
  fixed = CD3_abs_log ~ (ns(intSCT2_5, df = 1) + ns(t_from_dli, df = 1)) * ATG,
  #random = ~ intSCT2_5 + t_from_dli + tdep_dli_ind | IDAA,
  random = list(IDAA = pdDiag(~ ns(intSCT2_5, df = 2) + ns(t_from_dli, df = 1))),
  control = lmeControl(opt = "optim", msMaxIter = 200),
  data = data.frame(dat_long_evfree)
)

# Make raw plots

lmeFit_tdep <- lme(
  fixed = CD3_abs_log ~ ns(intSCT2_5, df = 3) * ATG,
  #random = ~ intSCT2_5 + t_from_dli + tdep_dli_ind | IDAA,
  random = list(IDAA = pdDiag(~ ns(intSCT2_5, df = 3))),
  control = lmeControl(opt = "optim", msMaxIter = 200),
  data = data.frame(dat_long_evfree)
)

summary(lmeFit_tdep)

newdat <- expand.grid(
  "ATG" = levels(dat_long_evfree$ATG),
  "VCMVPAT_pre" = levels(dat_long_evfree$CMV_PD),
  "intSCT2_5" = seq(0.1, 24, by = 0.1),
  "t_from_dli" = 0
)

dli_times <- seq(3, 10)
names(dli_times) <- paste0("DLI given at t = ", dli_times)

dli_preds <- lapply(dli_times, function(t_dli) {

  newdat$t_from_dli <- (newdat$intSCT2_5 >= t_dli) * (newdat$intSCT2_5 - t_dli)
  newdat$tdep_dli_ind <- as.numeric(newdat$intSCT2_5 >= t_dli)

  cbind.data.frame(
    newdat,
    "preds" = predict(lmeFit_tdep, newdata = newdat, level = 0L),
    "t_dli" = t_dli
  )
})

preds_df <- rbindlist(dli_preds, idcol = "DLI_time")
preds_df[, "DLI_time" := as.factor(DLI_time)]

ggplot(data = preds_df, aes(x = intSCT2_5, y = preds, col = ATG)) +
  #geom_vline(aes(xintercept = t_dli), linetype = "dotted") +
  geom_line(data = preds_df[preds_df$tdep_dli_ind == 0, ], size = 1.25) +
  geom_line(data = preds_df[preds_df$tdep_dli_ind == 1, ], size = 1.25) +
  scale_y_continuous(
    breaks = log(c(1, 5, 25, 100, 500, 2000)),
    labels = c(1, 5, 25, 100, 500, 2000)
  ) +
  labs(x = "Time since HSCT (months)", y = "CD3 cell counts") #+
 # facet_wrap(~ DLI_time)



# Testing subset with DLI -------------------------------------------------


dat_long_dli <- dat_long[uDLI_s == "uDLI"]
dat_wide_dli <- dat_wide[uDLI_s == "uDLI"]

table(dat_wide_dli$endpoint6_s)

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

# Model with DLI as terminating time-point?
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


# JM fitting --------------------------------------------------------------



tmat <- trans.comprisk(K = 2, names = levels(dat_wide_prepped$cr_new))
tmat

covs <- c("CMV_PD", "hirisk", "ATG", paste0(c("CD4", "CD8", "CD19", "CD3", "NK"), "_last"))
dat_wide_prepped[, ':=' (
  GVHD_ind = as.numeric(cr_new == "GVHD"),
  REL_NRF_ind = as.numeric(cr_new == "REL_NRF")
)]

msdat <- msprep(
  time = c(NA, "endpoint6", "endpoint6"),
  status = c(NA, "GVHD_ind", "REL_NRF_ind"),
  data = data.frame(dat_wide_prepped),
  trans = tmat,
  keep = covs,
  id = "IDAA"
)

msdat_expand <- mstate::expand.covs(
  msdat,
  covs,
  append = TRUE,
  longnames = FALSE
)

mod_comp <- coxph(
  Surv(Tstart, Tstop, status) ~
    ATG.1 + CD3_last.1 + #GVHD
    ATG.2 + hirisk.2 + CD3_last.2 +  # REL and NRF
    strata(trans),
  cluster = IDAA,
  model = TRUE,
  x = TRUE,
  data = msdat_expand
)

summary(mod_comp)

# Model longitudinal
lmeFit <- lme(
  fixed = CD3_abs_log ~ ns(intSCT2_5, 3) * ATG + CMV_PD,
  #random = ~ ns(intSCT2_5, 2) | IDAA,
  random = list(IDAA = pdDiag(~ ns(intSCT2_5, 3))),
  control = lmeControl(opt = "optim", msMaxIter = 200),
  data = data.frame(dat_long_prepped)
)

summary(lmeFit)


# Try joint model
jm_fit <- jointModel(
  lmeObject = lmeFit,
  survObject = mod_comp,
  CompRisk = TRUE,
  timeVar = "intSCT2_5",
  method = "spline-PH-aGH",
  interFact = list("value" = ~ strata(trans) - 1),
  control = list("iter.EM" = 300)
)

summary(jm_fit)


newdat_jm <- expand.grid(
  "ATG" = levels(dat_long_evfree$ATG),
  "CMV_PD" = levels(dat_long_evfree$CMV_PD),
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
    alpha = 0.5,
    col = NA
  ) +
  geom_line(size = 1.5, alpha = 0.75) +
  #geom_text(aes(label = label), hjust = 0, na.rm = TRUE, fontface = "bold") +
  facet_wrap(facets = ~ CMV_PD) +
  labs(x = "Time since alloHCT (months)", y = "CD8 cell counts") +
  scale_y_continuous(
    breaks = log(c(5, 25, 100, 500, 1500)),
    labels = c(5, 25, 100, 500, 1500)
  ) +
  xlim(c(0, 24)) +
  geom_vline(xintercept = 6, linetype = "dashed")

