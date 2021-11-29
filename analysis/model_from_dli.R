# Script with model from DLI (without including current val at DLI)
# .. for testing


tar_load(c("datasets", "dat_merged", "reference_values"))
theme_set(theme_bw(base_size = 14))

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


# Prep models -------------------------------------------------------------


tmat <- trans.comprisk(K = 2, names = c("DLI", "GVHD", "REL_NRF"))
tmat

covs <- c("CMV_PD", "hirisk", "ATG", "DLI_type")
          #"CD3_predli", paste0(c("CD4", "CD8", "CD19", "CD3", "NK"), "_last")

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


# Survival submodels
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

# Submodel
lmeFit <- lme(
  fixed = CD3_abs_log ~ intSCT2_5_reset + DLI_type + ATG, # cmv_pdnot plannen in analysis but apparently important
  #random = ~ ns(intSCT2_5, 2) | IDAA,   #intSCT2_5
  random = list(IDAA = pdDiag(~ intSCT2_5_reset)),
  control = lmeControl(opt = "optim", msMaxIter = 200),
  data = data.frame(dat_long_prepped)
)

summary(mod_comp)
summary(lmeFit)



# Checks on fit of longitudinal part --------------------------------------


df_preds_lme <- cbind(
  dat_long_prepped,
  "pred" = predict(lmeFit, newdata = data.frame(dat_long_prepped))
)

ggplot(
  data = df_preds_lme[, .SD[.N > 1], by = IDAA],
  aes(intSCT2_5_reset, CD3_abs_log)
) +
  geom_point(alpha = 0.8, col = "blue") +
  geom_line(aes(y = pred, group = IDAA), size = 1) +
  theme(legend.position = "none") +
  facet_wrap(~ IDAA)


# Marginal
newdat_jm <- expand.grid(
  "ATG" = levels(dat_long_prepped$ATG),
  "DLI_type" = levels(dat_long_prepped$DLI_type),
  "intSCT2_5_reset" = seq(0, 12, by = 0.1)
)

df_preds_marg <- cbind(
  newdat_jm,
  "pred" = predict(lmeFit, newdata = newdat_jm, level = 0L)
)

df_preds_marg |>
  ggplot(aes(intSCT2_5_reset, pred)) +
  geom_line(aes(col = DLI_type, group = DLI_type)) +
  facet_wrap(~ ATG)

# Slopes JM ---------------------------------------------------------------


dform <- list(
  fixed = ~ 0 + dns(intSCT2_5_reset, 3),
  random = ~ 0 + dns(intSCT2_5_reset, 3),
  indFixed = c(2:4),
  indRandom = c(2:4)
)

jm_fit_both <- jointModel(
  lmeObject = lmeFit,
  survObject = mod_comp,
  CompRisk = TRUE,
  timeVar = "intSCT2_5_reset",
  method = "spline-PH-aGH",
  derivForm = dform,
  interFact = list(
    "value" = ~ strata(trans) - 1,
    "slope" = ~ strata(trans) - 1
  ),
  parameterization = "both",
  control = list("iter.EM" = 200)
)


jm_fit_slopes <- jointModel(
  lmeObject = lmeFit,
  survObject = mod_comp,
  CompRisk = TRUE,
  timeVar = "intSCT2_5_reset",
  method = "spline-PH-aGH",
  derivForm = dform,
  interFact = list("slope" = ~ strata(trans) - 1),
  parameterization = "slope",
  control = list("iter.EM" = 200)
)

jm_fit_current <- jointModel(
  lmeObject = lmeFit,
  survObject = mod_comp,
  CompRisk = TRUE,
  timeVar = "intSCT2_5_reset",
  method = "spline-PH-aGH",
  interFact = list("value" = ~ strata(trans) - 1),
  parameterization = "value",
  control = list("iter.EM" = 200)
)

summary(jm_fit_current)
summary(jm_fit_both)
summary(jm_fit_slopes)


# With Bayesian software --------------------------------------------------


jm_slopes <- jm(
  Surv_object = mod_comp,
  Mixed_objects = lmeFit,
  data_Surv = msdat_expand,
  time_var = "intSCT2_5_reset",
  id_var = "IDAA",
  functional_forms = ~ slope(CD3_abs_log,):strata(trans)
)

jm_slopes |> coef()
jm_fit_slopes$coefficients$Dalpha
ggtraceplot(jm_slopes, "alphas")

jm_both <- jm(
  Surv_object = mod_comp,
  Mixed_objects = lmeFit,
  data_Surv = msdat_expand,
  time_var = "intSCT2_5_reset",
  id_var = "IDAA",
  functional_forms = ~ (value(CD3_abs_log) + slope(CD3_abs_log)):strata(trans),
  n_iter = 5000L,
  n_burnin = 2500L
)
jm_both
ggtraceplot(jm_both, "alphas", grid = TRUE, gridrows = 2, gridcols = 2,
            theme = "moonlight")


# Try just lines ----------------------------------------------------------



lmeFit <- lme(
  fixed = CD3_abs_log ~ intSCT2_5_reset + ATG + DLI_type,
  random = ~ intSCT2_5_reset | IDAA,   #intSCT2_5
  #random = list(IDAA = pdDiag(~ intSCT2_5_reset)),
  control = lmeControl(opt = "optim", msMaxIter = 200),
  data = data.frame(dat_long_prepped)
)

lmeFit <- lme(
  fixed = CD3_abs_log ~ ns(intSCT2_5_reset, 2) + ATG + DLI_type,
  #random = ~ ns(intSCT2_5, 2) | IDAA,   #intSCT2_5
  random = list(IDAA = pdDiag(~ ns(intSCT2_5_reset, 2))),
  control = lmeControl(opt = "optim", msMaxIter = 200),
  data = data.frame(dat_long_prepped)
)



df_preds_lme <- cbind(
  dat_long_prepped,
  "pred" = predict(lmeFit, newdata = data.frame(dat_long_prepped))
)

ggplot(
  data = df_preds_lme[, .SD[.N > 1], by = IDAA],
  aes(intSCT2_5_reset, CD3_abs_log)
) +
  geom_point(alpha = 0.8, col = "blue") +
  geom_line(aes(y = pred, group = IDAA), size = 1) +
  theme(legend.position = "none") +
  facet_wrap(~ IDAA * cr_new)

ggplot(
  data = df_preds_lme,#[, .SD[.N > 1], by = IDAA],
  aes(intSCT2_5_reset, CD3_abs_log)
) +
  #geom_point(alpha = 0.8, col = "blue") +
  geom_point(aes(y = pred, group = IDAA, col = IDAA), size = 1) +
  geom_line(aes(y = pred, group = IDAA, col = IDAA), size = 0.5) +
  facet_wrap(~ cr_new) +
  theme(legend.position = "none")


# Use prediction
new_dat_test <- copy(dat_wide_prepped)
long_test <- expand.grid("IDAA" = new_dat_test$IDAA, "intSCT2_5_reset" = seq(0, 10, by = 0.1))
newz <- merge(long_test, dat_wide_prepped)

df_predz <- cbind(
  newz,
  "pred" = predict(lmeFit, newdata = newz)
)

ggplot(
  data = df_predz,#[, .SD[.N > 1], by = IDAA],
  aes(intSCT2_5_reset, pred)
) +
  #geom_point(alpha = 0.8, col = "blue") +
  geom_point(data = dat_long_prepped,
             aes(y = CD3_abs_log, group = IDAA, col = IDAA), size = 1) +
  geom_line(aes(y = pred, group = IDAA, col = IDAA), size = 0.5) +
  facet_wrap(~ cr_new) +
  theme(legend.position = "none")

ggplot(
  data = df_predz,#[, .SD[.N > 1], by = IDAA],
  aes(intSCT2_5_reset, pred)
) +
  #geom_point(alpha = 0.8, col = "blue") +
  geom_point(data = dat_long_prepped,
             aes(y = CD3_abs_log, group = IDAA, col = IDAA), size = 1) +
  geom_line(aes(y = pred, group = IDAA, col = IDAA), size = 0.5) +
  facet_wrap(~ IDAA) +
  theme(legend.position = "none")

ggplot(
  data = df_predz,#[, .SD[.N > 1], by = IDAA],
  aes(intSCT2_5_reset, pred)
) +
  #geom_point(alpha = 0.8, col = "blue") +
  geom_point(data = dat_long_prepped,
             aes(y = CD3_abs_log, group = IDAA, col = IDAA), size = 1) +
  #geom_line(aes(y = pred, group = IDAA, col = IDAA), size = 0.5) +
  geom_smooth(data = dat_long_prepped,
              aes(y = CD3_abs_log, group = IDAA, col = IDAA),
              method = "lm", formula = y ~ x, se = FALSE) +
  facet_wrap(IDAA ~ cr_new) +
  theme(legend.position = "none")




jm_slopes <- jm(
  Surv_object = mod_comp,
  Mixed_objects = lmeFit,
  data_Surv = msdat_expand,
  time_var = "intSCT2_5_reset",
  id_var = "IDAA",
  functional_forms = ~ slope(CD3_abs_log):strata(trans)
)

jm_current <- jm(
  Surv_object = mod_comp,
  Mixed_objects = lmeFit,
  data_Surv = msdat_expand,
  time_var = "intSCT2_5_reset",
  id_var = "IDAA",
  functional_forms = ~ value(CD3_abs_log):strata(trans)
)

# Hint of current value better for rel_nrf, but slope for DLI
# But too few measurement with those early events?

jm_slopes
jm_current


dform_lin <- list(
  fixed = ~ 1,
  random = ~ 1,
  indFixed = 2,
  indRandom = 2
)


jm_lin_slopes <- jointModel(
  lmeObject = lmeFit,
  survObject = mod_comp,
  CompRisk = TRUE,
  timeVar = "intSCT2_5_reset",
  method = "spline-PH-aGH",
  derivForm = dform_lin,
  interFact = list("slope" = ~ strata(trans) - 1),
  parameterization = "slope",
  control = list("iter.EM" = 200)
)

jm_lin_slopes |>  summary()
jm_lin_slopes$coefficients$Dalpha
jm_slopes |>  coef()


msdat_expand$trans1 <- as.numeric(msdat_expand$trans == 1)
msdat_expand$trans2 <- as.numeric(msdat_expand$trans == 2)


# Jm specific

jm_specific <- jm(
  Surv_object = mod_comp,
  Mixed_objects = lmeFit,
  data_Surv = msdat_expand,
  time_var = "intSCT2_5_reset",
  id_var = "IDAA",
  functional_forms = ~ slope(CD3_abs_log):trans1 + value(CD3_abs_log):trans2,
  n_iter = 5000L,
  n_burnin = 2500L
)

jm_specific
ggtraceplot(jm_specific, "alphas", grid = TRUE, gridrows = 1, gridcols = 2,
            theme = "moonlight")
