# Only with CD3 in NMA cohort

# Load objects
tar_load(c("datasets", "dat_merged", "reference_values"))
theme_set(theme_bw(base_size = 14))

dat_wide <- datasets$wide
dat_long <- datasets$long

# Other references:
# Slope and dform:
# - https://github.com/graemeleehickey/comprisk/blob/7e558447541f8734d94135cd2e509528dd53d650/rizopoulos2012.R
# - https://github.com/drizopoulos/website/blob/fb292ed3ba777e8734906c84a35340996368d272/static/courses/Int/Solutions_SACEMA_2015.R
# Tests on CR and multistate:
# - https://github.com/drizopoulos/JMbayes2/blob/29a4f72e192d7847686432501f6116b381e32764/Development/Dev_Local_GP/MS_CR/implementation_example_v001.R
# - https://github.com/drizopoulos/JMbayes2/blob/29a4f72e192d7847686432501f6116b381e32764/Development/Dev_Local_GP/MS_CR/Competing_Risks_Reproduce.R


# Raw plots ---------------------------------------------------------------


dat_evfree <- dat_long[endpoint6_s == "cens"] # also "7 days after cellular intervention" ?

dat_evfree |>
  ggplot(aes(intSCT2_5, CD3_abs_log)) +
  geom_vline(aes(xintercept = uDLI), linetype = "dashed") +
  geom_point(
    size = 1.5,
    pch = 21,
    alpha = 0.8,
    col = "#359fda",
    fill = colorspace::lighten("#359fda", amount = 0)
  ) +
  xlim(c(0, 24)) +
  facet_wrap(~ IDAA) +
  scale_y_continuous(
    breaks = log(c(1, 5, 25, 100, 500, 2000)),
    labels = c(1, 5, 25, 100, 500, 2000)
  ) +
  labs(x = "Time since HSCT (months)", y = "CD3 cell counts")


# Pre-DLI model  ----------------------------------------------------------


# Take code from here:
#https://git.lumc.nl/hp-ldw-eb/ImmuneReconstJM/-/blob/9a94053fe0dfbd3ae1b6d9798fa69584f67ed6b3/analysis/model-msprep.R

dat_long_predli <- copy(dat_long)
dat_wide_predli <- copy(dat_wide)

# Limit to first 12 months after SCT, keep only cell measures pre-endpoint
dat_long_predli[endpoint5 >= 9, ':=' (endpoint5 = 9, endpoint5_s = "cens")]
dat_long_predli <- dat_long_predli[intSCT2_5 < endpoint5]
dat_wide_predli[endpoint5 >= 9, ':=' (endpoint5 = 9, endpoint5_s = "cens")]

# Manage modified t-cell products
dat_wide_predli[endpoint_specify5 == "modified T-cell product", endpoint5_s := "cens"]
dat_long_predli[endpoint_specify5 == "modified T-cell product", endpoint5_s := "cens"]
table(dat_wide_predli$endpoint5_s)

# Prep msdata
event_names <- levels(dat_wide_predli$endpoint5_s)
tmat <- trans.comprisk(K = 4, names = c("event_free", event_names[-1]))

ind_cols <- paste0("ind_", event_names[-1])
dat_wide_predli[, (ind_cols) := lapply(
  event_names[-1], function(col) as.numeric(endpoint5_s == col)
)]
covs <- c("SCT_May2010", "hirisk", "ATG", "CMV_PD")

msdat_predli <- msprep(
  time = c(NA, rep("endpoint5", times = length(ind_cols))),
  status = c(NA, ind_cols),
  data = data.frame(dat_wide_predli),
  trans = tmat,
  keep = covs,
  id = "IDAA"
)

msdat_predli_expand <- expand.covs(
  msdat_predli,
  covs,
  append = TRUE,
  longnames = FALSE
)

coxCRfit <- coxph(
  Surv(time, status) ~
    SCT_May2010.1 + hirisk.1 + # cellular intervention
    ATG.2 + hirisk.2 + # relapse
    ATG.3 + # NRF other; then no covars for gvhd (only current val in JM)
    strata(trans) + cluster(IDAA),
  data = msdat_predli_expand,
  x = TRUE,
  model = TRUE
)

basehaz(coxCRfit, centered = FALSE) |>
  ggplot(aes(time, hazard)) +
  geom_step(size = 1) +
  facet_wrap(~ strata, labeller = labeller(strata = c(
    "trans=1" = "1: DLI",
    "trans=2" = "2: Relapse",
    "trans=3" = "3: NRF",
    "trans=4" = "4: GVHD"
  ))) +
  labs(x = "Time since alloSCT (months)", y = "Baseline hazard")

lmeFit_df4 <- lme(
  CD4_abs_log ~ ns(intSCT2_5, 4) * ATG + CMV_PD,
  random = list(IDAA = pdDiag(form = ~ ns(intSCT2_5, 4))),
  #random = ~ ns(intSCT2_5, 2) | IDAA,
  control = lmeControl(opt = 'optim'),
  data = dat_long_predli
)
lmeFit_df3 <- lme(
  CD4_abs_log ~ ns(intSCT2_5, 3) * ATG + CMV_PD,
  random = list(IDAA = pdDiag(form = ~ ns(intSCT2_5, 3))),
  #random = ~ ns(intSCT2_5, 2) | IDAA,
  control = lmeControl(opt = 'optim'),
  data = dat_long_predli,
)
anova(update(lmeFit_df3, method = "ML"), update(lmeFit_df, method = "ML"))

# Predict mixed model
cbind.data.frame(newdat_jm, "pred" = predict(lmeFit, level = 0L, newdata = newdat_jm)) |>
  ggplot(aes(x = intSCT2_5, y = pred, group = ATG, col = ATG)) +
  geom_line(size = 1.5, alpha = 0.75) +
  #geom_text(aes(label = label), hjust = 0, na.rm = TRUE, fontface = "bold") +
  facet_wrap(facets = ~ CMV_PD) +
  labs(x = "Time since alloHCT (months)", y = "CD4 cell counts") +
  scale_y_continuous(
    breaks = log(c(5, 25, 100, 500, 1500)),
    labels = c(5, 25, 100, 500, 1500)
  ) +
  xlim(c(0, 9)) +
  geom_vline(xintercept = 6, linetype = "dashed")


# Plot the raw data
dat_long_predli |>
  ggplot(aes(intSCT2_5, CD3_abs_log)) +
  geom_line(aes(group = IDAA), alpha = 0.8, col = "gray") +
  geom_smooth(se = FALSE, method = "lm", formula = y ~ ns(x, 4)) +
  facet_grid(ATG ~ CMV_PD) +
  scale_y_continuous(
    breaks = log(c(5, 25, 100, 500, 1500)),
    labels = c(5, 25, 100, 500, 1500)
  )


jm_predli <- jointModel(
  lmeObject = lmeFit_df4,
  survObject = coxCRfit,
  timeVar = "intSCT2_5",
  method = "spline-PH-aGH",# or piecewise
  CompRisk = TRUE,
  interFact = list("value" = ~ strata(trans) - 1),
  iter.EM = 200
)

summary(jm_predli)

newdat_jm <- expand.grid(
  "ATG" = levels(dat_wide_predli$ATG),
  "CMV_PD" = levels(dat_wide_predli$CMV_PD),
  "intSCT2_5" = seq(0.1, 9, by = 0.1)
)

predict(
  jm_predli,
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
  xlim(c(0, 9)) #+
  #geom_vline(xintercept = 6, linetype = "dashed")


# Merge NRF other and DLI attempt  ----------------------------------------


dat_long_predli <- copy(dat_long)
dat_wide_predli <- copy(dat_wide)

# Limit to first 12 months after SCT, keep only cell measures pre-endpoint
dat_long_predli[endpoint5 >= 9, ':=' (endpoint5 = 9, endpoint5_s = "cens")]
dat_long_predli <- dat_long_predli[intSCT2_5 < endpoint5]
dat_wide_predli[endpoint5 >= 9, ':=' (endpoint5 = 9, endpoint5_s = "cens")]

# Combine endpoints
dat_long_predli[, cr_new := factor(
  fcase(
    endpoint5_s %in% c("cell_interv", "NRF_other"), "cell_or_NRF",
    endpoint5_s == "REL", "REL",
    endpoint5_s == "NRF_gvhd", "GVHD",
    endpoint5_s == "cens", "event_free"
  ), levels = c("event_free", "GVHD", "REL", "cell_or_NRF")
)]
dat_wide_predli[, cr_new := factor(
  fcase(
    endpoint5_s %in% c("cell_interv", "NRF_other"), "cell_or_NRF",
    endpoint5_s == "REL", "REL",
    endpoint5_s == "NRF_gvhd", "GVHD",
    endpoint5_s == "cens", "event_free"
  ), levels = c("event_free", "GVHD", "REL", "cell_or_NRF")
)]

table(dat_wide_predli$cr_new)

# Prep msdata
event_names <- levels(dat_wide_predli$cr_new)
tmat <- trans.comprisk(K = 3, names = c("event_free", event_names[-1]))

ind_cols <- paste0("ind_", event_names[-1])
dat_wide_predli[, (ind_cols) := lapply(
  event_names[-1], function(col) as.numeric(cr_new == col)
)]
covs <- c("SCT_May2010", "hirisk", "ATG", "CMV_PD")

msdat_predli <- msprep(
  time = c(NA, rep("endpoint5", times = length(ind_cols))),
  status = c(NA, ind_cols),
  data = data.frame(dat_wide_predli),
  trans = tmat,
  keep = covs,
  id = "IDAA"
)

msdat_predli_expand <- expand.covs(
  msdat_predli,
  covs,
  append = TRUE,
  longnames = FALSE
)

coxCRfit <- coxph(
  Surv(time, status) ~ #.1 is gvhd: no events
    ATG.2 + hirisk.2 + # relapse
    ATG.3 + SCT_May2010.3 + hirisk.3 + # NRF other/cell_intervention
    strata(trans) + cluster(IDAA),
  data = msdat_predli_expand,
  x = TRUE,
  model = TRUE
)

lmeFit <- lme(
  CD4_abs_log ~ ns(intSCT2_5, 4) * ATG + CMV_PD,
  random = list(IDAA = pdDiag(form = ~ ns(intSCT2_5, 4))),
  #random = ~ ns(intSCT2_5, 3) | IDAA,
  control = lmeControl(opt = 'optim'),
  data = dat_long_predli
)

jm_predli <- jointModel(
  lmeObject = lmeFit,
  survObject = coxCRfit,
  timeVar = "intSCT2_5",
  method = "spline-PH-aGH",# or piecewise
  CompRisk = TRUE,
  interFact = list("value" = ~ strata(trans) - 1),
  iter.EM = 200
)

summary(jm_predli)

#Plot
newdat_jm <- expand.grid(
  "ATG" = levels(dat_wide_predli$ATG),
  "CMV_PD" = levels(dat_wide_predli$CMV_PD),
  "intSCT2_5" = seq(0.1, 9, by = 0.1)
)

predict(
  jm_predli,
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
  labs(x = "Time since alloHCT (months)", y = "CD4 cell counts") +
  scale_y_continuous(
    breaks = log(c(5, 25, 100, 500, 1500)),
    labels = c(5, 25, 100, 500, 1500)
  ) +
  xlim(c(0, 20)) +
  geom_vline(xintercept = 6, linetype = "dashed")

# Post-DLI model (clock reset) --------------------------------------------

dat_long_dli <- dat_long[uDLI_s == "uDLI"]
dat_wide_dli <- dat_wide[uDLI_s == "uDLI"]
table(dat_wide_dli$endpoint6_s)


# #Add predicted value to dataset
newdat_predval <- cbind.data.frame(
  "ATG" = dat_wide_dli$ATG,
  "CMV_PD" = dat_wide_dli$CMV_PD,
  "intSCT2_5" = dat_wide_dli$uDLI
)

dat_wide_dli[, "CD3_predli" := predict(
  jm_predli,
  newdata = newdat_predval,
  type = "Marginal",
  idVar = "IDAA",
  returnData = FALSE,
  interval = "none"
)]

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

dat_long_dli[uDLI < intSCT2_5] |>
  ggplot(aes(intSCT2_5 - uDLI, CD3_abs_log)) +
  geom_line(aes(group = IDAA, col = ATG), size = 1.25, alpha = 0.7) +
  facet_grid(ATG ~ CMV_PD) +
  scale_y_continuous(
    breaks = log(c(1, 5, 25, 100, 500, 2000)),
    labels = c(1, 5, 25, 100, 500, 2000)
  ) +
  xlim(c(0, 18)) +
  labs(x = "Time since DLI", y = "CD3 cell counts")


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

# Do from time zero

table(dat_wide_dli$cr_new)

# First: predict current val at DLI time using previous model
#...
# dat_wide_prepped is after adding predicted value at DLI
dat_wide_prepped <- dat_wide_dli

dat_long_prepped <- dat_long_dli[uDLI < intSCT2_5]
dat_wide_prepped <- dat_wide_prepped[IDAA %in% unique(dat_long_prepped$IDAA)]

# Then go on


tmat <- trans.comprisk(K = 2, names = c("DLI", "GVHD", "REL_NRF"))
tmat

covs <- c("CMV_PD", "hirisk", "ATG",
          "DLI_type", "CD3_predli"
          #paste0(c("CD4", "CD8", "CD19", "CD3", "NK"), "_last")
          )

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

# Covar contrasts
table(dat_wide_prepped$DLI_type, dat_wide_prepped$cr_new)
table(dat_wide_prepped$ATG, dat_wide_prepped$cr_new)
table(dat_wide_prepped$hirisk, dat_wide_prepped$cr_new)

# Model with DLI_type.2 complains
mod_comp <- coxph(
  Surv(time, status) ~
    DLI_type.1 + #CD3_predli.1  + # GVHD  # + ATG.1 + CD3_last.1 +
    hirisk.2 + #CD3_predli.2 + # REL and NRF # + ATG.2 + CD3_last.2 +
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
  control = list("iter.EM" = 200)#, "lng.in.kn" = 4)
)

summary(jm_fit)


newdat_jm <- expand.grid(
  "ATG" = levels(dat_long_prepped$ATG),
  "CMV_PD" = levels(dat_long_prepped$CMV_PD),
  "intSCT2_5" = seq(6, 18, by = 0.1)
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
  #facet_wrap(facets = ~ CMV_PD) +
  labs(x = "Time since alloSCT (months)", y = "CD3 cell counts") +
  scale_y_continuous(
    breaks = log(c(5, 25, 100, 500, 1500)),
    labels = c(5, 25, 100, 500, 1500)
  ) +
  xlim(c(0, 18)) +
  geom_vline(xintercept = 6, linetype = "dashed")



# With slope --------------------------------------------------------------


dform <- list(
  fixed = ~ 0 + dns(intSCT2_5, 2),
  random = ~ 0 + dns(intSCT2_5, 2),
  indFixed = c(2, 3),
  indRandom = c(2, 3)
)

jm_fit_slopes <- jointModel(
  lmeObject = lmeFit,
  survObject = mod_comp,
  CompRisk = TRUE,
  timeVar = "intSCT2_5",
  method = "spline-PH-aGH",
  derivForm = dform,
  interFact = list(
    #"value" = ~ strata(trans) - 1,
    "slope" = ~ strata(trans) - 1
  ),
  parameterization = "slope",
  control = list("iter.EM" = 200)#, "lng.in.kn" = 4)
)

summary(jm_fit)
summary(jm_fit_slopes)

fitted(jm_fit_slopes, process = "Longitudinal", type = "Slope")

# Test with slope param ---------------------------------------------------


lmeFit <- lme(
  fixed = CD3_abs_log ~ intSCT2_5 + ATG,
  #random = ~ ns(intSCT2_5, 2) | IDAA,   #intSCT2_5
  random = list(IDAA = pdDiag(~ intSCT2_5)), #ns(intSCT2_5, 2)
  control = lmeControl(opt = "optim", msMaxIter = 200),
  data = data.frame(dat_long_prepped)
)

dform <- list(
  fixed = ~ 1,
  random = ~ 1,
  indFixed = c(2),
  indRandom = c(2)
)

dform <- list(
  fixed = ~ 0 + dns(intSCT2_5, 2),
  random = ~ 0 + dns(intSCT2_5, 2),
  indFixed = c(2, 3),
  indRandom = c(2, 3)
)

# Or in msdat directly
mod_comp$model$trans1 <- as.numeric(mod_comp$model$`strata(trans)` == "trans=1")
mod_comp$model$trans2 <- as.numeric(mod_comp$model$`strata(trans)` == "trans=2")

jm_fit_slopes <- jointModel(
  lmeObject = lmeFit,
  survObject = mod_comp,
  CompRisk = TRUE,
  timeVar = "intSCT2_5",
  method = "spline-PH-aGH",
  derivForm = dform,
  interFact = list(
    #"value" = ~ strata(trans) - 1,
    "slope" = ~ strata(trans) - 1
  ),
  parameterization = "slope",
  control = list("iter.EM" = 200)#, "lng.in.kn" = 4)
)

summary(jm_fit_slopes)

summary(jm_slopes)

jm_fit_both_gvhd <- jointModel(
  lmeObject = lmeFit,
  survObject = mod_comp,
  CompRisk = TRUE,
  timeVar = "intSCT2_5",
  method = "spline-PH-aGH",
  derivForm = dform,
  interFact = list(
    "value" = ~ trans1 + trans2 - 1,
    "slope" = ~ trans1 + trans2 - 1
  ),
  parameterization = "both",
  control = list("iter.EM" = 200)#, "lng.in.kn" = 4)
)

summary(jm_fit_both_gvhd)

summary(jm_fit)$AIC
summary(jm_fit_both_gvhd)$AIC


# Try with jmbayes2


msdat_expand$trans1 <- as.numeric(msdat_expand$trans == 1)
msdat_expand$trans2 <- as.numeric(msdat_expand$trans == 2)

jm_currval <- jm(
  Surv_object = mod_comp,
  Mixed_objects = lmeFit,
  data_Surv = msdat_expand,
  time_var = "intSCT2_5",
  id_var = "IDAA",
  functional_forms = ~ value(CD3_abs_log):trans # no need for -1?
)

jm_ideal <- jm(
  Surv_object = mod_comp,
  Mixed_objects = lmeFit,
  data_Surv = msdat_expand,
  time_var = "intSCT2_5",
  id_var = "IDAA",
  functional_forms = ~ value(CD3_abs_log):(trans1 + trans2) +
    slope(CD3_abs_log):trans1 - 1
)

jm_slopes <- jm(
  Surv_object = mod_comp,
  Mixed_objects = lmeFit,
  data_Surv = msdat_expand,
  time_var = "intSCT2_5",
  id_var = "IDAA",
  #functional_forms = ~ slope(CD3_abs_log) + slope(CD3_abs_log):trans2
  functional_forms = ~ slope(CD3_abs_log):strata(trans)
)

jm_fit_slopes |>  summary()
jm_slopes
jm_fit_slopes$coefficients$Dalpha
coefficients(jm_slopes)

jm_fit_slopes <- jointModel(
  lmeObject = lmeFit,
  survObject = mod_comp,
  CompRisk = TRUE,
  timeVar = "intSCT2_5",
  method = "spline-PH-aGH",
  derivForm = dform,
  interFact = list("slope" = ~ trans1 + trans2 - 1),
  parameterization = "slope",
  control = list("iter.EM" = 200)#, "lng.in.kn" = 4)
)

summary(jm_slopes)
jm_fit_slopes |>  summary()

predict(object = jm_slopes, type = "mean_subject")


pat <- data.frame(dat_long_prepped[IDAA == "6271"])
t0 <- 7
#pato <- pat[pat$end < t0, ]

mscopy <- msdat_expand
mscopy

pato$status <- 0
pato$intSCT2_5 <- t0
pato$DLI_type.1 <- pato$DLI_type
pato$hirisk.2 <- pato$hirisk
pato$trans1 <- 1
pato$trans2 <- 1

predict(jm_slopes, newdata = pato, type = "mean_subject")


# Try with one data for everything
dat_long_prepped



dat_long_prepped[, ':=' (
  GVHD_ind = as.numeric(cr_new == "GVHD"),
  REL_NRF_ind = as.numeric(cr_new == "REL_NRF")
)]

msdat <- msprep(
  time = c(NA, "endpoint6", "endpoint6"),
  status = c(NA, "GVHD_ind", "REL_NRF_ind"),
  start = list(state = rep(1, nrow(dat_long_prepped)), time = dat_long_prepped$uDLI),
  data = data.frame(dat_long_prepped),
  trans = tmat,
  #keep = covs,
  id = "IDAA"
)

#???


# Try with JMbayes --------------------------------------------------------

# Tutorials:
#- https://github.com/SaraBaart/Joint-Model-Tutorial/blob/main/Analyses-JM-Tutorial.R
# - https://github.com/drizopoulos/JMbayes/issues/39
library(JMbayes)
multMixedFit <- mvglmer(list(CD3_abs_log ~ intSCT2_5 + ATG  +
                               (intSCT2_5 | IDAA)),
                        data = dat_long_prepped, families = list(gaussian))

# Fit the joint model
JM1 <- mvJointModelBayes(
  multMixedFit, mod_comp,
  timeVar = "intSCT2_5",
  Interactions = list("CD3_abs_log" = strata(trans) - 1),
  multiState = TRUE,
  idVar_MultiState = "IDAA",
  #control = list(equal.strata.knots = TRUE,
   #              equal.strata.bound.knots = TRUE),
  data_MultiState = msdat_expand#,
  # Formulas = list(
  #   "CD3_abs_log" = list(
  #     fixed = ~ 1,
  #     random = ~ 1,
  #     indFixed = c(2),
  #     indRandom = c(2),
  #     name = "slope"
  #   )
  # )
)




# USE SAME DATA FOR BOTH...



# Attempt with same data --------------------------------------------------





