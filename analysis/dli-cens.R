##*********************************************************##
## Testing DLI as competing event or not (pre simulations) ##
##*********************************************************##

# Run code until JM_dat_wide

# Fit mixed model (fit simple linear one for now)
lmeFit <- nlme::lme(
  CD8_abs_log ~ ns(intSCT2_5, 3) * ATG + VCMVPAT_pre,
  random = ~ ns(intSCT2_5, 3) | IDAA,
  control = lmeControl(opt = 'optim'),
  data = dat
)

# DLI as competing event --------------------------------------------------


event_names <- levels(JM_dat_wide$endpoint5_s)
tmat <- trans.comprisk(
  K = 4,
  names = c("event_free", event_names[-1])
)

tmat

# Make numeric indicators
ind_cols <- paste0("ind_", event_names[-1])
JM_dat_wide[, (ind_cols) := lapply(
  event_names[-1], function(col) as.numeric(endpoint5_s == col)
)]

covs <- c("SCTyear_2010", "hirisk", "ATG", "VCMVPAT_pre")

JM_msdat <- msprep(
  time = c(NA, rep("endpoint5", times = length(ind_cols))),
  status = c(NA, ind_cols),
  data = data.frame(JM_dat_wide),
  trans = tmat,
  keep = covs,
  id = "IDAA"
)

JM_msdat_expand <- expand.covs(
  JM_msdat,
  covs,
  append = TRUE,
  longnames = FALSE
)

# Cannot do left truncation
coxCRfit <- coxph(
  Surv(Tstop, status) ~
    SCTyear_2010.1 + hirisk.1 + # cellular intervention
    ATG.2 + hirisk.2 + # relapse
    ATG.3 + # NRF other; then no covars for gvhd (only current val in JM)
    strata(trans) + cluster(IDAA),
  data = JM_msdat_expand,
  x = TRUE,
  model = TRUE
)

JMfit <- JM::jointModel(
  lmeObject = lmeFit,
  survObject = coxCRfit,
  timeVar = "intSCT2_5",
  method = "spline-PH-aGH",# or piecewise
  CompRisk = TRUE,
  interFact = list("value" = ~ strata(trans) - 1),
  iter.EM = 200
)


# DLI as censoring --------------------------------------------------------



event_names <- event_names[!(event_names %in% c("cell_interv"))]
tmat <- trans.comprisk(
  K = 3,
  names = c("event_free", event_names[-1])
)

tmat

# Make numeric indicators
JM_dat_wide[, endpoint5_s_dliCens := endpoint5_s]
JM_dat_wide[endpoint5_s == "cell_interv", endpoint5_s_dliCens := "cens"]
JM_dat_wide[, endpoint5_s_dliCens := droplevels(endpoint5_s_dliCens)]


ind_cols <- paste0("ind_dli_", event_names[-1])
JM_dat_wide[, (ind_cols) := lapply(
  event_names[-1], function(col) as.numeric(endpoint5_s_dliCens == col)
)]

covs <- c("hirisk", "ATG")

JM_msdat_dli <- msprep(
  time = c(NA, rep("endpoint5", times = length(ind_cols))),
  status = c(NA, ind_cols),
  data = data.frame(JM_dat_wide),
  trans = tmat,
  keep = covs,
  id = "IDAA"
)

JM_msdat_expand_dli <- expand.covs(
  JM_msdat_dli,
  covs,
  append = TRUE,
  longnames = FALSE
)

# Cannot do left truncation
coxCRfit_dli <- coxph(
  Surv(Tstop, status) ~
    ATG.1 + hirisk.1 +
    ATG.2 +
    strata(trans) + cluster(IDAA),
  data = JM_msdat_expand_dli,
  x = TRUE,
  model = TRUE
)

JMfit_dli <- JM::jointModel(
  lmeObject = lmeFit,
  survObject = coxCRfit_dli,
  timeVar = "intSCT2_5",
  method = "spline-PH-aGH",# or piecewise
  CompRisk = TRUE,
  interFact = list("value" = ~ strata(trans) - 1),
  iter.EM = 200
)

# Compare models ----------------------------------------------------------


round(summary(JMfit)$`CoefTable-Long`, 3)
round(summary(JMfit_dli)$`CoefTable-Long`, 3)
round(summary(JMfit)$`CoefTable-Event`, 3)[seq_len(9),]
round(summary(JMfit_dli)$`CoefTable-Event`, 3)[seq_len(9),]
