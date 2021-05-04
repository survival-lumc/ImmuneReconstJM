# Test of bivariate JMbayes2
library(mstate)
library(JMbayes2)

tar_load(dat_merged)

datasets <- prepare_JM_data(
  merged_data = dat_merged,
  c("CD8_abs_log", "CD4_abs_log"),
  admin_cens_time = 7
)

fform <- ~ trans2 + trans3 + trans4 - 1

# - Prep mstate
JM_dat_wide <- data.table::copy(datasets$wide)
event_names <- levels(JM_dat_wide$endpoint5_s)
tmat <- mstate::trans.comprisk(
  K = 4,
  names = c("event_free", event_names[-1])
)

# Make numeric indicators
ind_cols <- paste0("ind_", event_names[-1])
JM_dat_wide[, (ind_cols) := lapply(
  event_names[-1], function(col) as.numeric(endpoint5_s == col)
)]

covs <- c("SCTyear_2010", "hirisk", "ATG", "VCMVPAT_pre")

JM_msdat <- mstate::msprep(
  time = c(NA, rep("endpoint5", times = length(ind_cols))),
  status = c(NA, ind_cols),
  data = data.frame(JM_dat_wide),
  trans = tmat,
  keep = covs,
  id = "IDAA"
)

JM_msdat_expand <- mstate::expand.covs(
  JM_msdat,
  covs,
  append = TRUE,
  longnames = FALSE
)

# Cannot do left truncation
coxCRfit <- survival::coxph(
  Surv(Tstart, Tstop, status) ~
    SCTyear_2010.1 + hirisk.1 + # cellular intervention
    ATG.2 + hirisk.2 + # relapse
    ATG.3 + # NRF other; then no covars for gvhd (only current val in JM)
    strata(trans) + cluster(IDAA),
  data = JM_msdat_expand,
  x = TRUE,
  model = TRUE
)

# Prep functional form
coxCRfit$model$trans2 <- as.numeric(coxCRfit$model$`strata(trans)` == "trans=2")
coxCRfit$model$trans3 <- as.numeric(coxCRfit$model$`strata(trans)` == "trans=3")
coxCRfit$model$trans4 <- as.numeric(coxCRfit$model$`strata(trans)` == "trans=4")

# Prep long models
fm1 <- do.call(
  nlme::lme,
  list(
    "fixed" = CD8_abs_log ~ splines::ns(intSCT2_5, 3) * ATG + VCMVPAT_pre,
    #random = list(IDAA = pdDiag(form = ~ ns(intSCT2_2, 3))),
    "random" = ~ splines::ns(intSCT2_5, 3) | IDAA,
    "control" = nlme::lmeControl(opt = "optim"),
    "data" = datasets$long
  )
)

fm2 <- do.call(
  nlme::lme,
  list(
    "fixed" = CD4_abs_log ~ splines::ns(intSCT2_5, 3) * ATG + VCMVPAT_pre,
    #random = list(IDAA = pdDiag(form = ~ ns(intSCT2_2, 3))),
    "random" = ~ splines::ns(intSCT2_5, 3) | IDAA,
    "control" = nlme::lmeControl(opt = "optim"),
    "data" = datasets$long
  )
)
trans2 + trans3 + trans4 - 1

# Run joint fit
CR_forms <- list(
  "CD4_abs_log" = ~ value(CD4_abs_log):(strata(trans) - 1),
  "CD8_abs_log" = ~ value(CD8_abs_log):(strata(trans) - 1)
)

jFit_CR <- jm(coxCRfit, list(fm1, fm2), time_var = "intSCT2_5",
              functional_forms = CR_forms,
              n_iter = 10000L, n_burnin = 5000L, n_thin = 2L)

mmat <- coxCRfit$model
model.matrix(~ ATG.2:(strata(trans) - 1), data = mmat)

# Pretty horrid non-covergence
ggtraceplot(jFit_CR, "betas")
