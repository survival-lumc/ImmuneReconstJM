# Test of bivariate CD4 and CD8 model -------------------------------------


library(data.table)
library(JM)
library(JMbayes2)
#library(ggsankey)
library(tidyverse) # reduce to only ggplot2 later
library(mstate)

coxCRfit <- survival::coxph(
  Surv(Tstart, Tstop, status) ~
    hirisk.1 + DLI.1 + ATG.1 + # Relapse submodel
    DLI.2 + ATG.2 + # GVHD submodel
    DLI.3 + ATG.3 + # NRF_other submodel
    strata(trans),
  data = dli_msdata,
  x = TRUE,
  model = TRUE,
  id = IDAA
)

# Prep long models
fm1 <- do.call(
  nlme::lme,
  list(
    "fixed" = CD8_abs_log ~ ns(intSCT2_5, 3) * ATG + VCMVPAT_pre,
    #random = list(IDAA = pdDiag(form = ~ ns(intSCT2_2, 3))),
    "random" = ~ ns(intSCT2_5, 3) | IDAA,
    "control" = nlme::lmeControl(opt = "optim"),
    "data" = datasets$long
  )
)

fm2 <- do.call(
  nlme::lme,
  list(
    "fixed" = CD4_abs_log ~ ns(intSCT2_5, 3) * ATG + VCMVPAT_pre,
    #random = list(IDAA = pdDiag(form = ~ ns(intSCT2_2, 3))),
    "random" = ~ ns(intSCT2_5, 3) | IDAA,
    "control" = nlme::lmeControl(opt = "optim"),
    "data" = datasets$long
  )
)

# Run joint fit
CR_forms <- list(
  "CD4_abs_log" = ~ DLI * (value(CD4_abs_log) + value(CD4_abs_log):(strata(trans) - 1)),
  "CD8_abs_log" = ~ DLI * (value(CD8_abs_log) + value(CD8_abs_log):(strata(trans) - 1))
)

CR_forms <- list(
  "CD4_abs_log" = ~ value(CD4_abs_log) * (strata(trans) - 1),
  "CD8_abs_log" = ~ value(CD8_abs_log) * (strata(trans) - 1)
)

jFit_CR <- jm(
  Surv_object = coxCRfit,
  Mixed_objects = list(fm1, fm2),
  time_var = "intSCT2_5",
  functional_forms = CR_forms,
  n_iter = 6000L,
  n_burnin = 2000L,
  priors = list("penalty_alphas" = "ridge")
)

coef(jFit_CR)
CD4_all_dli_avfform$coefficients$gammas
CD4_all_dli_avfform$coefficients$alpha
CD8_all_dli_avfform$coefficients$alpha


# Pretty horrid non-covergence
ggtraceplot(jFit_CR, "alphas", grid = TRUE, gridcols = 3)
