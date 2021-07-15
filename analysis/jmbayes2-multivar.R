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
  "CD4_abs_log" = ~ value(CD4_abs_log) * DLI * (strata(trans) - 1),
  "CD8_abs_log" = ~ value(CD8_abs_log) * DLI * (strata(trans) - 1)
)

# Add to targets pipeline
jFit_CR <- jm(
  Surv_object = coxCRfit,
  Mixed_objects = list(fm1, fm2),
  time_var = "intSCT2_5",
  functional_forms = CR_forms,
  data_Surv = coxCRfit$model, # try with this
  n_iter = 6000L,
  n_burnin = 2000L,
  priors = list("penalty_alphas" = "ridge", penalty_gammas = "ridge")
)

# Maybe also uncorrelated random effects

coef(jFit_CR)
CD4_all_dli$coefficients$gammas
CD4_all_dli$coefficients$alpha
CD8_all_dli$coefficients$alpha


# Pretty horrid non-covergence
ggtraceplot(jFit_CR, "gammas", grid = TRUE, gridcols = 3)
ggtraceplot(jFit_CR, "alphas", grid = TRUE, gridcols = 4)



# Attempt wit data_surv ---------------------------------------------------


coxCRfit_bis <- survival::coxph(
  Surv(Tstart, Tstop, status) ~
    hirisk.1 + ATG.1 + # Relapse submodel
    DLI.2 + ATG.2 + # GVHD submodel
    ATG.3 + # NRF_other submodel
    strata(trans),
  data = dli_msdata,
  x = TRUE,
  model = TRUE,
  id = IDAA
)

coxCRfit_bis <- survival::coxph(
  Surv(Tstart, Tstop, status) ~
    hirisk.1 + ATG.1 + DLI.1 + # Relapse submodel
    DLI.2 + ATG.2 + # GVHD submodel
    ATG.3 + DLI.3 + # NRF_other submodel
    strata(trans),
  data = dli_msdata,
  x = TRUE,
  model = TRUE,
  id = IDAA
)

dli_msdata$trans1 <- as.numeric(dli_msdata$trans == 1)
dli_msdata$trans2 <- as.numeric(dli_msdata$trans == 2)
dli_msdata$trans3 <- as.numeric(dli_msdata$trans == 3)

# Compare this one to JM estimate
CR_forms_bis <- list(
  "CD4_abs_log" = ~ value(CD4_abs_log):(trans1 + trans2 + trans3) +
    value(CD4_abs_log):trans2:DLI.2 - 1,
  "CD8_abs_log" = ~ value(CD8_abs_log):(trans1 + trans2 + trans3) +
    value(CD8_abs_log):trans2:DLI.2 - 1
)

CR_forms_av <- list(
  "CD4_abs_log" = ~ value(CD4_abs_log):(trans1 + trans2 + trans3) - 1,
  "CD8_abs_log" = ~ value(CD8_abs_log):(trans1 + trans2 + trans3) - 1
)

CR_forms <- list(
  "CD4_abs_log" = ~ value(CD4_abs_log):(trans1 + trans2 + trans2:DLI.2 + trans3) - 1,
  "CD8_abs_log" = ~ value(CD8_abs_log):(trans1 + trans2 + trans2:DLI.2 + trans3) - 1
)

CR_forms <- list(
  "CD4_abs_log" = ~ value(CD4_abs_log):(
    trans1 + trans1:DLI.1 + trans2 + trans2:DLI.2 + trans3 + trans3:DLI.3
  ) - 1,
  "CD8_abs_log" = ~ value(CD8_abs_log):(
    trans1 + trans1:DLI.1 + trans2 + trans2:DLI.2 + trans3 + trans3:DLI.3
  ) - 1
)

# Trans1 value, trans3 value and trans2 pre and post DLI
assoc_priors <- lapply(CR_forms, function(fform) {
  cr_form_terms <- terms(fform)
  n_alphas <- length(attr(cr_form_terms, "term.labels"))
  diag(5, n_alphas) # precision so sd is 1/5
})

# Or list of other length
n_assoc_terms <- sum(
  vapply(CR_forms, function(fform) {
    length(attr(terms(fform), "term.labels"))
  }, FUN.VALUE = integer(1L))
)

# https://github.com/drizopoulos/JMbayes2/blob/master/R/jm.R
# Add to targets pipeline
jFit_CR <- jm(
  Surv_object = coxCRfit_bis,
  Mixed_objects = list(fm1, fm2),
  time_var = "intSCT2_5",
  functional_forms = CR_forms,
  data_Surv = dli_msdata, # try with this
  n_iter = 12500L,
  n_burnin = 2500L,
  priors = list(
    # Local ridge priors for each alpha ~ N(0, 1/2)
    Tau_alphas = lapply(seq_len(n_assoc_terms), function(alpha) matrix(data = 2))
  )
)

# try 4 chains, longer burn-in
ggdensityplot(jFit_CR, "alphas", grid = TRUE, gridcols = 4)

coef(jFit_CR)
CD4_all_dli$coefficients$gammas

ggtraceplot(jFit_CR, "betas", grid = TRUE, gridcols = 5, gridrows = 5)
ggtraceplot(jFit_CR, "gammas", grid = TRUE, gridcols = 5)
ggtraceplot(jFit_CR, "alphas", grid = TRUE, gridcols = 4, gridrows = 3)


jFit_CR_bis <- jm(
  Surv_object = coxCRfit_bis,
  Mixed_objects = list(fm1, fm2),
  time_var = "intSCT2_5",
  functional_forms = CR_forms,
  data_Surv = dli_msdata
)

# Maybe much more informative priors?

jFit_CR_bis
coef(jFit_CR_bis)



CD4_all_dli_avfform$coefficients$gammas
CD4_all_dli_avfform$coefficients$alpha
CD8_all_dli_avfform$coefficients$alpha

CD4_all_dli$coefficients$gammas
CD4_all_dli$coefficients$alpha
CD8_all_dli$coefficients$alpha


# Pretty horrid non-covergence
ggtraceplot(jFit_CR_bis, "betas", grid = TRUE, gridcols = 4)
ggtraceplot(jFit_CR_bis, "gammas", grid = TRUE, gridcols = 4)
ggtraceplot(jFit_CR_bis, "alphas", grid = TRUE, gridcols = 4)

predict(jFit_CR_bis, newdata = datasets$long[IDAA == "221"],
        times = seq(7, 12, length.out = 51))

newdat <- expand.grid(
  "ATG" = levels(dat_wide$ATG),
  "VCMVPAT_pre" = levels(dat_wide$VCMVPAT_pre),
  "intSCT2_5" = seq(0.1, 24, by = 0.1)
)

modmat_newdat <- model.matrix(~ ns(intSCT2_5, 3) * ATG + VCMVPAT_pre, data = newdat)
betas <- drop(summary(jFit_CR_bis)$Outcome1$Mean[-10])
coefficients(jFit_CR_bis)
newdat$preds <- modmat_newdat %*% betas

newdat %>%
  ggplot(aes(intSCT2_5, preds)) +
  geom_line(aes(col = ATG)) +
  facet_wrap(~ VCMVPAT_pre)
