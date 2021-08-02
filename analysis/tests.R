kts <- c(3, 6, 8)
modo <- lme(
  fixed = CD3_abs_log ~ ns(intSCT2_5, df = 4) * ATG + VCMVPAT_pre,
  #random = ~ ns(intSCT2_5, 4) | IDAA,
  random = list(IDAA = pdDiag(~ ns(intSCT2_5, df = 4))),
  control = lmeControl(opt = "optim", msMaxIter = 100),
  data = data.frame(dat)
)

df <- cbind.data.frame(
  newdat,
  "preds" = predict(modo, newdata = newdat, level = 0L)
)

ggplot(df, aes(intSCT2_5, preds, col = ATG)) +
  geom_line() +
  facet_wrap(~ VCMVPAT_pre)


# Re-attempt on changepoint DLI
dat <- datasets$long

# Try again with time since interm


dat[, uDLI_t := ifelse(uDLI >= endpoint6, 1e5, uDLI), by = "IDAA"]
dat[, uDLI_ind := as.numeric(intSCT2_5 >= uDLI_t)]
dat[, uDLI_tsince := pmax(0, intSCT2_5 - uDLI_t)]


# Try mixed mod
lmeFit_tdep <- lme(
  fixed = CD4_abs_log ~ intSCT2_5 * ATG + VCMVPAT_pre + uDLI_tsince,
  random = ~ intSCT2_5 + uDLI_tsince | IDAA,
  control = lmeControl(opt = "optim", msMaxIter = 100),
  data = data.frame(dat)
)
lmeFit_tdep

t_interm <- 8
max_t <- 24
newdat <- expand.grid(
  "ATG" = levels(dat_wide$ATG),
  "VCMVPAT_pre" = levels(dat_wide$VCMVPAT_pre),
  "intSCT2_5" = seq(0.1, max_t, by = 0.1)
)
newdat$uDLI_t <- ifelse(t_interm >= max_t, 1e5, t_interm)
newdat$uDLI_ind <- as.numeric(newdat$intSCT2_5 >= newdat$uDLI_t)
newdat$uDLI_tsince <- pmax(0, newdat$intSCT2_5 - newdat$uDLI_t)
newdat


df <- cbind.data.frame(
  newdat,
  "preds" = predict(lmeFit_tdep, newdata = newdat, level = 0L)
)

ggplot(df, aes(intSCT2_5, preds, col = ATG)) +
  geom_line() +
  facet_wrap(~ VCMVPAT_pre)



# Try JMbayes 2 with more flexible spline and pdDiag ----------------------


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

fm1 <- do.call(
  nlme::lme,
  list(
    "fixed" = CD8_abs_log ~ ns(intSCT2_5, df = 4) * ATG + VCMVPAT_pre,
    "random" = list(IDAA = pdDiag(form = ~ ns(intSCT2_5, df = 4))),
    #random" = ~ ns(intSCT2_5, df = 4) | IDAA,
    "control" = nlme::lmeControl(opt = "optim"),
    "data" = datasets$long
  )
)

fm2 <- do.call(
  nlme::lme,
  list(
    "fixed" = CD4_abs_log ~ ns(intSCT2_5, df = 4) * ATG + VCMVPAT_pre,
    "random" = list(IDAA = pdDiag(form = ~ ns(intSCT2_5, df = 4))),
    #"random" = ~ ns(intSCT2_5, 3) | IDAA,
    "control" = nlme::lmeControl(opt = "optim"),
    "data" = datasets$long
  )
)

CR_forms <- list(
  "CD8_abs_log" = ~ value(CD8_abs_log):(
    trans1 + trans1:DLI.1 + trans2 + trans2:DLI.2 + trans3 + trans3:DLI.3
  ) - 1,
  "CD4_abs_log" = ~ value(CD4_abs_log):(
    trans1 + trans1:DLI.1 + trans2 + trans2:DLI.2 + trans3 + trans3:DLI.3
  ) - 1
)

n_assoc_terms <- sum(
  vapply(CR_forms, function(fform) {
    length(attr(terms(fform), "term.labels"))
  }, FUN.VALUE = integer(1L))
)

dli_msdata$trans1 <- as.numeric(dli_msdata$trans == 1)
dli_msdata$trans2 <- as.numeric(dli_msdata$trans == 2)
dli_msdata$trans3 <- as.numeric(dli_msdata$trans == 3)

jFit_CR <- jm(
  Surv_object = coxCRfit,
  Mixed_objects = list(fm1, fm2), #give names?
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

coef(jFit_CR)
jFit_CR
fm1$coefficients$fixed
ggtraceplot(jFit_CR, "betas", grid = TRUE, gridcols = 4, gridrows = 3)
ggtraceplot(jFit_CR, "alphas", grid = TRUE, gridcols = 4, gridrows = 3)

# Try basic predicting
dat_wide <- datasets$wide
newdat <- expand.grid(
  "ATG" = levels(dat_wide$ATG),
  "VCMVPAT_pre" = levels(dat_wide$VCMVPAT_pre),
  "intSCT2_5" = seq(0.1, 24, by = 0.1)
)

mean_coefs <- summary(jFit_CR)$Outcome2[, "Mean"]
coefs <- mean_coefs[-length(mean_coefs)]
spliney <- ns(datasets$long$intSCT2_5, 4)
bound_knots <- attr(spliney, "Boundary.knots")
int_knots <- attr(spliney, "knots")
X <- model.matrix(
  ~ ns(intSCT2_5, Boundary.knots = bound_knots, knots = int_knots) * ATG + VCMVPAT_pre,
  data = newdat
)

# Make predict function..

df <- cbind.data.frame(
  newdat,
  "preds" = X %*% coefs
)
4
ggplot(df, aes(intSCT2_5, preds, col = ATG)) +
  geom_line() +
  facet_wrap(~ VCMVPAT_pre)

jFit_CR$model_info$terms$terms_FE_noResp[[1]] # this is formula on its own
lol <- attr(jFit_CR$model_info$terms$terms_FE_noResp[[1]], "predvars")
heyz <- do.call(cbind.data.frame, eval(lol, envir = newdat))
as.character(lol)

# or just re-eval spliine - MODEL.FRAME WILL KEEP SPLINE LOCATIONS
hoze <- model.frame.default(jFit_CR$model_info$terms$terms_FE_noResp[[1]], data = newdat)
model.matrix.default(jFit_CR$model_info$terms$terms_FE_noResp[[1]], hoze)
