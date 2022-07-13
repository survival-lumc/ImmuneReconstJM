# ONLY IF TIME (for fun):
# Try dynamic prediction with new jmbayes2 (also to double check results)
# - need to use cr_setup() (interaction parametrisation)
# - Try bivariate CD4 and CD8, or just CD3

tar_load(c(dat_merged, NMA_preDLI_datasets))

# This was syntax with mstate parametrisation (just for guidance)
jm_mod <- jm(
  Surv_object = preDLI_cox,
  Mixed_objects = list(preDLI_CD3_long_reffcorr),
  time_var = "intSCT2_5",
  #"CD3_abs_log" = ~ value(CD3_abs_log) + slope(CD3_abs_log) +
  #  (value(CD3_abs_log) + slope(CD3_abs_log)):(strata(trans) - 1)
  functional_forms = list(
    "CD3_abs_log" = ~ value(CD3_abs_log) + value(CD3_abs_log):(strata(trans) - 1)
  ),
  id_var = "IDAA",
  n_iter = 13000L,
  n_burnin = 3000L
  #priors = list("penalty_alphas" = "horseshoe") #"horseshoe")
  # Numver of alphas
  #priors = list(
  #  Tau_alphas = lapply(seq_len(3), function(alpha) matrix(data = 0.75))
  #),
  #base_hazard_segments = 5
)

ggtraceplot(jm_mod, "alphas", gridrows = 1, gridcols = 3, grid = TRUE)
ggtraceplot(jm_mod, "betas", gridrows = 6, gridcols = 3, grid = TRUE)

gelman_diag(jm_slopes)



# Try dynamic pred --------------------------------------------------------

dat_long <- NMA_preDLI_datasets$long
dat_wide <- NMA_preDLI_datasets$wide

# crisk_setup(...)
# ...




# New set of test! --------------------------------------------------------

library(targets)
library(survival)
library(JM)
library(JMbayes2)

tar_load(
  c(
    NMA_preDLI_datasets,
    preDLI_JM_value_corr_CD4,
    preDLI_long_corr_CD4,
    preDLI_cox
  )
)

# Survival set-up
dat_long <- NMA_preDLI_datasets$long
dat_wide <- NMA_preDLI_datasets$wide

dat_wide_idCR <- crisk_setup(
  dat_wide,
  statusVar = "endpoint7_s",
  censLevel = "cens",
  nameStrata = "CR"
)

# For transition-specific covars:
dat_wide_idCR$ATG.1 <- with(dat_wide_idCR, (as.numeric(ATG) - 1L) * (CR == "gvhd"))
dat_wide_idCR$hirisk.1 <- with(dat_wide_idCR, (as.numeric(hirisk) - 1L) * (CR == "gvhd"))
dat_wide_idCR$ATG.2 <- with(dat_wide_idCR, (as.numeric(ATG) - 1L) * (CR == "relapse"))
dat_wide_idCR$hirisk.2 <- with(dat_wide_idCR, (as.numeric(hirisk) - 1L) * (CR == "relapse"))
dat_wide_idCR$ATG.3 <- with(dat_wide_idCR, (as.numeric(ATG) - 1L) * (CR == "other_nrf"))



# First just GVHD ---------------------------------------------------------


coxfit_gvhd <- coxph(Surv(endpoint7, endpoint7_s == "gvhd") ~ ATG + hirisk, data = dat_wide)
jm_gvhd <- jm(
  Surv_object = coxfit_gvhd,
  Mixed_objects = preDLI_long_corr_CD4,
  time_var = "intSCT2_7",
  functional_forms = list(CD4_abs_log = ~ value(CD4_abs_log))
)

# Some diagnostics
jm_gvhd

# Comparing parameters
summary(preDLI_JM_value_corr_CD4)$`CoefTable-Event`[1:8, ]
preDLI_long_corr_CD4 |> VarCorr()
?GLMMadaptive::mixed_model

adaptive_mixed <- mixed_model(
  fixed = CD4_abs_log ~ ns(intSCT2_7, 3) * hirisk * ATG + CMV_PD,
  random = ~ 0 + ns(intSCT2_7, 3) | IDAA,
  data = dat_long,
  family = gaussian()
)

library(brms)
library(splines)
library(rstan)
mmod_bayes <- brm(
  #CD4_abs_log ~ ns(intSCT2_7, 3) * hirisk * ATG + CMV_PD + (0 + ns(intSCT2_7, 3) | IDAA),
  formula = CD4_abs_log ~ intSCT2_7,
  data  = dat_long
)
summary(mmod_bayes)

?rstan
mtcars$mpg10 <- mtcars$mpg / 10
library(rstanarm)
fit <- stan_glm(
  mpg10 ~ wt + cyl + am,
  data = mtcars,
  QR = TRUE,
  # for speed of example only (default is "sampling")
  algorithm = "sampling"
)
plot(fit, plotfun = "trace")

fit_basic <- stan_lmer(
     formula = CD4_abs_log ~ ns(intSCT2_7, 3) + (0 + ns(intSCT2_7, 3) | IDAA),
     data = dat_long
  )

options(mc.cores = 3)
fit <- stan_lmer(
  formula = CD4_abs_log ~ ns(intSCT2_7, 3) * ATG + CMV_PD + (0 + ns(intSCT2_7, 3) | IDAA),
  data = dat_long,
  chains = 3,
  cores = 3
)


fit2 <- stan_lmer(
  formula = CD4_abs_log ~ ns(intSCT2_7, 3) * ATG * hirisk + CMV_PD + (0 + ns(intSCT2_7, 3) | IDAA),
  data = dat_long,
  chains = 3,
  cores = 3
)

plot(fit, plotfun = "trace")
fit$fitted.values
plot(dat_long$intSCT2_7, fit$fitted.values)

library("tidybayes")
library(modelr)
library(ggplot2)
library(splines)

dat_long |>
  data_grid(intSCT2_7 = seq(0, 6, by = 0.1)) |>
  add_epred_draws(fit_basic, ndraws = 100, re_formula = NA
  ) |>
  ggplot(aes(intSCT2_7, .epred, group = .draw)) +
  geom_line()

dat_long |>
  data_grid(intSCT2_7 = seq(0, 6, by = 0.1), IDAA = 221) |>
  add_epred_draws(fit, ndraws = 100
          #,re_formula =
                  ) |>
  ggplot(aes(intSCT2_7, .epred, group = .draw)) +
  geom_line()

# Now with complex one


dat_long |>
  data_grid(intSCT2_7 = seq(0, 6, by = 0.1), CMV_PD, ATG) |>
  add_epred_draws(fit, ndraws = 100, re_formula = NA) |>
  ggplot(aes(intSCT2_7, .epred)) +
  geom_line(
    aes(
      group = paste(.draw, CMV_PD, ATG),
      col = paste(ATG)
    )
  ) +
  facet_wrap(~ CMV_PD, nrow = 2)



dat_long |>
  data_grid(intSCT2_7 = seq(0, 6, by = 0.1), CMV_PD, ATG, hirisk) |>
  add_epred_draws(fit2, ndraws = 100, re_formula = NA) |>
  ggplot(aes(intSCT2_7, .epred)) +
  geom_line(
    aes(
      group = paste(.draw, CMV_PD, ATG, hirisk),
      col = paste(ATG, hirisk)#@,
    ),
    alpha = 0.5
  ) +
  facet_wrap(~ CMV_PD, nrow = 2)


dat_long |>
  data_grid(intSCT2_7 = seq(0, 6, by = 0.1), CMV_PD = "-/-", ATG = "UD(+ATG)", hirisk) |>
  add_epred_draws(fit2, ndraws = 100, re_formula = NA) |>
  ggplot(aes(intSCT2_7, .epred)) +
  geom_line(
    aes(
      group = paste(.draw, hirisk),
      col = paste(hirisk)#@,
    ),
    alpha = 0.25
  ) +
  theme_minimal(
  )


dat_long |>
  data_grid(intSCT2_7 = seq(0, 6, by = 0.25), CMV_PD = "-/-", ATG = "UD(+ATG)", hirisk) |>
  add_predicted_draws(fit2, re_formula = NA) |>
  #add_epred_draws(fit2, ndraws = 100, re_formula = NA) |>
  ggplot(aes(intSCT2_7, group = hirisk, col = hirisk)) +
  stat_lineribbon(aes(y = .prediction, col = hirisk), .width = c(.95),
                  color = "#08519C", alpha = 0.5) +
  scale_fill_brewer()




# Try a stanjm

# Only less than 20 mins for this beast!!
stan_gvhd <- stan_jm(
  formulaLong = CD4_abs_log ~ ns(intSCT2_7, 3) * ATG * hirisk + CMV_PD
  + (0 + ns(intSCT2_7, 3) | IDAA),
  dataLong = dat_long,
  formulaEvent = survival::Surv(endpoint7, endpoint7_s == "gvhd") ~ ATG + hirisk,
  dataEvent = dat_wide,
  time_var = "intSCT2_7",
  chains = 3,
  refresh = 20,
  seed = 12345,
  cores = 3
)

library(JM)
summary(preDLI_JM_value_corr_CD4)
summary(stan_gvhd)
stan_gvhd
VarCorr(stan_gvhd)
launch_shinystan(stan_gvhd)

ids <- sample(dat_long$IDAA, size = 16)
plot(posterior_traj(stan_gvhd, extrapolate = TRUE, ids = ids),
     plot_observed = TRUE)

# Pop level
dat_long$CMV_PD
ndL <- expand.grid(intSCT2_7 = seq(0.01, 6, by = 0.1),
                   CMV_PD = factor(levels(dat_long$CMV_PD),
                                   levels = levels(dat_long$CMV_PD)),
                   ATG = factor(levels(dat_long$ATG),
                                levels = levels(dat_long$ATG)),
                   hirisk = factor(levels(dat_long$hirisk),
                                   levels = levels(dat_long$hirisk)))
ndL$IDAA <- as.character(with(ndL, interaction(CMV_PD, ATG, hirisk))) #sample(dat_long$IDAA, size = nrow(ndL))
a <- posterior_traj(stan_gvhd, newdataLong = ndL, dynamic = FALSE)
plot(a)

a <- posterior_traj(stan_gvhd, newdataLong = cbind.data.frame(
  intSCT2_7 = seq(0.01, 6, by = 0.1),
  CMV_PD = "-/-",
  ATG = "UD(+ATG)",
  hirisk = "yes",
  IDAA = 1
), dynamic = F)
plot(a)

posterior_predict(
  stan_gvhd,
  newdata = cbind.data.frame(
    intSCT2_7 = seq(0.01, 6, by = 0.1),
    CMV_PD = "-/-",
    ATG = "UD(+ATG)",
    hirisk = "yes"
  ),
  re.form = NA
) |>
  dim()

p1 <- posterior_traj(stan_gvhd,
              newdataLong = ndL,
              dynamic = FALSE,   re.form = NA)

test <- posterior_traj(stan_gvhd,
               newdataLong = ndL,
               dynamic = FALSE)

posterior_traj(stan_gvhd,
               newdataLong = ndL,
               dynamic = FALSE)
test |>  ggplot(aes(intSCT2_7, yfit)) +
  geom_line(aes(group = IDAA, col = IDAA))
plot(test)


# Try also with slope!!!
options()


# Checks ------------------------------------------------------------------


# Check jmbayes2 example of univariat vs stan_jm and JM


pbc2.id$status2 <- as.numeric(pbc2.id$status != 'alive')
CoxFit <- coxph(Surv(years, status2) ~ sex, data = pbc2.id, x = TRUE)
fm1 <- lme(log(serBilir) ~ year * sex, data = pbc2, random = ~ year + 0 | id)
univar_JMbayes2 <- jm(CoxFit, fm1, time_var = "year")
#coef(univar_JMbayes2, "Event")

univar_JM <- jointModel(
  lmeObject = fm1,
  survObject = CoxFit,
  method = "spline-PH-aGH",
  timeVar = "year",
  parameterization = "value",
  numeriDeriv = "cd",
  lng.in.kn = 6
)
coef(univar_JMbayes2, "Event")
coef(univar_JM, "Event")


# Cause-specific GVHD example JM vs stan_jm vs JMbayes vs Jmbayes2 --------



library(targets)
library(rstanarm)
library(JM)
library(JMbayes)
library(JMbayes2)

options(mc.cores = 3)

# Load objects
tar_load(
  c(
    NMA_preDLI_datasets,
    preDLI_JM_value_corr_CD4,
    preDLI_long_corr_CD4,
    preDLI_cox
  )
)

# Separate datasets
dat_long <- NMA_preDLI_datasets$long
dat_wide <- NMA_preDLI_datasets$wide

# Sub models:
coxfit <- coxph(
  Surv(endpoint7, endpoint7_s == "gvhd") ~ ATG + hirisk,
  data = dat_wide,
  x = TRUE
)
# preDLI_long_corr_CD4 is longitudinal

## Method 1: JM
gvhd_JM <- jointModel(
  lmeObject = preDLI_long_corr_CD4,
  survObject = coxfit,
  method = "spline-PH-aGH",
  timeVar = "intSCT2_7",
  parameterization = "value",
  numeriDeriv = "cd",
  lng.in.kn = 6
)

## Method 2: JMbayes
gvhd_JMbayes <- jointModelBayes(
  lmeObject = preDLI_long_corr_CD4,
  survObject = coxfit,
  baseHaz = "regression-splines",
  lng.in.kn = 6,
  timeVar = "intSCT2_7",
  n.thin = 1,
  n.burnin = 10000 # make this looonger
)

summary(gvhd_JM)
summary(gvhd_JMbayes)
plot(gvhd_JMbayes)

## Method 3: JMbayes2 (see https://github.com/drizopoulos/JMbayes2/issues/24)
# .. not a not running long enough issues.. maybe splines problem?
gvhd_JMbayes2 <- jm(
  Surv_object = coxfit,
  Mixed_objects = preDLI_long_corr_CD4,
  time_var = "intSCT2_7",
  Bsplines_degree = 3L,
  base_hazard_segments = 7L,
  n_burnin = 30000,
  n_iter = 50000,
  n_thin = 2
)

ggtraceplot(gvhd_JMbayes2, c("betas"), grid = TRUE, gridrows = 6, gridcols = 3)
ggtraceplot(gvhd_JMbayes2, c("gammas"), grid = TRUE, gridrows = 2)
ggtraceplot(gvhd_JMbayes2, c("alphas"))

coef(gvhd_JMbayes2)
  coef(gvhd_JMbayes2, process = "Event")

## Method 4: stan_jm
gvhd_stan <- stan_jm(
  formulaLong = CD4_abs_log ~ ns(intSCT2_7, 3) * ATG * hirisk + CMV_PD
  + (0 + ns(intSCT2_7, 3) | IDAA),
  dataLong = dat_long,
  formulaEvent = survival::Surv(endpoint7, endpoint7_s == "gvhd") ~ ATG + hirisk,
  dataEvent = dat_wide,
  time_var = "intSCT2_7",
  chains = 3,
  refresh = 20,
  seed = 12345,
  cores = 3
)


# Try with random intercepts too!
gvhd_stan_randint <- stan_jm(
  formulaLong = CD4_abs_log ~ ns(intSCT2_7, 3) * ATG * hirisk + CMV_PD
  + (ns(intSCT2_7, 3) | IDAA),
  dataLong = dat_long,
  formulaEvent = survival::Surv(endpoint7, endpoint7_s == "gvhd") ~ ATG + hirisk,
  dataEvent = dat_wide,
  time_var = "intSCT2_7",
  chains = 3,
  refresh = 20,
  seed = 12345,
  cores = 3
)

gvhd_stan_randint
VarCorr(gvhd_stan_randint)
summary(gvhd_stan_randint)


traceplot(gvhd_stan_randint)

#
ids <- sample(dat_long$IDAA, size = 16)
plot(posterior_traj(gvhd_stan_randint, extrapolate = TRUE, ids = ids),
     plot_observed = TRUE)

test <- posterior_traj(gvhd_stan_randint, dynamic = FALSE, re.form = NA, newdataLong = ndL)
plot(test)
test |>
  ggplot(aes(intSCT2_7, yfit, group = IDAA)) +
  geom_line(aes(col = IDAA))

pp_check(gvhd_stan_randint, m = 1)
plot(gvhd_stan_randint, plotfun = "mcmc_trace", regex_pars = "^Event")
plot(gvhd_stan_randint, regex_pars = "^Long")
plot(gvhd_stan_randint, plotfun = "mcmc_trace", regex_pars = "^Long")
plot(gvhd_stan_randint, plotfun = "mcmc_hist", pars = "Assoc|Long1|etavalue")

bayesplot::available_mcmc()
gvhd_stan_randint$coefficients

# One more time but with rand int
# it fookin works!! Try competing risks one now?
gvhd_JMbayes2_randint <- jm(
  Surv_object = coxfit,
  Mixed_objects = update(preDLI_long_corr_CD4, random = ~ ns(intSCT2_7, 3) | IDAA),
  time_var = "intSCT2_7",
  n_burnin = 10000,
  n_iter = 20000
)

gvhd_JMbayes2_randint
ggtraceplot(gvhd_JMbayes2_randint, c("betas"), grid = TRUE, gridrows = 6, gridcols = 3)
ggtraceplot(gvhd_JMbayes2_randint, c("gammas"), grid = TRUE, gridrows = 2)
ggtraceplot(gvhd_JMbayes2_randint, c("alphas"))



# Comprirsk JMbayes2 with rand int ----------------------------------------


# Survival set-up
dat_long <- NMA_preDLI_datasets$long
dat_wide <- NMA_preDLI_datasets$wide

dat_wide_idCR <- crisk_setup(
  dat_wide,
  statusVar = "endpoint7_s",
  censLevel = "cens",
  nameStrata = "CR"
)

# For transition-specific covars:
dat_wide_idCR$ATG.1 <- with(dat_wide_idCR, (as.numeric(ATG) - 1L) * (CR == "gvhd"))
dat_wide_idCR$hirisk.1 <- with(dat_wide_idCR, (as.numeric(hirisk) - 1L) * (CR == "gvhd"))
dat_wide_idCR$ATG.2 <- with(dat_wide_idCR, (as.numeric(ATG) - 1L) * (CR == "relapse"))
dat_wide_idCR$hirisk.2 <- with(dat_wide_idCR, (as.numeric(hirisk) - 1L) * (CR == "relapse"))
dat_wide_idCR$ATG.3 <- with(dat_wide_idCR, (as.numeric(ATG) - 1L) * (CR == "other_nrf"))

CoxFit <-  coxph(
  Surv(endpoint7, status2) ~
    ATG.1 + hirisk.1 +
    ATG.2 + hirisk.2 +
    ATG.3 + strata(CR),
  data = dat_wide_idCR,
  x = TRUE
)
lmeFit <- update(preDLI_long_corr_CD4, random = ~ ns(intSCT2_7, 3) | IDAA)

gvhd_JMbayes2_randint_CR <- jm(
  Surv_object = CoxFit,
  Mixed_objects = lmeFit,
  time_var = "intSCT2_7",
  functional_forms = list("CD4_abs_log" = ~ value(CD4_abs_log):CR),
  n_burnin = 10000, # needs 10000 more..
  n_iter = 20000
)

gvhd_JMbayes2_randint_CR
ggtraceplot(gvhd_JMbayes2_randint_CR, c("betas"), grid = TRUE, gridrows = 6, gridcols = 3)
ggtraceplot(gvhd_JMbayes2_randint_CR, c("gammas"), grid = TRUE, gridrows = 2)
ggtraceplot(gvhd_JMbayes2_randint_CR, c("alphas"))


# Even more basic example -------------------------------------------------


library(JM)
library(JMbayes)
library(JMbayes2)
library(rstanarm)
options(mc.cores = 3)

head(aids)
lmeFit <- lme(log(serBilir) ~ ns(year, 2) * drug, random = ~ 0 + ns(year, 2) | id, data = pbc2)
coxFit <- coxph(Surv(years, status2) ~ drug, data = pbc2.id, x = TRUE)


# JM
aids_JM <- jointModel(
  lmeObject = lmeFit,
  survObject = coxFit,
  method = "spline-PH-aGH",
  timeVar = "year"#3,
  #lng.in.kn = 6
)

# JMbayes
aids_JMbayes <- jointModelBayes(
  lmeObject = lmeFit,
  survObject = coxFit,
  timeVar = "year",
  baseHaz = "regression-splines"#,
  #lng.in.kn = 6
)

# JMbayes2
aids_JMbayes2 <- jm(
  Mixed_objects = lmeFit,
  Surv_object = coxFit,
  time_var = "year"#,
  #Bsplines_degree = 3L,
  #base_hazard_segments = 7L,
)

aids_stan <- stan_jm(
  formulaLong = log(serBilir) ~ ns(year, 2) * drug + (0 + ns(year, 2) | id),
  dataLong = pbc2,
  formulaEvent = survival::Surv(years, status2) ~ drug,
  dataEvent = pbc2.id,
  time_var = "year",
  chains = 3,
  refresh = 50,
  cores = 3
)

# https://github.com/drizopoulos/JMbayes2/blob/f4998f398b863969c125199c090c32c1fd7c70ee/R/jm.R#L130

#coef(aids_JM, "Event")
#coef(aids_JMbayes, "Event")
#coef(aids_JMbayes2, "Event")
cbind(
  "lme" = fixef(lmeFit),
  "JM" = fixef(aids_JM),
  "stan" = fixef(aids_stan)$Long1,
  "JMbayes" = fixef(aids_JMbayes),
  "JMbayes2" = fixef(aids_JMbayes2)$`log(serBilir)`
)


