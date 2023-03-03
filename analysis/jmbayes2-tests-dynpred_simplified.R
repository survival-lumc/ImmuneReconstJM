source(here::here("packages.R"))
library(JMbayes2)

tar_load(
  c(
    NMA_preDLI_datasets,
    preDLI_long_corr_CD4,
    preDLI_long_corr_CD8,
    preDLI_cox
  )
)

dat_wide_pre <- NMA_preDLI_datasets$wide
dat_long_pre <- NMA_preDLI_datasets$long

dat_cr <- crisk_setup(
  dat_wide_pre,
  statusVar = "endpoint7_s",
  censLevel = "cens",
  nameStrata = "CR"
)

CoxFit_CR <- coxph(
  Surv(endpoint7, status2) ~ (ATG + hirisk) * strata(CR),
  data = dat_cr
)

# For transition-specific covars:
dat_cr$ATG.1 <- with(dat_cr, (as.numeric(ATG) - 1L) * (CR == "gvhd"))
dat_cr$hirisk.1 <- with(dat_cr, (as.numeric(hirisk) - 1L) * (CR == "gvhd"))
dat_cr$ATG.2 <- with(dat_cr, (as.numeric(ATG) - 1L) * (CR == "relapse"))
dat_cr$hirisk.2 <- with(dat_cr, (as.numeric(hirisk) - 1L) * (CR == "relapse"))
dat_cr$ATG.3 <- with(dat_cr, (as.numeric(ATG) - 1L) * (CR == "other_nrf"))

# To use
CoxFit <-  coxph(
  Surv(endpoint7, status2) ~
    ATG.1 + hirisk.1 +
    ATG.2 + hirisk.2 +
    ATG.3 + strata(CR),
  data = dat_cr,
  x = TRUE
)


coef(preDLI_cox)
coef(CoxFit_CR)

CR_forms <- list(
  "CD4_abs_log" = ~ value(CD4_abs_log):CR,
  "CD8_abs_log" = ~ value(CD8_abs_log):CR
)

# Remove also the hirisk three-way interaction
preDLI_long_corr_CD4_randint <- update(
  preDLI_long_corr_CD4,
  fixed. = . ~ ns(intSCT2_7, 3) * ATG + CMV_PD,
  random = ~ ns(intSCT2_7, 3) | IDAA
)

preDLI_long_corr_CD8_randint <- update(
  preDLI_long_corr_CD8,
  fixed. = . ~ ns(intSCT2_7, 3) * ATG + CMV_PD,
  random = ~ ns(intSCT2_7, 3) | IDAA
)

# Remember fixed intercept model is not possible!!
# Also make sure to add less basehaz knots
jFit_CR <- jm(
  CoxFit,
  list(preDLI_long_corr_CD4_randint),
  #list(preDLI_long_corr_CD4_randint, preDLI_long_corr_CD8_randint),
  time_var = "intSCT2_7",
  #functional_forms = list("CD4_abs_log" = ~ value(CD4_abs_log) * CR - 1),
  functional_forms = list(
    #"CD4_abs_log" = ~ value(CD4_abs_log):CR
    "CD4_abs_log" = ~ value(CD4_abs_log) + value(CD4_abs_log):CR - 1
  ),
  n_burnin = 10,
  n_iter = 25000,
  base_hazard_segments = 4,
  Bsplines_degree = 3,
  n_chains = 1#2
)
# see https://github.com/drizopoulos/JMbayes2/issues/31
summary(jFit_CR)
summary(jFit_CR)$Survival
summary(tar_read(preDLI_JM_value_corr_CD4))$`CoefTable-Event`[seq_len(8), ]
ggtraceplot(jFit_CR, c("betas"), grid = TRUE, gridrows = 6, gridcols = 3)
ggtraceplot(jFit_CR, c("gammas"), grid = TRUE, gridrows = 6, gridcols = 3)
ggtraceplot(jFit_CR, c("alphas"), grid = TRUE, gridrows = 2, gridcols = 3)


summary(jFit_CR)
