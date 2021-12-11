tar_load(c("datasets", "dat_merged", "reference_values"))
theme_set(theme_bw(base_size = 14)) # add also to ebmt script

# Prepare wide and long datasets
dat_wide <- datasets$wide
dat_long <- datasets$long

dat_wide[, earlylow_DLI := factor(earlylow_DLI, labels = c("no", "yes"))]
dat_long[, earlylow_DLI := factor(earlylow_DLI, labels = c("no", "yes"))]

# Try the one dataset approach??

# Question 1 --------------------------------------------------------------

# 6 months post HSCT
admin_cens <- 6
dat_long_q1 <- copy(dat_long)
dat_wide_q1 <- copy(dat_wide)

dat_wide_q1[endpoint7 >= admin_cens, ':=' (
  endpoint7 = admin_cens,
  endpoint7_s = "cens"
)]

dat_long_q1[endpoint7 >= admin_cens, ':=' (
 endpoint7 = admin_cens,
 endpoint7_s = "cens"
)]

# Keep measurements prior to endpoint - first measurments at endpoint are taken JUST  before
dat_long_q1[intSCT2_5 == endpoint7, intSCT2_5 := intSCT2_5 - 0.01]
dat_long_q1 <- dat_long_q1[intSCT2_5 < endpoint7]
dat_wide_q1 <- dat_wide_q1[IDAA %in% unique(dat_long_q1$IDAA)]

# Only 7 patients have hd dli as endpoint (rest are modif t-cell prods), all close to 6 mos
#.. just censor them? Mix with nrf_other?
dat_wide[endpoint7_s == "hd_dli" & endpoint_specify7 != "modified T-cell product",
         c("endpoint7_s", "endpoint7")][endpoint7 < 6]

# Make raw plots by ITT
dat_long_q1 |>
  ggplot(aes(intSCT2_5, CD3_abs_log)) +
  geom_line(aes(group = IDAA)) +
  facet_grid(CMV_PD * ATG ~ hirisk) +
  theme(legend.position = "none") +
  scale_y_continuous(
    breaks = log(c(1, 5, 25, 100, 500, 2000)),
    labels = c(1, 5, 25, 100, 500, 2000)
  ) +
  labs(x = "Time since alloSCT (months)", y = "CD3 cell counts")


# Prep submodels ----------------------------------------------------------

dat_wide_q1[endpoint_specify7 == "modified T-cell product", endpoint7_s := "cens"]
dat_long_q1[endpoint_specify7 == "modified T-cell product", endpoint7_s := "cens"]
dat_wide_q1[endpoint7_s == "hd_dli", endpoint7_s := "cens"]
dat_long_q1[endpoint7_s == "hd_dli", endpoint7_s := "cens"]

dat_wide_q1[, cr_new := droplevels(endpoint7_s)]
dat_long_q1[, cr_new := droplevels(endpoint7_s)]


table(dat_wide_q1$cr_new)

tmat <- trans.comprisk(K = 3, names = c("gvhd", "relapse", "other_nrf"))
covs <- c("CMV_PD", "hirisk", "ATG", "earlylow_DLI")

msdat <- msprep(
  time = c(NA, rep("endpoint7", 3)),
  status = with(
    dat_wide_q1, cbind(
      NA,
      1 * (cr_new == "gvhd"),
      1 * (cr_new == "relapse"),
      1 * (cr_new == "other_nrf")
    )
  ),
  data = data.frame(dat_wide_q1),
  trans = tmat,
  keep = covs,
  id = "IDAA"
)

table(dat_wide_q1$endpoint_specify7, dat_wide_q1$cr_new)

msdat_expand <- expand.covs(msdat, covs, append = TRUE, longnames = FALSE)

# Survival submodel
mod_comp <- coxph(
  Surv(time, status) ~ hirisk.1 + # gvhd, normally without
    ATG.2 + hirisk.2 + # relapse
    ATG.3 + # other_nrf
    strata(trans),
  cluster = IDAA,
  model = TRUE,
  x = TRUE,
  data = msdat_expand
)

mod_comp
# mod_comp <- coxph(
#   Surv(time, status) ~ hirisk.1 + # gvhd
#     ATG.2 + hirisk.2 + # relapse
#     ATG.3 + hirisk.3 + # other_nrf
#     strata(trans),
#   cluster = IDAA,
#   model = TRUE,
#   x = TRUE,
#   data = msdat_expand
# )

mod_comp$model$trans1 <- as.numeric(mod_comp$model$`strata(trans)` == "trans=1")
mod_comp$model$trans2 <- as.numeric(mod_comp$model$`strata(trans)` == "trans=2")
mod_comp$model$trans3 <- as.numeric(mod_comp$model$`strata(trans)` == "trans=3")


bh <- data.table(basehaz(mod_comp, centered = FALSE))

bh[order(strata, time)][, .(
  "time" = time[-length(time)],
  "hazard" = diff(hazard)
), by = strata]

bh|>
  ggplot(aes(time, hazard)) +
  geom_step(size = 1) +
  facet_wrap(~ strata, labeller = labeller(strata = c(
    "trans=1" = "1: gvhd",
    "trans=2" = "2: relapse",
    "trans=3" = "3: other_nrf"
  ))) +
  labs(x = "Time since alloSCT (months)", y = "Baseline hazard")

summary(mod_comp)

lmeFit <- lme( # * hirisk model
  fixed = CD3_abs_log ~ ns(intSCT2_5, 3) * hirisk + ATG + CMV_PD, # should be DLI type and ATG
  #random = ~ intSCT2_5_reset | IDAA,
  random = list(IDAA = pdDiag(~ ns(intSCT2_5, 3))),
  control = lmeControl(opt = "optim", msMaxIter = 200),
  data = data.frame(dat_long_q1)
)

summary(lmeFit)

# Check lmeFit - marginal
newdat_jm <- expand.grid(
  "ATG" = levels(dat_long_q1$ATG),
  "CMV_PD" = levels(dat_long_q1$CMV_PD),
  "hirisk" = levels(dat_long_q1$hirisk),
  "intSCT2_5" = seq(0, 6, by = 0.05)
)

df_preds_marg <- cbind(
  newdat_jm,
  "pred" = predict(lmeFit, newdata = newdat_jm, level = 0L)
)

df_preds_marg |>
  ggplot(aes(intSCT2_5, pred,
             col = interaction(hirisk, ATG),
             group = interaction(hirisk, ATG))) +
  geom_line(size = 1.5) +
  facet_wrap(~ CMV_PD)
  #geom_line(aes(col = hirisk, group = hirisk)) +
  #facet_wrap(CMV_PD ~ ATG)



summary(lmeFit)

dform <- list(
  fixed = ~ 0 + dns(intSCT2_5, 3),
  random = ~ 0 + dns(intSCT2_5, 3),
  indFixed = c(2:4),
  indRandom = c(2:4)
)

# Covnerges with current values - 8 minutes w/ slopes
jm_predli <- jointModel(
  lmeObject = lmeFit,
  survObject = mod_comp,
  timeVar = "intSCT2_5",
  method = "spline-PH-aGH",# or piecewise
  derivForm = dform,
  parameterization = "both", # slope
  interFact = list(
    "value" = ~ strata(trans) - 1,
    "slope" = ~ strata(trans) - 1
  ),
  # interFact = list(
  #   "slope" = ~ trans1 + trans1:hirisk.1 + #need to add hirisk to actual model
  #     trans2 + trans2:hirisk.2 +
  #     trans3 + trans3:hirisk.3 - 1
  # ),
  CompRisk = TRUE,
  iter.EM = 500 # and check and jmbayes2
)

# Interaction with hirisk?

summary(jm_predli)

jm_predli$coefficients
mod_comp$coefficients

reffs <- coefficients(jm_predli)
jm_predli
# Do da slope?? # eventually get baseline hazards


#

predict(
  jm_predli,
  newdata = newdat_jm,
  type = "Marginal",
  idVar = "IDAA",
  returnData = TRUE,
  interval = "confidence"
) |>
  ggplot(aes(x = intSCT2_5, y = pred,
             group = interaction(ATG, hirisk), col = interaction(ATG, hirisk))) +
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
  xlim(c(0, 6))

# Three-way interaction actually converges!


# Model from DLI ----------------------------------------------------------

# Get predicted random effects and add to data?
# Try adding offset in simple mixed model?

dat_long_dli <- copy(dat_long)[uDLI_s == "uDLI"]
dat_wide_dli <- copy(dat_wide)[uDLI_s == "uDLI"]

dat_wide_dli[#endpoint6 >= 18 &
  endpoint6_s %in% c("gvhd", "other_nrf")][, c("endpoint6", "endpoint6_s")]

dat_wide_dli$endpoint6_s |>  table()

admin_cens_dli <- 18
dat_wide_dli[endpoint6 >= admin_cens_dli, ':=' (
  endpoint6 = admin_cens_dli,
  endpoint6_s = "cens"
)]

dat_wide_dli$endpoint6_s |>  table()


dat_long_dli[endpoint6 >= admin_cens_dli, ':=' (
  endpoint6 = admin_cens_dli,
  endpoint6_s = "cens"
)]

table(dat_wide_dli$endpoint6_s)

dat_wide_dli <- cbind(
  dat_wide_dli,
  cd3_pred_tdli = predict(
    jm_predli,
    newdata = dat_wide_dli[, c("uDLI", "hirisk", "ATG", "CMV_PD")][, intSCT2_5 := uDLI][],
    type = "Marginal"
  )
)
dat_wide_dli$endpoint6_s |>  table()


# Define the competing endpoints
dat_long_dli[, cr_new := factor(
  fcase(
    endpoint6_s %in% c("relapse", "other_nrf"), "rel_other_nrf",
    endpoint6_s == "gvhd", "gvhd",
    endpoint6_s == "cens", "cens"
  ), levels = c("cens", "gvhd", "rel_other_nrf")
)]
dat_wide_dli[, cr_new := factor(
  fcase(
    endpoint6_s %in% c("relapse", "other_nrf"), "rel_other_nrf",
    endpoint6_s == "gvhd", "gvhd",
    endpoint6_s == "cens", "cens"
  ), levels = c("cens", "gvhd", "rel_other_nrf")
)]

table(dat_wide_dli$cr_new)

# Make sure to clock-reset measurements!
dat_long_prepped <- dat_long_dli[uDLI < intSCT2_5]
dat_long_prepped[, intSCT2_5_reset := intSCT2_5 - uDLI]
dat_wide_prepped <- dat_wide_dli[IDAA %in% unique(dat_long_prepped$IDAA)]

#
tmat <- trans.comprisk(K = 2, names = c("first_dli", "gvhd", "rel_other_nrf"))
covs <- c("CMV_PD", "hirisk", "ATG", "earlylow_DLI", "cd3_pred_tdli")

msdat <- msprep(
  time = c(NA, rep("endpoint6", 2)),
  status = with(
    dat_wide_prepped, cbind(
      NA,
      1 * (cr_new == "gvhd"),
      1 * (cr_new == "rel_other_nrf")
    )
  ),
  data = data.frame(dat_wide_prepped),
  start = list(state = rep(1, nrow(dat_wide_prepped)), time = dat_wide_prepped$uDLI),
  trans = tmat,
  keep = covs,
  id = "IDAA"
)

msdat_expand <- expand.covs(msdat, covs, append = TRUE, longnames = FALSE)

#rel_other_nrf
coxph(
  Surv(endpoint6 - uDLI, cr_new == "rel_other_nrf") ~ earlylow_DLI +
    offset(cd3_pred_tdli),
  data = data.frame(dat_wide_prepped)
)

# Survival submodel (excluding pre-DLI info)
mod_comp <- coxph(
  Surv(time, status) ~
    earlylow_DLI.1 + cd3_pred_tdli.1 + # gvhd
    hirisk.2 + cd3_pred_tdli.2 + # rel_other_nrf
    strata(trans),
  cluster = IDAA,
  model = TRUE,
  x = TRUE,
  data = msdat_expand
)

summary(mod_comp)

# All 14 get early low-dose
table(dat_wide_prepped[cr_new == "gvhd"]$earlylow_DLI)
table(dat_wide_prepped[cr_new == "rel_other_nrf"]$cd3_pred_tdli)
table(dat_wide_prepped[cr_new == "rel_other_nrf"]$earlylow_DLI)

# is uDLI only low-dose??

# Include previous CD3 here? And CMV?
lmeFit <- lme(
  fixed = CD3_abs_log ~ ns(intSCT2_5_reset, 2) + earlylow_DLI + ATG,
  #random = ~ ns(intSCT2_5_reset, 2)| IDAA,
  random = list(IDAA = pdDiag(~ ns(intSCT2_5_reset, 2))),
  control = lmeControl(opt = "optim", msMaxIter = 200),
  data = data.frame(dat_long_prepped)
)

summary(lmeFit)


dform <- list(
  fixed = ~ 0 + dns(intSCT2_5_reset, 2),
  random = ~ 0 + dns(intSCT2_5_reset, 2),
  indFixed = c(2, 3),
  indRandom = c(2, 3)
)

jm_postdli_slope <- jointModel(
  lmeObject = lmeFit,
  survObject = mod_comp,
  CompRisk = TRUE,
  timeVar = "intSCT2_5_reset",
  method = "spline-PH-aGH",
  derivForm = dform,
  interFact = list(
    "value" = ~ strata(trans) - 1#,
    #"slope" = ~ strata(trans) - 1
  ),
  parameterization = "value",
  control = list("iter.EM" = 200)
)

table(dat_wide_prepped$cr_new, dat_wide_prepped$earlylow_DLI)
table(dat_wide_prepped$cr_new)
summary(jm_postdli_slope)

# To-do
# - Stijn: try with Cd8
# - try fitting all with jmbayes2
# - What is interpretation with only slope?
