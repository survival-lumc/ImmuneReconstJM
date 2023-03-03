tar_load(
  c(
    preDLI_JM_value_corr_CD3_inter,
    preDLI_JM_value_corr_CD4_inter,
    preDLI_JM_value_corr_CD8_inter,

  )
)

summary(preDLI_JM_value_corr_CD4_inter)
summary(preDLI_JM_value_corr_CD3_inter)
summary(preDLI_JM_value_corr_CD8_inter)

mod_gvh_rd <- coxph(
  Surv(endpoint7, endpoint7_s == "gvhd") ~ 1,
  data = dat_wide_pre[ATG == "RD"],
  x = TRUE, model = TRUE
)

mod_gvh_ud <- coxph(
  Surv(endpoint7, endpoint7_s == "gvhd") ~ 1,
  data = dat_wide_pre[ATG == "UD(+ATG)"],
  x = TRUE, model = TRUE
)

# RUn longitudinal in each subset
tar_load(preDLI_long_corr_CD4)
long_rd <- update(
  preDLI_long_corr_CD4,
  fixed = . ~ ns(intSCT2_7, 3) * hirisk + CMV_PD,
  data = dat_long_pre[ATG == "RD"]
)

long_ud <- update(
  preDLI_long_corr_CD4,
  fixed = . ~ ns(intSCT2_7, 3) * hirisk + CMV_PD,
  data = dat_long_pre[ATG == "UD(+ATG)"]
)

mod_rd <- jointModel(
  lmeObject = long_rd,
  survObject = mod_gvh_rd,
  method = "spline-PH-aGH",
  timeVar = "intSCT2_7",
  parameterization = "value",
  #interFact = list("value" = ~ strata(trans) - 1),
  iter.EM = 500L,
  #iter.qN = 1000L,
  lng.in.kn = 3L,
  numeriDeriv = "cd",
  eps.Hes = 1e-04,
  verbose = TRUE
)

mod_rd |>  summary()

mod_ud <- jointModel(
  lmeObject = long_ud,
  survObject = mod_gvh_ud,
  method = "spline-PH-aGH",
  timeVar = "intSCT2_7",
  parameterization = "value",
  #interFact = list("value" = ~ strata(trans) - 1),
  iter.EM = 500L,
  #iter.qN = 1000L,
  lng.in.kn = 3L,
  numeriDeriv = "cd",
  eps.Hes = 1e-04,
  verbose = TRUE
)



# Trying tdep Cox ---------------------------------------------------------


tar_load(c(NMA_preDLI_datasets))

dat_long <- data.frame(
  NMA_preDLI_datasets$long[, c(
    "IDAA",
    paste0("CD", c(3, 4, 8), "_abs_log"),
    "intSCT2_7"
  )]
)

dat_wide <- data.frame(
  NMA_preDLI_datasets$wide[, c(
    "IDAA",
    "hirisk",
    "ATG",
    "CMV_PD",
    "endpoint7_s",
    "endpoint7"
  )]
)

temp <- tmerge(
  dat_wide,
  dat_wide,
  id = IDAA,
  gvhd_ind = event(endpoint7, as.numeric(endpoint7_s == "gvhd")),
  relapse_ind = event(endpoint7, as.numeric(endpoint7_s == "relapse")),
  other_nrf_ind = event(endpoint7, as.numeric(endpoint7_s == "other_nrf"))
)
temp_updt <- tmerge(
  temp, dat_long,
  id = IDAA,
  log_CD3 = tdc(intSCT2_7, CD3_abs_log, init = log(0.25)), #add rest
  log_CD4 = tdc(intSCT2_7, CD4_abs_log, init = log(0.25)),
  log_CD8 = tdc(intSCT2_7, CD8_abs_log, init = log(0.25))
)
head(temp_updt, 10)

mod_cd4 <- coxph(
  Surv(tstart, tstop, other_nrf_ind) ~ ATG + hirisk + log_CD4 + log_CD8,
  data = temp_updt
)
coef(mod_cd4)

mod_cd8 <- coxph(
  Surv(tstart, tstop, gvhd_ind) ~ ATG + hirisk + log_CD8,
  data = temp_updt
)

summary(tar_read(preDLI_JM_value_corr_CD4))$`CoefTable-Event`[c(1, 2, 6), ]


temp_updt |>
  ggplot(aes(tstop, log_CD4, group = IDAA)) +
  geom_step() +
  facet_grid(hirisk * ATG ~ CMV_PD) +
  theme_light()


#NOCB

temp_updt <- tmerge(
  temp, dat_long,
  id = IDAA,
  log_CD3 = tdc(intSCT2_7, CD3_abs_log), #add rest
  log_CD4 = tdc(intSCT2_7, CD4_abs_log),
  log_CD8 = tdc(intSCT2_7, CD8_abs_log)
)
