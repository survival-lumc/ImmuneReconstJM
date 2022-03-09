# Do DLI testing here and run for tomorrow
tar_meta(fields = c(name, time, seconds, warnings)) |>  View()
tar_load(
  c(
    dat_merged,
    preDLI_CD4_JM_value_corr,
    preDLI_CD4_JM_value_indep,
    preDLI_CD3_JM_value_corr,
    preDLI_CD3_JM_value_indep,
    preDLI_CD8_JM_value_corr,
    preDLI_CD8_JM_value_indep,
  )
)

anova(preDLI_CD3_JM_value_indep, preDLI_CD3_JM_value_corr)
anova(preDLI_CD4_JM_value_indep, preDLI_CD4_JM_value_corr)
anova(preDLI_CD8_JM_value_indep, preDLI_CD8_JM_value_corr)

# Check the correlated ones
preDLI_CD3_JM_value_corr |>  summary()
preDLI_CD4_JM_value_corr |>  summary() # CD4 has issues with correlated reffs
preDLI_CD8_JM_value_corr |>  summary()
plot(preDLI_CD4_JM_value_indep, 1)
plot(preDLI_CD4_JM_value_corr, 1)



# Post-DLI testing --------------------------------------------------------

theme_set(theme_bw(base_size = 14))
tar_load(NMA_postDLI_datasets)

dat_long <- NMA_postDLI_datasets$long
dat_wide <- NMA_postDLI_datasets$wide
#ids_oneobs <- dat_long[, .N, by = IDAA][N == 1][["IDAA"]]

ggplot(
  dat_long,
  aes(intDLI1, CD8_abs_log)
) +
  geom_point(
    size = 3.5, pch = 21,
    alpha = 0.8,
    col = "#359fda",
    fill = colorspace::lighten("#359fda", amount = 0.3)
  ) +
  facet_wrap(IDAA ~ sec_endpoint2_s)


# Try some model
postDLI_cox <- run_postDLI_cox(
  Surv(time, status) ~ ATG.1 + strata(trans), #ATG.2
  dat_wide = NMA_postDLI_datasets$wide#[!(IDAA %in% ids_oneobs)]
)
summary(postDLI_cox)

dat_long_bis <- cbind(dat_long, ns(dat_long$intDLI1, 3))
setnames(dat_long_bis, old = c("1", "2", "3"), new = paste0("basisf_", 1:3))

postDLI_CD4 <- run_preDLI_longitudinal(
  cell_line = "CD4_abs_log",
  form_fixed = "(basisf_1 + basisf_2 + basisf_3) * ATG + CMV_PD", # remove cmvpd
  form_random = list(IDAA = pdDiag(~ basisf_2)),
  #form_random = ~ ns(intDLI1, 3) | IDAA,
  dat = dat_long_bis#[!(IDAA %in% ids_oneobs)]
)

VarCorr(postDLI_CD4)

# Long CD3 - change function name
postDLI_CD4 <- run_preDLI_longitudinal(
  cell_line = "CD4_abs_log",
  form_fixed = "ns(intDLI1, 1) * ATG + CMV_PD", # remove cmvpd
  #form_random = list(IDAA = pdDiag(~ ns(intDLI1, 1))),
  form_random = ~ ns(intDLI1, 1) | IDAA,
  dat = NMA_postDLI_datasets$long#[!(IDAA %in% ids_oneobs)]
)


plot(postDLI_CD4)
VarCorr(postDLI_CD4)
dat_long[, preds_ind := fitted(postDLI_CD4)]

ggplot(
  dat_long,
  aes(intDLI1, CD4_abs_log)
) +
  geom_point(
    size = 3.5, pch = 21,
    alpha = 0.8,
    col = "#359fda",
    fill = colorspace::lighten("#359fda", amount = 0.3)
  ) +
  geom_line(aes(y = preds_ind), size = 1.25) +
  facet_wrap(~ IDAA)


# try JM
postDLI_CD4_JM <- jointModel(
  lmeObject = postDLI_CD4,
  survObject = postDLI_cox,
  CompRisk = TRUE,
  parameterization = "value",
  interFact = list("value" = ~ strata(trans) - 1),
  method = "spline-PH-aGH",
  timeVar = "intDLI1",
  lng.in.kn = 2,
  #iter.EM = 1000,
  verbose = TRUE
)

summary(postDLI_CD4_JM)
postDLI_CD4_JM$iters
matrixcalc::is.positive.definite(postDLI_CD4_JM$Hessian)
plot(postDLI_CD4_JM, 1)

newdat_postDLI <- expand.grid(
  "ATG" = levels(dat_long$ATG),
  "CMV_PD" = levels(dat_long$CMV_PD),
  "intDLI1" = seq(0, 3, by = 0.01)
)


predict(
  postDLI_CD4_JM,
  newdata = newdat_postDLI,
  interval = "confidence",
  type = "Marginal",
  returnData = TRUE
) |>
  ggplot(aes(intDLI1, pred, col = ATG, group = ATG)) +
  geom_ribbon(
    aes(ymin = low, ymax = upp),
    fill = "gray",
    alpha = 0.5,
    col = NA
  ) +
  geom_line(size = 1.25) +
  scale_y_continuous(
    breaks = log(c(5, 25, 100, 500, 1500)),
    labels = c(5, 25, 100, 500, 1500)
  ) +
  coord_cartesian(ylim = c(0, 8)) +
  facet_wrap(~ CMV_PD)


dat_long[, "fitted_JM" := fitted(
  postDLI_CD4_JM,
  process = "Longitudinal",
  type = "Subject"
)]


ggplot(
  dat_long,
  aes(intDLI1, CD4_abs_log)
) +
  geom_point(
    size = 3.5, pch = 21,
    alpha = 0.8,
    col = "#359fda",
    fill = colorspace::lighten("#359fda", amount = 0.3)
  ) +
  geom_line(aes(y = fitted_JM), size = 1.25) +
  geom_point(aes(y = fitted_JM), size = 1.5) +
  facet_wrap(~ IDAA)


endpoints_df <- dat_long[order(intDLI1), .SD[.N], by = IDAA]
dat_long[, .(.N, "endp" = unique(sec_endpoint2_s)), by = IDAA]

ggplot(dat_long, aes(intDLI1, fitted_JM, group = IDAA)) +
  geom_line(alpha = 0.75, size = 1, aes(col = IDAA)) +
  facet_grid(sec_endpoint2_s ~ ATG) +
  geom_point(data = endpoints_df, aes(sec_endpoint2, CD4_abs_log, shape = sec_endpoint2_s)) +
  theme(legend.position = "none")



# Testing basis 4 ---------------------------------------------------------


mods_cd4 <- lapply(1:3, function(df) {

  form_random_indep <- as.formula(paste0("~ ns(intDLI1, ", df, ")"))
  form_random_corr <- as.formula(paste0("~ ns(intDLI1, ", df, ") | IDAA"))

  mod_indep <- run_preDLI_longitudinal(
    cell_line = "CD4_abs_log",
    form_fixed = "ns(intDLI1, 2) * ATG + CMV_PD", # remove cmvpd
    form_random = list(IDAA = pdDiag(form_random_indep)),
    dat = NMA_postDLI_datasets$long#[!(IDAA %in% ids_oneobs)]
  )

  mod_corr <- run_preDLI_longitudinal(
    cell_line = "CD4_abs_log",
    form_fixed = "ns(intDLI1, 2) * ATG + CMV_PD", # remove cmvpd
    form_random = form_random_corr,
    dat = NMA_postDLI_datasets$long#[!(IDAA %in% ids_oneobs)]
  )

  list("indep" = mod_indep, "corr" = mod_corr)
})

names(mods_cd4) <- as.character(1:3)
mods_cd4$`1`$indep

# Compare independent mods
anova(mods_cd4$`1`$indep, mods_cd4$`2`$indep, mods_cd4$`3`$indep)
anova(mods_cd4$`1`$corr, mods_cd4$`2`$corr, mods_cd4$`3`$corr)


mods_cd3 <- lapply(1:3, function(df) {

  form_random_indep <- as.formula(paste0("~ ns(intDLI1, ", df, ")"))
  form_random_corr <- as.formula(paste0("~ ns(intDLI1, ", df, ") | IDAA"))

  mod_indep <- run_preDLI_longitudinal(
    cell_line = "CD3_abs_log",
    form_fixed = "ns(intDLI1, 2) * ATG + CMV_PD", # remove cmvpd
    form_random = list(IDAA = pdDiag(form_random_indep)),
    dat = NMA_postDLI_datasets$long#[!(IDAA %in% ids_oneobs)]
  )

  mod_corr <- run_preDLI_longitudinal(
    cell_line = "CD3_abs_log",
    form_fixed = "ns(intDLI1, 2) * ATG + CMV_PD", # remove cmvpd
    form_random = form_random_corr,
    dat = NMA_postDLI_datasets$long#[!(IDAA %in% ids_oneobs)]
  )

  list("indep" = mod_indep, "corr" = mod_corr)
})
names(mods_cd3) <- as.character(1:3)


mods_cd8 <- lapply(1:3, function(df) {

  form_random_indep <- as.formula(paste0("~ ns(intDLI1, ", df, ")"))
  form_random_corr <- as.formula(paste0("~ ns(intDLI1, ", df, ") | IDAA"))

  mod_indep <- run_preDLI_longitudinal(
    cell_line = "CD8_abs_log",
    form_fixed = "ns(intDLI1, 2) * ATG + CMV_PD", # remove cmvpd
    form_random = list(IDAA = pdDiag(form_random_indep)),
    dat = NMA_postDLI_datasets$long#[!(IDAA %in% ids_oneobs)]
  )

  mod_corr <- run_preDLI_longitudinal(
    cell_line = "CD8_abs_log",
    form_fixed = "ns(intDLI1, 2) * ATG + CMV_PD", # remove cmvpd
    form_random = form_random_corr,
    dat = NMA_postDLI_datasets$long#[!(IDAA %in% ids_oneobs)]
  )

  list("indep" = mod_indep, "corr" = mod_corr)
})
names(mods_cd8) <- as.character(1:3)


cd8_a <-  run_preDLI_longitudinal(
  cell_line = "CD8_abs_log",
  form_fixed = "ns(intDLI1, 1) * ATG + CMV_PD", # remove cmvpd
  #form_random = ~ ns(intDLI1, 1) | IDAA,
  form_random = list(IDAA = pdDiag(~ 1)),
  dat = NMA_postDLI_datasets$long
)

cd8_c <-  run_preDLI_longitudinal(
  cell_line = "CD8_abs_log",
  form_fixed = "ns(intDLI1, 1) * ATG + CMV_PD", # remove cmvpd
  #form_random = ~ ns(intDLI1, 1) | IDAA,
  form_random = list(IDAA = pdDiag(~ ns(intDLI1, 1))),
  dat = NMA_postDLI_datasets$long
)

anova(cd8_a, cd8_c)

anova(cd8_noint, mods_cd8$`2`$indep)

#form_random = list(IDAA = pdDiag(~ ns(intDLI1, 1))),
#form_random = ,

anova(mods_cd8$`1`$indep, mods_cd8$`2`$indep, mods_cd8$`3`$indep)
anova(mods_cd8$`1`$corr, mods_cd8$`2`$corr, mods_cd8$`3`$corr)

anova(mods_cd3$`1`$indep, mods_cd3$`2`$indep, mods_cd3$`3`$indep)
anova(mods_cd3$`1`$corr, mods_cd3$`2`$corr, mods_cd3$`3`$corr)

anova(mods_cd4$`1`$indep, mods_cd4$`2`$indep, mods_cd4$`3`$indep)
anova(mods_cd4$`1`$corr, mods_cd4$`2`$corr, mods_cd4$`3`$corr)
# Use only two dfs



# Tests -------------------------------------------------------------------


library(ggplot2)
library(data.table)

tar_load(c(postDLI_cox, postDLI_CD4_long_corr, NMA_postDLI_datasets))
tar_load(c(postDLI_cox, postDLI_CD8_long_corr, NMA_postDLI_datasets))

# CD8 only random int
cd8_long <- run_preDLI_longitudinal(
  cell_line = "CD8_abs_log",
  form_fixed = "ns(intDLI1, 2) * ATG + CMV_PD", # remove cmvpd
  form_random = ~ ns(intDLI1, 2) | IDAA,
  #form_random = list(IDAA = pdDiag(~ ns(intDLI1, 2))),
  dat = NMA_postDLI_datasets$long
)

dat_long <- NMA_postDLI_datasets$long

postDLI_CD8_JM_value_corr <- jointModel(
  lmeObject = cd8_long,
  survObject = postDLI_cox,
  CompRisk = TRUE,
  parameterization = "value",
  interFact = list("value" = ~ strata(trans) - 1),
  method = "spline-PH-aGH",
  timeVar = "intDLI1",
  lng.in.kn = 3,
  #iter.EM = 500,
  verbose = TRUE
)


tar_load(NMA_postDLI_datasets)



newdat_postDLI <- expand.grid(
  "ATG" = levels(dat_long$ATG),
  "CMV_PD" = levels(dat_long$CMV_PD),
  "intDLI1" = seq(0, 3, by = 0.01)
)


predict(
  postDLI_CD8_JM_value_corr,
  newdata = newdat_postDLI,
  interval = "confidence",
  type = "Marginal",
  returnData = TRUE
) |>
  ggplot(aes(intDLI1, pred, col = ATG, group = ATG)) +
  geom_ribbon(
    aes(ymin = low, ymax = upp),
    fill = "gray",
    alpha = 0.5,
    col = NA
  ) +
  geom_line(size = 1.25) +
  scale_y_continuous(
    breaks = log(c(5, 25, 100, 500, 1500)),
    labels = c(5, 25, 100, 500, 1500)
  ) +
  coord_cartesian(ylim = c(0, 8)) +
  facet_wrap(~ CMV_PD)


dat_long[, "fitted_JM" := fitted(
  postDLI_CD8_JM_value_corr,
  process = "Longitudinal",
  type = "Subject"
)]


ggplot(
  dat_long,
  aes(intDLI1, CD8_abs_log)
) +
  geom_point(
    size = 3.5, pch = 21,
    alpha = 0.8,
    col = "#359fda",
    fill = colorspace::lighten("#359fda", amount = 0.3)
  ) +
  geom_line(aes(y = fitted_JM), size = 1.25) +
  geom_point(aes(y = fitted_JM), size = 1.5) +
  facet_wrap(~ IDAA)7



# Plot all three cell lines -----------------------------------------------


dat_long |>
  melt.data.table(
    measure.vars = patterns("*_log$"),
    variable.name = "cell_line",
    value.name = "count"
  ) |>
  ggplot(aes(intDLI1, count, col = cell_line)) +
  geom_point(aes(shape = cell_line), size = 4, alpha = 0.75) +
  facet_wrap(sec_endpoint2_s ~ IDAA) +
  scale_y_continuous(
    breaks = log(c(5, 25, 100, 500, 1500)),
    labels = c(5, 25, 100, 500, 1500)
  ) +
  theme_bw(base_size = 14)


dat_long |>
  melt.data.table(
    measure.vars = patterns("*_log$"),
    variable.name = "cell_line",
    value.name = "count"
  ) |>
  ggplot(aes(intDLI1, count, col = ATG, group = IDAA)) +
  geom_line(size = 1) +
  scale_y_continuous(
    breaks = log(c(5, 25, 100, 500, 1500)),
    labels = c(5, 25, 100, 500, 1500)
  ) +
  facet_grid(sec_endpoint2_s ~ cell_line) +
  theme_bw(base_size = 14)
