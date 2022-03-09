# Thinking about longitudinal submodel ------------------------------------


# First check convergence errors Tar meta
tar_meta() |> View()

tar_load(
  c(
    NMA_preDLI_datasets,
    preDLI_CD3_long,
    preDLI_CD8_jointModel_both,
    preDLI_CD3_jointModel_corr
  )
)

dat_wide <- NMA_preDLI_datasets$wide
dat_long <- NMA_preDLI_datasets$long

# Source support functions
source("data-raw/prepare-raw-data.R")
source("R/modelling-helpers.R")
source("R/plotting-helpers.R")

# For marginal fit
newdat_preDLI <- expand.grid(
  "ATG" = levels(dat_long$ATG),
  "CMV_PD" = levels(dat_long$CMV_PD),
  "hirisk" = levels(dat_long$hirisk),
  "intSCT2_5" = seq(0, 6, by = 0.05)
)

set.seed(684646)
ID_samps <- droplevels(sample(unique(dat_long$IDAA), size = 56, replace = FALSE))


# First checks of improvement with three-way interaction ------------------


mod_main <- run_preDLI_longitudinal(
  cell_line = "CD4_abs_log",
  form_fixed = "ns(intSCT2_5, 3) + ATG + hirisk + CMV_PD",
  form_random = list("IDAA" = pdDiag(~ ns(intSCT2_5, 3))),
  dat = NMA_preDLI_datasets$long
)

mod_twoway <- run_preDLI_longitudinal(
  cell_line = "CD4_abs_log",
  form_fixed = "ns(intSCT2_5, 3) * hirisk + ATG + CMV_PD",
  form_random = list("IDAA" = pdDiag(~ ns(intSCT2_5, 3))),
  dat = NMA_preDLI_datasets$long
)

mod_twoway_ATG <- run_preDLI_longitudinal(
  cell_line = "CD4_abs_log",
  form_fixed = "ns(intSCT2_5, 3) * ATG + hirisk + CMV_PD",
  form_random = list("IDAA" = pdDiag(~ ns(intSCT2_5, 3))),
  dat = NMA_preDLI_datasets$long
)

mod_threeway <- run_preDLI_longitudinal(
  cell_line = "CD4_abs_log",
  form_fixed = "ns(intSCT2_5, 3) * hirisk * ATG + CMV_PD",
  form_random = list("IDAA" = pdDiag(~ ns(intSCT2_5, 3))),
  dat = NMA_preDLI_datasets$long
)
VarCorr(mod_threeway)



anova(update(mod_main, method = "ML"), update(mod_twoway_ATG, method = "ML"))
anova(update(mod_main, method = "ML"), update(mod_twoway, method = "ML"))
anova(update(mod_twoway_ATG, method = "ML"), update(mod_threeway, method = "ML"))


# Plot marginal of two way - you see divergence but it is non-signif
cbind(
  newdat_preDLI,
  "pred" = predict(mod_twoway, newdata = newdat_preDLI, level = 0L)
) |>
  ggplot(
    aes(
      x = intSCT2_5,
      y = pred,
      group = interaction(ATG, hirisk),
      col = ATG
    )
  ) +
  geom_line(
    aes(linetype = hirisk),
    size = 1.5,
    alpha = 0.75
  ) +
  # add repel labels
  facet_wrap(facets = ~ CMV_PD) +
  labs(
    x = "Time since alloHCT (months)",
    y = "CD3 cell counts",
    linetype = "ITT",
    col = "Donor type"
  ) +
  scale_y_continuous(
    breaks = log(c(5, 25, 100, 500, 1500)),
    labels = c(5, 25, 100, 500, 1500)
  ) +
  xlim(c(0, 6)) +
  scale_color_manual(
    labels = c("Related", "Unrelated (+ATG)"),
    values = c("brown", "darkblue")
  ) +
  scale_linetype_manual(
    labels = c("No", "Yes (high risk)"),
    values = c("solid", "dotdash")
  )

summary(mod_twoway)
summary(mod_threeway)

cbind(mod_threeway_noint$coefficients$fixed,
mod_threeway_corr_noint$coefficients$fixed)


# Thinking of random effect structures  -----------------------------------


# Stick with three-way interaction
VarCorr(mod_threeway)

# Better in this case to omit random intercept as
mod_threeway_noint <- run_preDLI_longitudinal(
  cell_line = "CD4_abs_log",
  form_fixed = "ns(intSCT2_5, 3) * hirisk * ATG + CMV_PD",
  form_random = list("IDAA" = pdDiag(~ 0 + ns(intSCT2_5, 3))),
  dat = NMA_preDLI_datasets$long
)

VarCorr(mod_threeway_noint)
anova(update(mod_threeway, method = "ML"), update(mod_threeway_noint, method = "ML"))

# Let's try correlated random effects now
mod_threeway_corr <- run_preDLI_longitudinal(
  cell_line = "CD4_abs_log",
  form_fixed = "ns(intSCT2_5, 3) * hirisk * ATG + CMV_PD",
  form_random = ~ ns(intSCT2_5, 3) | IDAA,
  dat = NMA_preDLI_datasets$long
)

# Correlation between intercept and reff for second basis fun is problematic
# .. almost collinear. If covariance too large relative to variances: Hessian
# not positive definite
vcor <- VarCorr(mod_threeway_corr)
vcor

# Go back
sds <- as.numeric(vcor[1:4, "StdDev"])
corr_mat <- diag(1, nrow = 4)
corr_red <- matrix(as.numeric(vcor[2:4, 3:5]), nrow = 3)
corr_mat[lower.tri(corr_mat)] <- corr_mat[upper.tri(corr_mat)] <-
  corr_red[lower.tri(corr_red, diag = TRUE)]
outer(sds, sds) * corr_mat

# And now with fixed int
mod_threeway_corr_noint <- run_preDLI_longitudinal(
  cell_line = "CD4_abs_log",
  form_fixed = "ns(intSCT2_5, 3) * hirisk * ATG + CMV_PD",
  form_random = ~ 0 + ns(intSCT2_5, 3) | IDAA,
  dat = NMA_preDLI_datasets$long
)
VarCorr(mod_threeway_corr_noint)

# Let's check how this one fits vs indep
cbind(
  dat_long,
  "corr" = fitted(mod_threeway_corr_noint),
  "indep" = fitted(mod_threeway_noint)
)[IDAA %in% ID_samps] |>
  melt.data.table(
    id.vars = c("IDAA", "intSCT2_5", "CD4_abs_log"),
    measure.vars = c("corr", "indep"),
    variable.name = "reffs_type",
    value.name = "pred"
  ) |>
  ggplot(aes(intSCT2_5, CD4_abs_log)) +
  geom_point(col = "blue") +
  geom_line(aes(y = pred, linetype = reffs_type), size = 1, alpha = 0.5) +
  facet_wrap(~ IDAA)



# Check the marginal fits of both
cbind(
  newdat_preDLI,
  "pred" = predict(mod_threeway_corr_noint, newdata = newdat_preDLI, level = 0L)
) |>
ggplot(
  aes(
    x = intSCT2_5,
    y = pred,
    group = interaction(ATG, hirisk),
    col = ATG
  )
) +
  geom_line(
    aes(linetype = hirisk),
    size = 1.5,
    alpha = 0.75
  ) +
  # add repel labels
  facet_wrap(facets = ~ CMV_PD) +
  labs(
    x = "Time since alloHCT (months)",
    y = "CD3 cell counts",
    linetype = "ITT",
    col = "Donor type"
  ) +
  scale_y_continuous(
    breaks = log(c(5, 25, 100, 500, 1500)),
    labels = c(5, 25, 100, 500, 1500)
  ) +
  xlim(c(0, 6)) +
  scale_color_manual(
    labels = c("Related", "Unrelated (+ATG)"),
    values = c("brown", "darkblue")
  ) +
  scale_linetype_manual(
    labels = c("No", "Yes (high risk)"),
    values = c("solid", "dotdash")
  )

cbind(
  newdat_preDLI,
  "pred" = predict(mod_threeway_noint, newdata = newdat_preDLI, level = 0L)
) |>
  ggplot(
    aes(
      x = intSCT2_5,
      y = pred,
      group = interaction(ATG, hirisk),
      col = ATG
    )
  ) +
  geom_line(
    aes(linetype = hirisk),
    size = 1.5,
    alpha = 0.75
  ) +
  # add repel labels
  facet_wrap(facets = ~ CMV_PD) +
  labs(
    x = "Time since alloHCT (months)",
    y = "CD3 cell counts",
    linetype = "ITT",
    col = "Donor type"
  ) +
  scale_y_continuous(
    breaks = log(c(5, 25, 100, 500, 1500)),
    labels = c(5, 25, 100, 500, 1500)
  ) +
  xlim(c(0, 6)) +
  scale_color_manual(
    labels = c("Related", "Unrelated (+ATG)"),
    values = c("brown", "darkblue")
  ) +
  scale_linetype_manual(
    labels = c("No", "Yes (high risk)"),
    values = c("solid", "dotdash")
  )


#
dat_long[, .SD[1], by = IDAA] |>
  melt.data.table(
    measure.vars = patterns("*abs_log$"),
    variable.name = "cell_line",
    value.name = "count"
  ) |>
  ggplot(aes(intSCT2_5, exp(count))) +
  geom_point() +
  xlim(c(0, 6)) +
  ylim(c(0, 1500)) +
  facet_wrap(~ cell_line) +
  theme_bw()


# "Flattening of trajectory: extra basis function" ------------------------


CD4_ns4 <- run_preDLI_longitudinal(
  cell_line = "CD4_abs_log",
  form_fixed = "ns(intSCT2_5, 4) * ATG * hirisk + CMV_PD", #* hirisk + CMV_PD",
  form_random = list("IDAA" = pdDiag(~ 0 + ns(intSCT2_5, 4))),
  dat = NMA_preDLI_datasets$long
)

VarCorr(CD4_ns4)


cbind(
  newdat_preDLI,
  "pred" = predict(CD4_ns4, newdata = newdat_preDLI, level = 0L)
) |>
  ggplot(
    aes(
      x = intSCT2_5,
      y = pred,
      group = interaction(ATG, hirisk),
      col = ATG
    )
  ) +
  geom_line(
    aes(linetype = hirisk),
    size = 1.5,
    alpha = 0.75
  ) +
  # add repel labels
  facet_wrap(facets = ~ CMV_PD) +
  labs(
    x = "Time since alloHCT (months)",
    y = "CD3 cell counts",
    linetype = "ITT",
    col = "Donor type"
  ) +
  scale_y_continuous(
    breaks = log(c(5, 25, 100, 500, 1500)),
    labels = c(5, 25, 100, 500, 1500)
  ) +
  xlim(c(0, 6)) +
  scale_color_manual(
    labels = c("Related", "Unrelated (+ATG)"),
    values = c("brown", "darkblue")
  ) +
  scale_linetype_manual(
    labels = c("No", "Yes (high risk)"),
    values = c("solid", "dotdash")
  )


# Some simplified jm test -------------------------------------------------

# (Skip for this meeting, go to CD8)

system.time({
  jm1 <- jointModel(
    lmeObject = CD4_simplif,
    survObject = preDLI_cox,
    CompRisk = TRUE,
    method = "spline-PH-aGH",
    timeVar = "intSCT2_5",
    parameterization = "value",
    #iter.EM = 2000,
    interFact = list("value" = ~ strata(trans) - 1)
  )
})

summary(jm1)

# 36 mins
system.time({
  jm2 <- update(
    jm1,
    parameterization = "both",
    iter.EM = 500,
    parameterization = "both",
    interFact = list("value" = ~ strata(trans) - 1, "slope" = ~ strata(trans) - 1),
    derivForm = list(
      fixed = ~ 0 + dns(intSCT2_5, 4) +
        dns(intSCT2_5, 4):as.numeric(ATG == "ALT+ATG") +
        dns(intSCT2_5, 4):as.numeric(hirisk == "yes") +
        dns(intSCT2_5, 4):as.numeric(ATG == "ALT+ATG"):as.numeric(hirisk == "yes"),
      random = ~ 0 + dns(intSCT2_5, 4),
      indFixed = c(2:5, 9:16, 18:21),
      indRandom = c(1:4)
    )
  )
})

anova(jm1, jm2)



newdat_preDLI <- expand.grid(
  "ATG" = levels(dat_long$ATG),
  "CMV_PD" = levels(dat_long$CMV_PD),
  "hirisk" = levels(dat_long$hirisk),
  "intSCT2_5" = seq(0, 6, by = 0.05)
)

dat_preds <- predict(
  jm1,
  newdata = newdat_preDLI,
  type = "Marginal",
  idVar = "IDAA",
  returnData = TRUE,
  interval = "confidence"
)

# This plot already answers first q
ggplot(
  data = dat_preds,
  aes(
    x = intSCT2_5,
    y = pred,
    group = interaction(ATG, hirisk),
    col = ATG
  )
) +
  geom_ribbon(
    aes(ymin = low, ymax = upp),
    fill = "gray",
    alpha = 0.5,
    col = NA
  ) +
  geom_line(
    aes(linetype = hirisk),
    size = 1.5,
    alpha = 0.75
  ) +
  # add repel labels
  facet_wrap(facets = ~ CMV_PD) +
  labs(
    x = "Time since alloHCT (months)",
    y = "CD3 cell counts",
    linetype = "ITT",
    col = "Donor type"
  ) +
  scale_y_continuous(
    breaks = log(c(5, 25, 100, 500, 1500)),
    labels = c(5, 25, 100, 500, 1500)
  ) +
  xlim(c(0, 6)) +
  scale_color_manual(
    labels = c("Related", "Unrelated (+ATG)"),
    values = c("brown", "darkblue")
  ) +
  scale_linetype_manual(
    labels = c("No", "Yes (high risk)"),
    values = c("solid", "dotdash")
  )

summary(jm1)
summary(jm2)


# Look at CD8 joint model (which did fit) ---------------------------------


summary(preDLI_CD8_jointModel_both)
summary(preDLI_CD3_jointModel_corr)
preDLI_CD8_jointModel_both$coefficients$gammas
preDLI_CD3_jointModel_corr$coefficients$gammas

# Check marginal CD8 fit


dat_preds <- predict(
  preDLI_CD8_jointModel_both,
  newdata = newdat_preDLI,
  type = "Marginal",
  idVar = "IDAA",
  returnData = TRUE,
  interval = "confidence"
)

# This plot already answers first q
ggplot(
  data = dat_preds,
  aes(
    x = intSCT2_5,
    y = pred,
    group = interaction(ATG, hirisk),
    col = ATG
  )
) +
  geom_ribbon(
    aes(ymin = low, ymax = upp),
    fill = "gray",
    alpha = 0.5,
    col = NA
  ) +
  geom_line(
    aes(linetype = hirisk),
    size = 1.5,
    alpha = 0.75
  ) +
  # add repel labels
  facet_wrap(facets = ~ CMV_PD) +
  labs(
    x = "Time since alloHCT (months)",
    y = "CD3 cell counts",
    linetype = "ITT",
    col = "Donor type"
  ) +
  scale_y_continuous(
    breaks = log(c(5, 25, 100, 500, 1500)),
    labels = c(5, 25, 100, 500, 1500)
  ) +
  xlim(c(0, 6)) +
  scale_color_manual(
    labels = c("Related", "Unrelated (+ATG)"),
    values = c("brown", "darkblue")
  ) +
  scale_linetype_manual(
    labels = c("No", "Yes (high risk)"),
    values = c("solid", "dotdash")
  )


# Try power simulation for three-way inter in this mixed mod --------------


# Base it on number of subjects and number of measurements pp

#...
