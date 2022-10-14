# To-do here:
# - Comparison of marginal and subj specific trajetories JM vs mixed mod
# - Comparison of JM estimates vs time-dep Cox
# - Plotting of marginal cumulative hazards
# - Residuals plotting (marg and subj spec, see book)
# - Checking overall convergence and Hessian positive definiteness
# - Fitting with both lme and lmer (latter give singular fit quicker)
# - Lmer-based repca for justifying random effects
#   VarCorr() and summary(rePCA(...))
# - Maximum likelihood anovas() for value vs value + slope, maybs also
#   for the longitudinal models themselves



# Overall hessian/convergence/iters checks --------------------------------

# For optimisation
# # https://raw.githubusercontent.com/VivianePhilipps/marqLevAlgPaper/master/section6-comparison.R


# Not gvhd sensitivity for now
mod_names <- c(
  paste0("preDLI_JM_value_corr_CD", c(3, 4, 8)),
  paste0("preDLI_JM_both_corr_CD", c(3, 4, 8)),
  paste0("postDLI_JM_corr_CD", c(3, 4, 8))
)
names(mod_names) <- mod_names
mods <- lapply(mod_names, tar_read_raw)
meta_dat <- tar_meta()

# Statistics df
fit_stats <- cbind.data.frame(
  "Converged? (0 is good)" = sapply(mods, function(mod) mod$convergence),
  "# Iters" = sapply(mods, function(mod) mod$iters),
  "Hessian ok?" = sapply(mods, function(mod) try(matrixcalc::is.positive.definite(mod$Hessian))),
  "Time (mins)" = meta_dat[match(mod_names, meta_dat$name), ]$seconds / 60,
  "Loglik" =sapply(mods, function(mod) mod$logLik)
)

fit_stats
anova(mods$preDLI_JM_value_corr_CD3, mods$preDLI_JM_both_corr_CD3)
anova(mods$preDLI_JM_value_corr_CD4, mods$preDLI_JM_both_corr_CD4) # warning to check! - now Fixed
anova(mods$preDLI_JM_value_corr_CD8, mods$preDLI_JM_both_corr_CD8)


# Baseline hazard (might need refitting) ----------------------------------


tar_load(c(preDLI_cox, postDLI_cox, NMA_preDLI_datasets, NMA_postDLI_datasets))

basehaz(preDLI_cox, centered = FALSE) |>
  as.data.frame() |>
  ggplot(aes(time, hazard, col = strata)) +
  geom_step(size = 1.5)

basehaz(postDLI_cox, centered = FALSE) |>
  as.data.frame() |>
  ggplot(aes(time, hazard, col = strata)) +
  geom_step(size = 1.5)

# plot(object, 4)
# lines(with basehaz)


# Use cause-specific model for gvhd (pre and post)
dat_wide <- NMA_preDLI_datasets$wide
dat_long <- NMA_preDLI_datasets$long
GVHD_preDLI <- coxph(
  Surv(endpoint7, endpoint7_s == "relapse") ~ ATG + hirisk,
  data = dat_wide,
  x = TRUE
)

# Plot non-param basehaz
basehaz(GVHD_preDLI, centered = FALSE) |>
  ggplot(aes(time, hazard)) + geom_step()

jm_gvhd <- jointModel(
  lmeObject = preDLI_long_indep_CD4,
  survObject = GVHD_preDLI,
  parameterization = "value",
  method = "spline-PH-aGH",
  timeVar = "intSCT2_7",
  iter.EM = 50,
  iter.qN = 500,
  lng.in.kn = 3L, # Can vary this; with lapply
  verbose = TRUE,
  numeriDeriv = "cd", # more precise
  eps.Hes = 1e-04
)

# Now compare marginals
dat_wide <- dat_wide[order(dat_wide$endpoint7), ]
dat_wide$gvhd_ind <- with(dat_wide, as.numeric(endpoint7_s == "gvhd"))
dat_wide$marg_cumhaz <- mice::nelsonaalen(dat_wide, endpoint7, gvhd_ind)
lines(dat_wide$endpoint7, dat_wide$marg_cumhaz, type = "s")


# Try time-dependent covariate tmerge -------------------------------------





# Residual plots ----------------------------------------------------------



plotResid <- function (x, y, col.loess = "black", ...) {
  plot(x, y, ...)
  lines(lowess(x, y), col = col.loess, lwd = 2)
  abline(h = 0, lty = 3, col = "grey", lwd = 2)
}

# Make wrapper function with ggplot2 facets marginal,subject~ model
# feeds list of models
plot(mods$postDLI_JM_corr_CD8)

# Marginal residuals
resMargY.cd8 <- residuals(postDLI_JM_corr_CD8, process = "Longitudinal", type = "Marginal")
fitMargY.cd8 <- fitted(postDLI_JM_corr_CD8, process = "Longitudinal", type = "Marginal")
plotResid(fitMargY.cd8, resMargY.cd8, xlab = "Fitted Values", ylab = "Marginal Residuals")

# Subj-specific residuals
resSubjY.cd8 <- residuals(postDLI_JM_corr_CD8, process = "Longitudinal", type = "Subject")
fitSubjY.cd8 <- fitted(postDLI_JM_corr_CD8, process = "Longitudinal", type = "Subject")
plotResid(fitSubjY.cd8, resSubjY.cd8, xlab = "Fitted Values", ylab = "Subject Residuals")



# Demo of singularity with rePCA ------------------------------------------


mod <- lmer(
  CD4_abs_log ~ ns(intDLI1, 2) * ATG + CMV_PD + (ns(intDLI1, 2) | IDAA),
  data = NMA_postDLI_datasets$long
)

VarCorr(mod)
summary(rePCA(mod))

# Also using anovas between mixed models with lme..
