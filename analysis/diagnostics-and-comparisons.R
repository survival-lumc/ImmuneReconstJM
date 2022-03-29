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
anova(mods$preDLI_JM_value_corr_CD4, mods$preDLI_JM_both_corr_CD4) # warning to check!
anova(mods$preDLI_JM_value_corr_CD8, mods$preDLI_JM_both_corr_CD8)


# Baseline hazard (might need refitting) ----------------------------------


tar_load(c(preDLI_cox, postDLI_cox))

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
