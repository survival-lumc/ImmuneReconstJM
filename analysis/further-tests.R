mod_names <- c(
  paste0("preDLI_JM_value_corr_CD", c(3, 4, 8)),
  paste0("preDLI_JM_value_indep_CD", c(3, 4, 8)),
  paste0("preDLI_JM_both_corr_CD", c(3, 4, 8)),
  paste0("preDLI_JM_both_indep_CD", c(3, 4, 8))
)
names(mod_names) <- mod_names
mod_names

mods <- lapply(mod_names, tar_read_raw)

# Check convergence
sapply(mods, function(mod) mod$convergence)

# Check iterations
sapply(mods, function(mod) mod$iters)

# Check hessian
sapply(mods, function(mod) try(matrixcalc::is.positive.definite(mod$Hessian)))

tar_meta(fields = c(name, time, seconds, warnings)) |>  View()

cbind.data.frame(
  "Converged? (0 is good)" = sapply(mods, function(mod) mod$convergence),

  # Check iterations
  "Iters" = sapply(mods, function(mod) mod$iters),

  # Check hessian
  "hessian good?" = sapply(mods, function(mod) try(matrixcalc::is.positive.definite(mod$Hessian)))
) |>  View()

# Do some likelihood ratios
anova(mods$preDLI_JM_value_indep_CD8, mods$preDLI_JM_both_indep_CD8)
anova(mods$preDLI_JM_value_indep_CD3, mods$preDLI_JM_both_indep_CD3)
anova(mods$preDLI_JM_value_indep_CD4, mods$preDLI_JM_both_indep_CD4)
summary(mods$preDLI_JM_both_indep_CD4)
summary(mods$preDLI_JM_both_indep_CD3)
summary(mods$preDLI_JM_both_indep_CD8)

summary(mods$preDLI_JM_value_indep_CD4)
summary(mods$preDLI_JM_value_indep_CD3)
summary(mods$preDLI_JM_value_indep_CD8)

mods$pre

pbc2.idCR <- crLong(pbc2.id, statusVar = "status",
                    censLevel = "alive", nameStrata = "CR")
coxFit4.pbc <-
  coxph(Surv(years, status2) ~ (drug + age) * CR + strata(CR),
        data = pbc2.idCR, x = TRUE)
coxFit4.pbc |> summary()


mods$preDLI_JM_value_indep_CD8$coefficients$alpha
mods$preDLI_JM_both_indep_CD8$coefficients$alpha
mods$preDLI_JM_both_indep_CD8$coefficients$Dalpha


mods$preDLI_JM_value_indep_CD4$coefficients$alpha
mods$preDLI_JM_both_indep_CD4$coefficients$alpha
mods$preDLI_JM_both_indep_CD4$coefficients$Dalpha


mods$preDLI_JM_value_indep_CD3$coefficients$alpha
mods$preDLI_JM_both_indep_CD3$coefficients$alpha
mods$preDLI_JM_both_indep_CD3$coefficients$Dalpha


# https://raw.githubusercontent.com/VivianePhilipps/marqLevAlgPaper/master/section6-comparison.R



tar_load(c(preDLI_long_corr_CD3, preDLI_cox, derivForm_preDLI))

tar_load(preDLI_JM_both_indep_CD3) #check loglike
preDLI_JM_both_indep_CD4$logLik

mod_qn <- jointModel(
  lmeObject = preDLI_long_corr_CD3,
  survObject = preDLI_cox,
  parameterization = "both",
  CompRisk = TRUE,
  method = "spline-PH-aGH",
  timeVar = "intSCT2_7",
  interFact = list("value" = ~ strata(trans) - 1, "slope" = ~ strata(trans) - 1),
  derivForm = derivForm_preDLI,
  iter.EM = 50,
  iter.qN = 1000,
  lng.in.kn = 2,
  verbose = TRUE
)
# -1747.65

matrixcalc::is.positive.definite(mod_qn$Hessian)
mod_qn$convergence
mod_qn$logLik
mod_qn$control$knots$`trans=1`

summary(mod_qn)



# Tests -------------------------------------------------------------------


run_preDLI_longitudinal(
  cell_line = "CD3_abs_log",
  form_fixed = "ns(intDLI1, 2) * ATG + CMV_PD", # remove cmvpd
  form_random = list(IDAA = pdDiag(~ ns(intDLI1, 2))),
  dat = NMA_postDLI_datasets$long
)

a <- run_preDLI_longitudinal(
  cell_line = "CD4_abs_log",
  form_fixed = "ns(intDLI1, 2) * ATG + CMV_PD", # remove cmvpd
  form_random = ~ ns(intDLI1, 2) | IDAA,
  #form_random = list(IDAA = pdDiag(~ ns(intDLI1, 2))),
  dat = NMA_postDLI_datasets$long
)

b <- run_preDLI_longitudinal(
  cell_line = "CD4_abs_log",
  form_fixed = "ns(intDLI1, 1) * ATG + CMV_PD", # remove cmvpd
  form_random = ~ ns(intDLI1, 1) | IDAA,
  #list(IDAA = pdDiag(~ ns(intDLI1, 1))),
  dat = NMA_postDLI_datasets$long
)
anova(update(a, method = "ML"), update(b, method = "ML"))




# Try pre-DLI corr, diff basehaz ------------------------------------------


# Also try just with value
tar_load(c(preDLI_long_corr_CD4, preDLI_cox, derivForm_preDLI))

basehaz(preDLI_cox, centered = FALSE) |>
  as.data.frame() |>
  ggplot(aes(time, hazard, col = strata)) +
  geom_step(size = 1.5)

mod_qn_value <- jointModel(
  lmeObject = preDLI_long_corr_CD4,
  survObject = preDLI_cox,
  parameterization = "value",
  CompRisk = TRUE,
  method = "spline-PH-aGH",
  timeVar = "intSCT2_7",
  interFact = list("value" = ~ strata(trans) - 1),
  iter.EM = 10L,
  iter.qN = 50L,
  lng.in.kn = 3L,
  verbose = TRUE,
  numeriDeriv = "cd",
  eps.Hes = 1e-04
)

t_seq <- seq(0.05, 6, by = 0.01)
bsf_basehaz <- splineDesign(
  knots = mod_qn_value$control$knots$`trans=1`,
  x = t_seq
)
plot(
  t_seq,
  bsf_basehaz %*% mod_qn_value$coefficients$gammas.bs[1:6]
)

tar_read(preDLI_JM_both_corr_CD4)

summary(mod_qn_value)
matrixcalc::is.positive.definite(mod_qn_value$Hessian, 1e2)
mod_qn_value$convergence
mod_qn_value$logLik
mod_qn_value$control$knots$`trans=1`
mod_qn_value$control$knots$`trans=2`
mod_qn_value$control$knots$`trans=3`

# Model only current value + slope
preDLI_cox_marg <- update(
  object = preDLI_cox,
  . ~ . - ATG.2 - ATG.3 - hirisk.1 - hirisk.2
)

mod_qn_both <- jointModel(
  lmeObject = preDLI_long_corr_CD4,
  survObject = preDLI_cox_marg,
  parameterization = "both",
  CompRisk = TRUE,
  method = "spline-PH-aGH",
  timeVar = "intSCT2_7",
  interFact = list("value" = ~ strata(trans) - 1, "slope" = ~ strata(trans) - 1),
  derivForm = derivForm_preDLI,
  iter.EM = 50,
  iter.qN = 1000,
  lng.in.kn = 3L, # 3 here? and 2 in other? or only 1 in other?
  verbose = TRUE,
  numeriDeriv = "cd", # more precise?
  eps.Hes = 1e-04
)

summary(mod_qn_both)
tar_read(preDLI_JM_both_corr_CD4) |>summary()
matrixcalc::is.positive.definite(mod_qn_both$Hessian)
mod_qn_both$convergence
mod_qn_both$logLik
mod_qn_both$control$knots$`trans=1`

mod_qn_both$coefficients$gammas.bs
mod_qn_both$x$W2 |>  dim()
t_seq <- seq(0.005, 6, by = 0.1)
spline_mat <- splineDesign(
  knots = mod_qn_both$control$knots$`trans=3`,
  x = t_seq
)
plot(t_seq, exp(spline_mat %*% mod_qn_both$coefficients$gammas.bs[13:18]))
tar_load(NMA_preDLI_datasets)
NMA_preDLI_datasets$wide * 3

mod_qn_both$ev <- eigen(mod_qn_both$Hessian, symmetric = TRUE, only.values = TRUE)$values
!all(ev >= -1e-06 * abs(ev[1]))
matrixcalc::is.positive.definite(mod_qn_both$Hessian, tol = 1e-6)
