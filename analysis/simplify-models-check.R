# Simpler JMs check

# To try for pre-DLI models:
# - bin relapse and NRF other together
# - remove three-way interaction with ATG
# - Remove CMV effect?


# Near-singular reffs.. ; try all models with independent reffs pre-DLI
# https://www.theanalysisfactor.com/wacky-hessian-matrix/
# https://cran.r-project.org/web/packages/glmmTMB/vignettes/troubleshooting.html
# Maybe simplify baseline hazard..
#https://www.youtube.com/watch?v=84LpYeyLvmY&t=18s


# Load objects ------------------------------------------------------------

tar_load(
  c(
    NMA_preDLI_datasets,
    preDLI_CD3_long,
    preDLI_cox
  )
)

dat_wide <- NMA_preDLI_datasets$wide
dat_long <- NMA_preDLI_datasets$long

dat_long |>
  ggplot(aes(intSCT2_5, CD3_abs_log, group = IDAA)) +
  geom_line(aes(col = interaction(ATG, hirisk),
                linetype = interaction(ATG, hirisk)),
            alpha = 0.7,
            size = 1.5) +
  theme_minimal() +
  facet_grid(ATG ~ hirisk)

dat_wide[, .(.N), by = c("ATG", "hirisk")]
  #theme(legend.position = "none")

# Simpler Cox model -------------------------------------------------------


cox_simplif <- {

  # Prepare data
  tmat <- trans.comprisk(K = 2, names = c("gvhd", "relapse_nrf"))
  covs <- c("CMV_PD", "hirisk", "ATG", "earlylow_DLI",
            "HLAmismatch_GvH", "relation", "SCT_May2010")

  msdat <- msprep(
    time = c(NA, rep("endpoint7", 2)),
    status = with(
      dat_wide, cbind(
        NA,
        1 * (endpoint7_s == "gvhd"),
        1 * (endpoint7_s %in% c("relapse", "other_nrf"))
      )
    ),
    data = data.frame(dat_wide),
    trans = tmat,
    keep = covs,
    id = "IDAA"
  )

  msdat_expand <- expand.covs(msdat, covs, append = TRUE, longnames = FALSE)

  coxph(
    formula = Surv(time, status) ~
      ATG.1 + hirisk.1 + # GVHD
      ATG.2 + hirisk.2 + strata(trans), # Relapse / NRF
    data = msdat_expand,
    cluster = msdat_expand$IDAA,
    model = TRUE,
    x = TRUE
  )
}


# Try random slopes but no random intercept?
preDLI_CD3_jointModel_value$coefficients$D
melt(
  data.table(preDLI_CD3_jointModel_value$coefficients$D, keep.rownames = TRUE),
  id.vars = "rn"
) |>
  ggplot(aes(rn, variable, fill = value)) +
  geom_tile()


# Simpler longitudinal mods -----------------------------------------------


CD3_simplif <- run_preDLI_longitudinal(
  cell_line = "CD3_abs_log",
  form_fixed = "ns(intSCT2_5, 3)", #* hirisk + CMV_PD",
  #form_random = ~ ns(intSCT2_5, 3) | IDAA,
  form_random = list("IDAA" = pdDiag(~ ns(intSCT2_5, 3))),
  dat = NMA_preDLI_datasets$long
)

CD3_simplif <- run_preDLI_longitudinal(
  cell_line = "CD3_abs_log",
  form_fixed = "ns(intSCT2_5, 3) * hirisk * ATG + CMV_PD", #* hirisk + CMV_PD",
  #form_random = ~ ns(intSCT2_5, 3) | IDAA,
  form_random = list("IDAA" = pdDiag(~ 0 + ns(intSCT2_5, 3))),
  dat = NMA_preDLI_datasets$long
)

CD3_simplif

# Variance of intercept is almost zero..

# This one fits no errors
CD3_simplif <- run_preDLI_longitudinal(
  cell_line = "CD3_abs_log",
  form_fixed = "ns(intSCT2_5, 3)", # hirisk * ATG + CMV_PD",* hirisk + CMV_PD",
  form_random = ~ 0 + ns(intSCT2_5, 3) | IDAA,
  #form_random = list("IDAA" = pdDiag(~ ns(intSCT2_5, 3))),
  dat = NMA_preDLI_datasets$long
)

# Also workss!
CD3_simplif <- run_preDLI_longitudinal(
  cell_line = "CD3_abs_log",
  form_fixed = "ns(intSCT2_5, 3) * hirisk", # hirisk * ATG + CMV_PD",* hirisk + CMV_PD",
  form_random = ~ 0 + ns(intSCT2_5, 3) | IDAA,
  #form_random = list("IDAA" = pdDiag(~ ns(intSCT2_5, 3))),
  dat = NMA_preDLI_datasets$long
)
summary(CD3_simplif)

# now obtain singular vcov mat..
CD3_simplif <- run_preDLI_longitudinal(
  cell_line = "CD3_abs_log",
  form_fixed = "ns(intSCT2_5, 3) * hirisk * ATG", # hirisk * ATG + CMV_PD",* hirisk + CMV_PD",
  form_random = ~ 0 + ns(intSCT2_5, 3) | IDAA,
  #form_random = list("IDAA" = pdDiag(~ ns(intSCT2_5, 3))),
  dat = NMA_preDLI_datasets$long
)
summary(CD3_simplif)

X <- data.frame(model.matrix(~ ns(intSCT2_5, 3) * hirisk * ATG, data = dat_long))
(X$ns.intSCT2_5..3.3)
dev.new()
cor(X[, -1]) |>  corrplot::corrplot()
cor(X[, -1]) |>  View()
# Do without rand intercept..

dat_long[, .SD[1], by = IDAA] |>
  ggplot(aes(intSCT2_5, CD3_abs_log)) +
  geom_point(aes(col = IDAA)) +
  theme(legend.position = "none")

summary(CD3_simplif)


#...
# make functino that creates fyncitonal slope form automatically


# JMs ---------------------------------------------------------------------

system.time({
  jm1 <- jointModel(
    lmeObject = CD3_simplif, #preDLI_CD3_long,
    survObject = preDLI_cox,#cox_simplif,
    CompRisk = TRUE,
    method = "spline-PH-aGH",
    timeVar = "intSCT2_5",
    parameterization = "value",
    #iter.EM = 2000,
    #equal.strata.knots = FALSE, # Different knots per strata
    #lng.in.kn = 1, # baseline haz is cubic spline (ord = 4), try only 4 knots instead of 5
    interFact = list("value" = ~ strata(trans) - 1)
  )
})

summary(jm1)
H <- jm1$Hessian
ev <- eigen(H, symmetric = TRUE, only.values = TRUE)$values
if (!all(ev >= -1e-06 * abs(ev[1])))
  warning("Hessian matrix at convergence is not positive definite.\n")

jm2 <- update(
  jm1,
  parameterization = "both",
  interFact = list("value" = ~ strata(trans) - 1, "slope" = ~ strata(trans) - 1),
  derivForm = list(
    fixed = ~ 0 + dns(intSCT2_5, 3) +
      dns(intSCT2_5, 3):as.numeric(hirisk == "yes") +
      dns(intSCT2_5, 3):as.numeric(ATG == "ALT+ATG") +
      dns(intSCT2_5, 3):as.numeric(hirisk == "yes"):as.numeric(ATG == "ALT+ATG"),
    random = ~ 0 + dns(intSCT2_5, 3),
    indFixed = c(2:4, 8:13, 15:17),
    indRandom = c(2:4)
  )
)



# Try CD4 -----------------------------------------------------------------



CD4_simplif <- run_preDLI_longitudinal(
  cell_line = "CD4_abs_log",
  form_fixed = "ns(intSCT2_7, 3) * hirisk * ATG + CMV_PD", #* hirisk + CMV_PD",
  #form_random = ~ ns(intSCT2_5, 3) | IDAA,
  form_random = list("IDAA" = pdDiag(~ ns(intSCT2_7, 3))), # 0
  dat = NMA_preDLI_datasets$long
)

VarCorr(
  run_preDLI_longitudinal(
    cell_line = "CD4_abs_log",
    form_fixed = "ns(intSCT2_7, 3) * hirisk * ATG + CMV_PD", #* hirisk + CMV_PD",
    #form_random = ~ ns(intSCT2_5, 3) | IDAA,
    form_random = list("IDAA" = pdDiag(~ 0 + ns(intSCT2_7, 3))), # 0
    dat = NMA_preDLI_datasets$long
  )
)

VarCorr(
  run_preDLI_longitudinal(
    cell_line = "CD4_abs_log",
    form_fixed = "ns(intSCT2_7, 3) * hirisk * ATG + CMV_PD", #* hirisk + CMV_PD",
    form_random = ~ 0 + ns(intSCT2_7, 3) | IDAA,
    #form_random = list("IDAA" = pdDiag(~ 0 + ns(intSCT2_7, 3))), # 0
    dat = NMA_preDLI_datasets$long
  )
)



# Investigation why does CD4 not converge ---------------------------------


tar_load(c(preDLI_CD4_long_indep, preDLI_cox, NMA_preDLI_datasets,
           preDLI_CD4_long_corr))

# This converges and no hessian issues..? With INDEP
mod <- jointModel(
  lmeObject = preDLI_CD4_long_corr,
  survObject = preDLI_cox,
  CompRisk = TRUE,
  method = "spline-PH-aGH",
  timeVar = "intSCT2_7",
  parameterization = "value",
  interFact = list("value" = ~ strata(trans) - 1),
  verbose = TRUE,
  iter.EM = 200
)

#plot log-likelihood
summary(mod)

# Issue with second basis function? Standardise log CD4 again?
long_stand <- run_preDLI_longitudinal(
  cell_line = "scale(CD4_abs_log)",
  form_fixed = "ns(intSCT2_7, 3) * hirisk * ATG + CMV_PD", #* hirisk + CMV_PD",
  #form_random = ~ ns(intSCT2_5, 3) | IDAA,
  form_random = list("IDAA" = pdDiag(~ 0 + ns(intSCT2_7, 3))), # 0
  dat = NMA_preDLI_datasets$long
)
mod$logLik
mod_both$loglik
mod_both <- jointModel(
  lmeObject = preDLI_CD4_long_indep,
  survObject = preDLI_cox,
  CompRisk = TRUE,
  method = "spline-PH-aGH",
  timeVar = "intSCT2_7",
  parameterization = "both",
  interFact = list("value" = ~ strata(trans) - 1, "slope" = ~ strata(trans) - 1),
  derivForm = list(
    fixed = ~ 0 + dns(intSCT2_7, 3) +
      dns(intSCT2_7, 3):as.numeric(hirisk == "yes") +
      dns(intSCT2_7, 3):as.numeric(ATG == "ALT+ATG") +
      dns(intSCT2_7, 3):as.numeric(hirisk == "yes"):as.numeric(ATG == "ALT+ATG"),
    random = ~ 0 + dns(intSCT2_7, 3),
    indFixed = c(2:4, 8:13, 15:17),
    indRandom = c(1:3)
  ),
  verbose = TRUE
)

# Standardise CD4s instead of log?

