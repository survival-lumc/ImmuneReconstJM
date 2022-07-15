library(reprex)

# Reprex 1: baseline hazards flexible -------------------------------------


library(JM)
library(JMbayes2)

# Suppose events only start after 3 years
ids_3y <- with(pbc2.id, id[which(years >= 3)])
pbc2.id_3y <- subset(pbc2.id, id %in% ids_3y)
pbc2_3y <- subset(pbc2, id %in% ids_3y)

# Fit submodels - second example from JM::jointModel()
lmeFit <- lme(log(serBilir) ~ year * drug, random = ~ year | id, data = pbc2_3y)
coxFit <- coxph(Surv(years, status2) ~ drug, data = pbc2.id_3y, x = TRUE)

# Fit with JM
pbc_JM <- jointModel(
  lmeObject = lmeFit,
  survObject = coxFit,
  method = "spline-PH-aGH",
  timeVar = "year",
  ord = 3L, # quadratic B-splines
  lng.in.kn = 9L # internal knots, so 10 "base_hazard_segments"
)

# Check knots (boundary knots repeated ord = 3L times)
pbc_JM$control$knots

# Fit JMbayes2 with defaults
pbc_JMbayes2_m1 <- jm(
  Mixed_objects = lmeFit,
  Surv_object = coxFit,
  time_var = "year",
  Bsplines_degree = 2L, # quadratic, the default
  base_hazard_segments = 10L # also default
)

# Check knots
pbc_JMbayes2_m1$control$knots

# Fit JMbayes2 with knots as formatted by JM
pbc_JMbayes2_m2 <- jm(
  Mixed_objects = lmeFit,
  Surv_object = coxFit,
  time_var = "year",
  knots = pbc_JM$control$knots,
  Bsplines_degree = 2L
)

# Compare coefficients
cbind.data.frame(
  "JM" = coef(pbc_JM, "Event"),
  "JMbayes2 (own knots)" = c(
    coef(pbc_JMbayes2_m1)$gammas,
    coef(pbc_JMbayes2_m1)$association
  ),
  "JMbayes2 (JM knots)" = c(
    coef(pbc_JMbayes2_m2)$gammas,
    coef(pbc_JMbayes2_m2)$association
  )
)

# Mention in text:
# Bsplines_degree = 3L ; knots at negative timepoints

# Reprex 2: fixed intercepts ----------------------------------------------


# https://github.com/drizopoulos/JMbayes2/blob/f4998f398b863969c125199c090c32c1fd7c70ee/R/jm.R#L131
# https://github.com/drizopoulos/JM/blob/038ff493f615f9d86f3f65a1a76c5d6606dd5359/R/jointModel.R#L335


set.seed(202207151)
library(JM)
library(JMbayes2)

# Fit submodel with random slopes but fixed intercept
lmeFit <- lme(log(serBilir) ~ ns(year, 2) * drug, random = ~ 0 + ns(year, 2) | id, data = pbc2)
coxFit <- coxph(Surv(years, status2) ~ drug, data = pbc2.id, x = TRUE)

# JM fit
pbc_JM <- jointModel(
  lmeObject = lmeFit,
  survObject = coxFit,
  method = "piecewise-PH-aGH",
  timeVar = "year",
  lng.in.kn = 10L
)

# JMbayes2 fit
pbc_JMbayes2 <- jm(
  Mixed_objects = lmeFit,
  Surv_object = coxFit,
  time_var = "year",
  Bsplines_degree = 0L,
  knots = list(c(0, pbc_JM$control$knots, max(pbc2.id$years) + 1))
)


# Compare coefficients
cbind.data.frame(
  "JM" = coef(pbc_JM, "Event"),
  "JMbayes2" = c(
    coef(pbc_JMbayes2)$gammas,
    coef(pbc_JMbayes2)$association
  )
)

# Using just linear random slope
update(
  object = pbc_JMbayes2,
  Mixed_objects = update(lmeFit, fixed. = . ~ year * drug, random = ~ 0 + year | id)
)



# The
library(rstanarm)
options(mc.cores = 3)

pbc_stan <- stan_jm(
  formulaLong = log(serBilir) ~ ns(year, 2) * drug + (0 + ns(year, 2) | id),
  dataLong = pbc2,
  formulaEvent = survival::Surv(years, status2) ~ drug,
  dataEvent = pbc2.id,
  time_var = "year",
  basehaz = "piecewise",
  basehaz_ops = list("knots" = pbc_JM$control$knots),
  chains = 3,
  cores = 3,
  refresh = 50
)
coef(pbc_stan)
pbc_stan$basehaz$knots
packageVersion("JM")
packageVersion("JMbayes2")


reprex({
  set.seed(202207151)
  library(JM)
  library(JMbayes2)

  # Fit submodel with random slopes but fixed intercept
  lmeFit <- lme(log(serBilir) ~ ns(year, 2) * drug, random = ~ 0 + ns(year, 2) | id, data = pbc2)
  coxFit <- coxph(Surv(years, status2) ~ drug, data = pbc2.id, x = TRUE)

  # JM fit
  pbc_JM <- jointModel(
    lmeObject = lmeFit,
    survObject = coxFit,
    method = "piecewise-PH-aGH",
    timeVar = "year",
    lng.in.kn = 10L
  )

  # JMbayes2 fit
  pbc_JMbayes2 <- jm(
    Mixed_objects = lmeFit,
    Surv_object = coxFit,
    time_var = "year",
    Bsplines_degree = 0L,
    knots = list(c(0, pbc_JM$control$knots, max(pbc2.id$years) + 1))
  )

  # Compare coefficients
  cbind.data.frame(
    "JM" = coef(pbc_JM, "Event"),
    "JMbayes2" = c(
      coef(pbc_JMbayes2)$gammas,
      coef(pbc_JMbayes2)$association
    )
  )

  # Using just linear random slope
  update(
    object = pbc_JMbayes2,
    Mixed_objects = update(lmeFit, fixed. = . ~ year * drug, random = ~ 0 + year | id)
  )
})
