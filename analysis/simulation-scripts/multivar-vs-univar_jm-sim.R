# Testing modelling of univariate JM when data-generating model is multivariate

#devtools::install_github("sambrilleman/simjm")
#devtools::install_github('graemeleehickey/joineRML')
library(simjm)
library(ggplot2)
library(JM)
library(JMbayes2)
library(joineRML)
library(survival)

simdat2 <- simjm(
  M = 2,
  n = 500,
  betaEvent_assoc = 0.3,
  balanced = TRUE,
  fixed_trajectory = "linear",
  random_trajectory = "linear",
  seed = 1984
)
attr(simdat2, "params")
table(simdat2$Event$status)

simdat2$Long1 |>
  ggplot(aes(tij, Yij_1, col = factor(id), group = id)) +
  geom_line() +
  theme(legend.position = "none")

simdat2$Long2 |>
  ggplot(aes(tij, Yij_2, col = id, group = id)) +
  geom_line() +
  theme(legend.position = "none")

# Merge longitudinal datasets
dat_long <- merge(simdat2$Long1, simdat2$Long2)
dat_long <- dat_long[order(dat_long$id, dat_long$tij), ]


# Prep models -------------------------------------------------------------


long1 <- lme(
  fixed = Yij_1 ~ tij + Z1 + Z2,
  random = ~ tij | id,
  data = dat_long,
  control = nlme::lmeControl(opt = "optim")
)

long2 <- lme(
  fixed = Yij_2 ~ tij + Z1 + Z2,
  random = ~ tij | id,
  data = dat_long,
  control = nlme::lmeControl(opt = "optim")
)

#
mod_ev <- coxph(Surv(eventtime, status) ~ Z1 + Z2, data = simdat2$Event, x = TRUE)

# True multivariate model
mod_multivar <- jm(
  Surv_object = mod_ev,
  Mixed_objects = list(long1, long2),
  time_var = "tij",
  id_var = "id"
)

mod_multivar

mod_long1 <- jointModel(
  lmeObject = long1,
  survObject = mod_ev,
  timeVar = "tij",
  method = "spline-PH-aGH",
  iter.EM = 200
)

summary(mod_long1)


# joinerml ----------------------------------------------------------------


# Takes aaages
mod_rml <- mjoint(
  formLongFixed = list(
    "long1" = Yij_1 ~ tij + Z1 + Z2,
    "long2" = Yij_2 ~ tij + Z1 + Z2
  ),
  formLongRandom = list(
    "long1" = ~ tij | id,
    "long2" = ~ tij | id
  ),
  formSurv = Surv(eventtime, status) ~ Z1 + Z2,
  data = dat_long,
  survData = simdat2$Event,
  timeVar = "tij"
)

# 25 minutes
mod_rml$comp.time
summary(mod_rml)
covmat <- mod_rml$coefficients$D
sds <- c(1.549, 0.069527, 1.5155, 0.071462)
Dinv <- solve(diag(sds))
Dinv %*% covmat %*% Dinv
