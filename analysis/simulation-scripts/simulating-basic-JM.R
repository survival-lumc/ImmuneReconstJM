# For JM in R http://www.drizopoulos.com/courses/Int/JMwithR_EMR_2017.pdf page 122
# From simsurv https://cran.r-project.org/web/packages/simsurv/vignettes/simsurv_usage.html#example-4-simulating-under-a-joint-model-for-longitudinal-and-survival-data
# From crowther article Simulating biologically plausible complex survival data
# Flexible parametric joint modelling oflongitudinal and survival data
# https://aosmith.rbind.io/2018/04/23/simulate-simulate-part-2/

library(lme4)
library(nlme)
library(survival)
library(JM)
library(ggplot2)
library(magrittr)
library(data.table)

set.seed(1960)



# Simple linear mixed mod  ------------------------------------------------


# - No covariates
# - No correlation between random effects


# Set sample size and measurement times
n <- 20
measurement_times <- 0:5

# Gather meta data
params <- data.table(
  "id" = 1:n,
  "B0" = 20, # Fixed effect, intercept
  "B1" = 1.5, # Fixed effect, slope
  "b0" = rnorm(n, mean = 0, sd = 1), # Draw random effects for intercept, ~ N(0, 1)
  "b1" = rnorm(n, mean = 0, sd = 0.5), # Draw random effects for slope, ~ N(0, 0.25) (0.25 = 0.5^2, variance)
  "sigma" = 1 # Variance of longitudinal measurements
)
params

# Correlated random effects? Specify var-covar matrix and simulate with MASS::mvrnorm()

# Simulate trajectories - maybe use matrices instead
list_trajectories <- with(
  data = params,
  expr = lapply(measurement_times, function(t) {
    y <- rnorm(n = n, mean = B0 + b0 + (B1 + b1) * t, sd = sqrt(sigma))
    df <- cbind.data.frame("id" = factor(id), "y" = y, "t" = t)
    return(df)
  })
)

df_trajectories <- rbindlist(list_trajectories)

# Fit model
mod_lmm <- nlme::lme(y ~ t, random = ~ t | id, data = df_trajectories)
summary(mod_lmm)
# Or: mod_lmm <- lme4::lmer(y ~ t + (t | id), data = df_trajectories)

cbind.data.frame(df_trajectories, "pred" = predict(mod_lmm, df_trajectories)) %>%
  ggplot(aes(t, y)) +
  geom_point(aes(col = id)) +
  geom_line(aes(y = pred, group = id, col = id), alpha = 0.5) +
  theme(legend.position = "none")


# Simulate basic joint model  ---------------------------------------------


# - One binary covariate (e.g. treatment) ~ Binom(0.5), effect on both long and surv submodels
# - Exponential baseline hazard
# - LMM submodel is thus yi(t) ~ (B0 + b0) + (B1 + b1) * t + beta * X
# - Surv submodel is hi(t) ~ h0(t)exp{gamma * X + alpha * yi(t)} # current value
# - Alpha is association parameter, for now simply linear


# Fix sample size (larger this time), generate treatment
n <- 500
X <- rbinom(n = n, size = 1, prob = 0.5)

# Same as before (but simplify a little), with addition of association and betas for X
params_JM <- data.table(
  "id" = 1:n,
  "X" = X,
  "intercept" = 20 + rnorm(n, mean = 0, sd = 1), # In one go here for intercept + slope
  "slope" = 1.5 + rnorm(n, mean = 0, sd = 0.5),
  "sigma" = 1,
  "alpha" = 0.15, # association
  "beta_X" = 0.5, # Effect of X on longitudinal submodel
  "gamma_X" = 0.25 # Effect of X on survival submodel
)

params_JM


# First - we generate failure times; we can get an expression for survival
# but cannot invert it - so we root solve!
# We can write: log{Hi(t)} = log{H0(t)} + gamma * X + yi(t)
# log{H0(t)} for exponential is log{lambda * t}, since hazard lambda is constant
# Thus Hi(t) = lambda * t * exp{gamma * X + yi(t)}
# We take exp{-Hi(t)} to get survival, and find root of Si(t) - U to get times!


# 1. Make function to invert - will take meta data as parameters
S_min_U <- function(t, lambda, alpha, intercept, slope, X, beta_X, gamma_X, U) {

  cumhaz <- lambda * t * exp(gamma_X * X + alpha * (intercept + slope * t + beta_X * X))
  return(exp(-cumhaz) - U)
}

# 2. Generate times by inverting - use vuniroot for vectors
params_JM[, "time" := rstpm2::vuniroot(
  f = S_min_U,
  U = runif(n = n),
  lambda = 1e-4, # Baseline hazard rate
  alpha = alpha,
  intercept = intercept,
  slope = slope,
  X = X,
  beta_X = beta_X,
  gamma_X = gamma_X,
  interval = cbind(rep(.Machine$double.eps, n), rep(1e10, n)),
  extendInt = "yes"
)$root]


# Check the times
hist(params_JM$time, breaks = 10)

# Add event indicator (no additional censoring)
params_JM[, "delta" := 1]


# 2. Simulate trajectories, and stop them if before time last measurement > event time
measurement_times <- 0:5

list_trajectories <- with(
  data = params_JM,
  expr = lapply(measurement_times, function(t) {
    y <- rnorm(n = n, mean = intercept + slope * t + beta_X * X, sd = sqrt(sigma))
    df <- cbind.data.frame("id" = id, "y" = y, "t" = t, X)
    return(df)
  })
)

df_trajectories <- rbindlist(list_trajectories)

# Keep before event - and keep patients with at least two or more measurements (optional, don't need to do this)
df_trajectories[, "time" := params_JM$time[match(id, params_JM$id)]]
ids_single_meas <- df_trajectories[t < time, .(n_measurements = .N), by = id][n_measurements >= 2][["id"]]
df_trajectories <- df_trajectories[t < time & id %in% ids_single_meas]


# Make plot
df_trajectories %>%
  ggplot(aes(t, y, col = factor(id))) +
  #geom_point(aes(max_t, y), col = "blue", shape = "square") +
  geom_point() +
  geom_line(aes(group = factor(id)), alpha = 0.5) +
  theme(legend.position = "none")


# 3. Create final datasets, and run JM
wide_dat <- with(
  data = params_JM,
  expr = data.table(id, time, delta, X)
)[id %in% ids_single_meas]

setorder(wide_dat, id)


long_dat <- with(
  data = df_trajectories,
  expr = data.table(id, y, "meas_time" = t, time, X)
)

setorder(long_dat, id)


# Check
View(wide_dat)
View(long_dat)
IDAA = pdDiag(~ ns(intSCT2_5_reset, 3))

# Run joint model
lmeFit <- lme(y ~ X + meas_time, random = ~ meas_time | id, data = long_dat)
coxFit <- coxph(Surv(time, delta) ~ X, data = wide_dat, x = TRUE)


jointFit <- jointModel(
  lmeObject = lmeFit,
  survObject = coxFit,
  timeVar = "meas_time",
  method = "spline-PH-aGH" # Spline-based baseline hazard
)

summary(jointFit)

# Different notes:
# - no event times at zero make sure
# - no single longitudinal measurements per person

# Do some slope tests -----------------------------------------------------

wide_dat$X <- factor(wide_dat$X, levels = 0:1, labels = c("no", "yes"))
long_dat$X <- factor(long_dat$X, levels = 0:1, labels = c("no", "yes"))
lmeFit <- lme(
  y ~ ns(meas_time, 2) * X,
  random = list(id = pdDiag(~ ns(meas_time, 2))),
  data = long_dat
)

jointFit <- jointModel(
  lmeObject = lmeFit,
  survObject = coxFit,
  parameterization = "both",
  derivForm =  list(
    fixed = ~ 0 + dns(meas_time, 2) +
      dns(meas_time, 2):as.numeric(X == "yes"),
    indFixed = c(2:3, 5:6),
    random = ~ 0 + dns(meas_time, 2),
    indRandom = c(2:3)
  ),
  timeVar = "meas_time",
  method = "spline-PH-aGH" # Spline-based baseline hazard
)

summary(jointFit)

# Somehow not working? # Make function to get slope
fitted(
  jointFit,
  process = "Longitudinal",
  type = "Slope"
)
object <- jointFit
long_dat[, "currval" := fitted(
  jointFit,
  process = "Longitudinal",
  type = "Subject"
)]
long_dat[, "slope" := ff]

get_int_tangent <- function(x, y, slope) {
  m <- slope * x
  int <- -x * slope + y
  return(int)
}

long_dat[, int_tang := get_int_tangent(
  x = meas_time, y = currval, slope = slope
)]

long_dat[, ':=' (
  start_x = meas_time - 0.2,
  start_y = int_tang + slope * (meas_time - 0.25),
  end_x = meas_time + 0.2,
  end_y = int_tang + slope * (meas_time + 0.25)
)]

long_dat[id %in% sample(id, replace = FALSE, size = 32)] |>
  ggplot(aes(meas_time, y, group = id)) +
  geom_point() +
  geom_line(aes(y = currval)) +
  geom_segment(
    aes(x = start_x, y = start_y, xend = end_x, yend = end_y),
    size = 1,
    arrow = arrow(length = unit(0.03, "npc"))
  ) +
  facet_wrap(~ id)

# Try now with splines
lmeFit <- lme(y ~ X + meas_time, random = ~ list(id = pdDiag(~ meas_time)), data = long_dat)
coxFit <- coxph(Surv(time, delta) ~ X, data = wide_dat, x = TRUE)




# Real data with JM -------------------------------------------------------


#aids <- JM::aids
#aid.id <- JM::aids.id
pbc2 <- pbc2
pbc2.idCR <- crLong(pbc2.id, statusVar = "status",
                    censLevel = "alive", nameStrata = "CR")



lmeFit.pbc <- lme(log(serBilir) ~ drug * (year + I(year^2)),
                  random = ~ year + I(year^2) | id, data = pbc2)

coxFit4.pbc <- coxph(Surv(years, status2) ~ (drug + age) * CR + strata(CR),
        data = pbc2.idCR, x = TRUE)

jointFit11.pbc <- jointModel(
  lmeFit.pbc, coxFit4.pbc,
  timeVar = "year", method = "spline-PH-aGH", CompRisk = TRUE,
  interFact = list(value = ~ CR, data = pbc2.idCR)
)

summary(jointFit11.pbc)

dform <- list(fixed = ~ I(2*year) * drug, indFixed = 3:6,
              random = ~ I(2*year), indRandom = 2:3)

jointFit12.pbc <- jointModel(
  lmeFit.pbc, coxFit4.pbc,
  timeVar = "year", method = "spline-PH-aGH", CompRisk = TRUE,
  parameterization = "both",
  derivForm = dform,
  interFact = list(value = ~ CR, slope = ~ CR, data = pbc2.idCR)
)

summary(jointFit12.pbc)

# Try the tangentsss
object <- jointFit12.pbc
pb2 <- setDT(pbc2)
pbc2[, "currval" := fitted(
  jointFit12.pbc,
  process = "Longitudinal",
  type = "Subject"
)]
pbc2[, "slope" := ff]

get_int_tangent <- function(x, y, slope) {
  m <- slope * x
  int <- -x * slope + y
  return(int)
}

pbc2[, int_tang := get_int_tangent(
  x = year, y = currval, slope = slope
)]

delta_tan <- 0.5
pbc2[, ':=' (
  start_x = year - delta_tan,
  start_y = int_tang + slope * (year - delta_tan),
  end_x = year + delta_tan,
  end_y = int_tang + slope * (year + delta_tan)
)]

pbc2[id %in% sample(id, replace = FALSE, size = 32)] |>
  ggplot(aes(year, log(serBilir), group = id)) +
  geom_point(size = 0.75) +
  geom_line(aes(y = currval), size = 0.5) +
  geom_segment(
    aes(x = start_x, y = start_y, xend = end_x, yend = end_y),
    size = 1,
    arrow = arrow(length = unit(0.03, "npc"))
  ) +
  facet_wrap(~ id)



# Same shabang but with splines -------------------------------------------


pbc2 <- pbc2
pbc2.idCR <- crLong(pbc2.id, statusVar = "status",
                    censLevel = "alive", nameStrata = "CR")



lmeFit.pbc <- lme(
  log(serBilir) ~ ns(year, 3) * drug,
  random = list(id = pdDiag(~ ns(year, 3))),
  data = pbc2
)

coxFit4.pbc <- coxph(Surv(years, status2) ~ (drug + age) * CR + strata(CR),
                     data = pbc2.idCR, x = TRUE)

dform <- list(
  fixed = ~ 0 + dns(year, 3) + dns(year, 3):as.numeric(drug == "D-penicil"),
  indFixed = c(2:4, 6:8),
  random = ~ 0 + dns(year, 3),
  indRandom = 2:4
)

jointFit12.pbc <- jointModel(
  lmeFit.pbc, coxFit4.pbc,
  timeVar = "year", method = "spline-PH-aGH", CompRisk = TRUE,
  parameterization = "both",
  derivForm = dform,
  interFact = list(value = ~ CR, slope = ~ CR, data = pbc2.idCR)
)

summary(jointFit12.pbc)

# Try the tangentsss
object <- jointFit12.pbc
pb2 <- setDT(pbc2)
pbc2[, "currval" := fitted(
  jointFit12.pbc,
  process = "Longitudinal",
  type = "Subject"
)]
pbc2[, "slope" := ff]

get_int_tangent <- function(x, y, slope) {
  m <- slope * x
  int <- -x * slope + y
  return(int)
}

pbc2[, int_tang := get_int_tangent(
  x = year, y = currval, slope = slope
)]

delta_tan <- 0.33
pbc2[, ':=' (
  start_x = year - delta_tan,
  start_y = int_tang + slope * (year - delta_tan),
  end_x = year + delta_tan,
  end_y = int_tang + slope * (year + delta_tan)
)]

set.seed(08908)

pbc2[id %in% sample(id, replace = FALSE, size = 36)] |>
  ggplot(aes(year, log(serBilir), group = id)) +
  geom_point(size = 0.75) +
  geom_line(aes(y = currval), size = 0.5) +
  geom_segment(
    aes(x = start_x, y = start_y, xend = end_x, yend = end_y),
    size = 1,
    arrow = arrow(length = unit(0.03, "npc"))
  ) +
  facet_wrap(~ id)


pbc2[id %in% sample(id, replace = FALSE, size = 9)] |>
  ggplot(aes(year, log(serBilir), group = id)) +
  geom_point(size = 0.75) +
  geom_line(aes(y = currval), size = 0.5) +
  geom_segment(
    aes(x = start_x, y = start_y, xend = end_x, yend = end_y),
    size = 1,
    arrow = arrow(length = unit(0.03, "npc"))
  ) +
  facet_wrap(~ id)



# Test out Gauss-quadrature here by hand ----------------------------------


