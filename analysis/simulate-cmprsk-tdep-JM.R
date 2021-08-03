# To fact check some results from our analysis


# Check whether time-dep assoc interaction is recoverable -----------------

# (no comprisks yet)

set.seed(1991)
n <- 5000
n_measurement_pp <- 10

# Intercept makes big difference in root solving
individual_params <- data.table(
  "id" = seq_len(n),
  "X" = rbinom(n = n, size = 1, prob = 0.5),
  "intercept" = -20 + rnorm(n, mean = 0, sd = 0.25), # In one go here for intercept + slope
  "slope" = 1 + rnorm(n, mean = 0, sd = 0.25),
  "t_dli" = runif(n, min = 0, max = 50) # Intermediate event (DLI) time
)

long_params <- list("beta_X" = 0.25, "sigma" = 1)
surv_params <- list("lambda" = 0.1, "gamma_X" = 0.25, "alpha" = 1, "alpha_dli" = -0.1)


#NO CLUE WHY IT DOES NOT WORK for alpha > 1

# Simulate survival times -------------------------------------------------


SminU_noDLI <- function(t, X, intercept, slope, long_params, surv_params, U) {
  y <- intercept + slope * t + long_params$beta_X * X
  cumhaz <- surv_params$lambda * t * exp(
    surv_params$gamma_X * X + surv_params$alpha * y
  )
  return(cumhaz + log(U))
}

SminU_postDLI <- function(t, X, intercept, slope, long_params, surv_params, U) {
  y <- intercept + slope * t + long_params$beta_X * X
  cumhaz <- surv_params$lambda * t * exp(
    surv_params$gamma_X * X + surv_params$alpha * y + 0.25
  )
  return(cumhaz + log(U))
}


# Time to direct death
individual_params[, "t_direct_death" := rstpm2::vuniroot(
  f = SminU_noDLI,
  U = runif(n = n),
  intercept = intercept,
  slope = slope,
  X = X,
  long_params = long_params,
  surv_params = surv_params,
  interval = cbind(rep(.Machine$double.eps, n), rep(1e10, n))
)$root]

# Time to death after dli
individual_params[, "t_dli_death" := rstpm2::vuniroot(
  f = SminU_postDLI,
  U = runif(n = n),
  intercept = intercept,
  slope = slope,
  X = X,
  long_params = long_params,
  surv_params = surv_params,
  interval = cbind(rep(.Machine$double.eps, n), rep(1e10, n)),
  extendInt = "yes"
)$root]

# Prepare for Cox model
individual_params[, ':=' (
  ind_dli = as.numeric(t_dli < t_direct_death),
  time_death = ifelse(t_dli < t_direct_death, t_dli + t_dli_death, t_direct_death),
  ind_death = 1
)]

# Set DLI to large number if not happened
individual_params[, "t_dli" := ifelse(ind_dli == 0, 1e5, t_dli)]
individual_params[, "time_death" := ifelse(time_death <= 0, 0.0001, time_death)]
hist(individual_params$time_death)

# Prepare long data -------------------------------------------------------


# Simulate measurment times, time of DLI, and binary covar
measure_times <- individual_params[, .(
  t_meas = runif(n_measurement_pp, min = 0, max = 100)
), by = "id"][order(id, t_meas)]

dat_long <- merge(measure_times, individual_params, by = "id")

# Simulate trajectories
dat_long[, ':=' (
  y = rnorm(
    n = .N,
    mean = intercept + slope * t_meas + long_params$beta_X * X,
    sd = sqrt(long_params$sigma)),
  id = factor(id)
)]

# Clean trajectories
dat_long <- dat_long[t_meas < time_death]
dat_long <- dat_long[, .SD[.N >= 2], by = "id"]
dat_long[, id := droplevels(id)]

dat_wide <- data.table::dcast(
  data = dat_long,
  formula = id + X + t_dli + time_death + ind_death ~ .,
  fun = length
)
data.table::setnames(dat_wide, old = ".", new = "n_measurements")

df_surv <- tmerge(
  data1 = data.frame(dat_wide),
  data2 = data.frame(dat_wide),
  id = id,
  status = event(time_death, ind_death),
  dli = tdc(t_dli)
)

lmeFit <- lme(y ~ X + t_meas, random = list(id = pdDiag(~ t_meas)), data = dat_long)
coxFit <- coxph(
  Surv(tstart, tstop, status) ~ X + dli,
  data = df_surv,
  model = TRUE,
  x = TRUE,
  cluster = id
)

JM_mod <- jointModel(
  lmeObject = lmeFit,
  survObject = coxFit,
  timeVar = "t_meas",
  method = "spline-PH-aGH",
  #interFact = list(value = ~ dli),
  iter.EM = 200
)

JM_mod |> summary()
