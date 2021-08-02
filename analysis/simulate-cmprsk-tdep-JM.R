# To fact check some results from our analysis

set.seed(1980)
n <- 200
n_measurement_pp <- 6

params <- data.table(
  "id" = seq_len(n),
  "B0" = 20, # Fixed effect, intercept
  "B1" = 1.5, # Fixed effect, slope
  "b0" = rnorm(n, mean = 0, sd = 1), # Draw random effects for intercept, ~ N(0, 1)
  "b1" = rnorm(n, mean = 0, sd = 0.5), # Draw random effects for slope, ~ N(0, 0.25) (0.25 = 0.5^2, variance)
  "sigma" = 1 # Variance of longitudinal measurements
)

# Simulate measurment times, time of DLI, and binary covar
measure_times <- params[, .(
  t_meas = runif(n_measurement_pp, min = 0, max = 12),
  t_dli = runif(1, min = 5, max = 7),
  X = rbinom(1, 1, prob = 0.5)
), by = "id"][order(id, t_meas)]

dat_long <- merge(measure_times, params, by = "id")

# Simulate trajectories
dat_long[, ':=' (
  y = rnorm(n = .N, mean = B0 + b0 + (B1 + b1) * t_meas, sd = sqrt(sigma)),
  id = factor(id)
)]

# model on log scale

