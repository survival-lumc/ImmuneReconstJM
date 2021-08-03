set.seed(49849)
n <- 2000
n_measurement_pp <- 10

# Intercept makes big difference in root solving
individual_params <- data.table(
  "id" = seq_len(n),
  "X" = rbinom(n = n, size = 1, prob = 0.5),
  "intercept" = -20 + rnorm(n, mean = 0, sd = 0.25), # In one go here for intercept + slope
  "slope" = 1 + rnorm(n, mean = 0, sd = 0.25),
  "t_dli" = runif(n, min = 0, max = 20) # Intermediate event (DLI) time
)



SminU <- function(t, X, intercept, slope, long_params, surv_params, U,
                  t_dli) {
  yt <- intercept + slope * t + long_params$beta_X * X
  yt_dli <- intercept + slope * t_dli + long_params$beta_X * X

  cond <- -log(U) < surv_params$lambda *
    exp(surv_params$gamma_X * X + surv_params$alpha * yt_dli) * t_dli

  haz <- Vectorize(function(s) {

    res <- if (cond) {
      surv_params$lambda * exp(surv_params$gamma_X * X + surv_params$alpha * yt)
    } else {
      surv_params$lambda * exp(surv_params$gamma_X * X +
                                 surv_params$alpha * yt +
                                 surv_params$gamma_dli)
    }
    res
  })

  #cat(cond)
  cumhaz <- if (cond) {
    surv_params$lambda * exp(surv_params$gamma_X * X + surv_params$alpha * yt) * t
  } else {
    surv_params$lambda * exp(surv_params$gamma_X * X +
                               surv_params$alpha * yt +
                               #(surv_params$alpha + surv_params$alpha_dli) * yt +
                               surv_params$gamma_dli) * t
    # surv_params$lambda * exp(surv_params$gamma_X * X + surv_params$alpha * yt) *
    #   (t_dli + exp(surv_params$gamma_dli) * t - exp(surv_params$gamma_dli) * t_dli)
  }

  cumhaz + log(U)
  # c(
  #   cond,
  #   cumhaz,
  #   integrate(haz, lower = 0, upper = t)$value
  # )
}

lapply(seq_len(25), function(i) {
  row <- as.data.frame(individual_params[i, ])
  SminU(
    t = 10,
    X = row[["X"]],
    intercept = row[["intercept"]],
    slope = row[["slope"]],
    long_params = long_params,
    surv_params = surv_params,
    t_dli = row[["t_dli"]],
    U = runif(1)
  )
})


individual_params[, time_death := apply(.SD, 1, function(row) {

  uniroot(
    f = SminU,
    interval = c(0, 1e3),
    extendInt = "yes",
    X = row["X"],
    intercept = row["intercept"],
    slope = row["slope"],
    long_params = long_params,
    surv_params = surv_params,
    t_dli = row["t_dli"],
    U = runif(1)
  )$root

})]

individual_params[, "time_death" := ifelse(time_death <= 0, 0.0001, time_death)]
individual_params[, ':=' (
  ind_dli = as.numeric(t_dli < time_death),
  ind_death = 1
)]
individual_params[, "t_dli" := ifelse(ind_dli == 0, 1e5, t_dli)]
table(individual_params$ind_dli)

# Intermediate test
df_surv <- tmerge(
  data1 = data.frame(individual_params),
  data2 = data.frame(individual_params),
  id = id,
  status = event(time_death, ind_death),
  dli = tdc(t_dli)
)

coxph(
  Surv(tstart, tstop, status) ~ X + dli,
  data = df_surv,
  model = TRUE,
  x = TRUE,
  cluster = id
)


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
coxFit$coefficients
JM_mod <- jointModel(
  lmeObject = lmeFit,
  survObject = coxFit,
  timeVar = "t_meas",
  method = "spline-PH-aGH",
  #interFact = list(value = ~ dli),
  iter.EM = 200
)

summary(JM_mod)


# Functions-based ---------------------------------------------------------

# Set-up parameters
long_params <- list(
  "beta_X" = 0.25,
  "sigma" = 1,
  "B0" = -20,
  "B1" = 1,
  "var_b0" = 0.5,
  "var_b1" = 0.5
)
surv_params <- list(
  "lambda" = 0.1,
  "gamma_X" = 0.25,
  "gamma_dli" = 0.25,
  "alpha" = 0.25,
  "alpha_dli" = -0.25
)

generate_dats <- function(n = 2000,
                          n_measures = 10,
                          long_params,
                          surv_params) {

  # Individual parameters
  individual_params <- data.table(
    "id" = seq_len(n),
    "X" = rbinom(n = n, size = 1, prob = 0.5),
    "intercept" = long_params$B0 + rnorm(n, mean = 0, sd = sqrt(long_params$var_b0)),
    "slope" = long_params$B1 + rnorm(n, mean = 0, sd = sqrt(long_params$var_b1)),
    "t_dli" = runif(n, min = 0, max = 20)
  )

}

dat_wide <- ...
dat_long <- ...
