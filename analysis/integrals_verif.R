long_params <- list("beta_X" = 0.25, "sigma" = 1)
surv_params <- list("lambda" = 0.1, "gamma_X" = 0.25,
                    "gamma_dli" = 0.1,
                    "alpha" = 0.25, "alpha_dli" = -0.25)

t <- 1

h <- Vectorize(function(t,
                        X = 1) {

  yi_t <- -20 + 1 * t + long_params$beta_X * X
  surv_params$lambda * exp(surv_params$gamma_X * X + surv_params$alpha * yi_t)
})

H <- Vectorize(function(t,
                        X = 1) {

  y_part <- surv_params$alpha * (-20 + long_params$beta_X * X) + surv_params$gamma_X * X
  surv_params$lambda * exp(y_part) * (exp(surv_params$alpha * 1 * t) - 1) /
    (surv_params$alpha * 1)
})

h(t)
integrate(h, 0, t)$value
H(t)

#--

t <- 50
h <- Vectorize(function(t,
                        X = 1,
                        pre = TRUE) {

  yi_t <- -20 + 1 * t + long_params$beta_X * X
  haz <- if (pre) {
    surv_params$lambda * exp(
      surv_params$gamma_X * X + surv_params$alpha * yi_t
    )
  } else {
    surv_params$lambda * exp(
      surv_params$gamma_X * X + surv_params$alpha * yi_t + surv_params$gamma_dli
    )
  }
  haz
})

H <- Vectorize(function(t,
                        X = 1,
                        Z = 1,
                        t0) {

  y_part <- surv_params$alpha * (-20 + long_params$beta_X * X) + surv_params$gamma_X * X
  lhs <- surv_params$lambda * exp(y_part) * (exp(surv_params$alpha * 1 * t0) - 1)
  rhs <- surv_params$lambda * exp(y_part + surv_params$gamma_dli) * (
    exp(surv_params$alpha * 1 * t) - exp(surv_params$alpha * 1 * t0)
  )
  (lhs + rhs)/ (surv_params$alpha * 1)
})


t0 <- 29
integrate(h, 0, t0, pre = TRUE)$value + integrate(h, t0, t, pre = FALSE)$value
H(t, t0 = t0)



# Try now proper simulation -----------------------------------------------


set.seed(1235)

# Set parameters
long_params <- list(
  "beta_X" = 0.25,
  "sigma" = 1
)
surv_params <- list(
  "lambda" = 0.01,
  "gamma_X" = 0.25,
  "gamma_dli" = 0, #-0.5,
  "alpha" = 0.25,
  "alpha_dli" = -0.2
)

n <- 20000
n_measurement_pp <- 15

# Intercept makes big difference in root solving
individual_params <- data.table(
  "id" = seq_len(n),
  "X" = rbinom(n = n, size = 1, prob = 0.5),
  "intercept" = 5 + rnorm(n, mean = 0, sd = 0.25), # In one go here for intercept + slope,
  "slope" = 1 + rnorm(n, mean = 0, sd = 0.25),
  "t_dli" = runif(n, min = 5, max = 8) # Intermediate event (DLI) time
)

SminU <- function(t, X, intercept, slope, long_params, surv_params, U, t_dli){

  # Pre-work
  y_part <- surv_params$alpha * (intercept + long_params$beta_X * X) + surv_params$gamma_X * X

  # Condition
  cond <- -log(U) < surv_params$lambda * exp(y_part) * (
    exp(surv_params$alpha * slope * t_dli) - 1
  ) / (surv_params$alpha * slope)

  # Cumulative hazards to use
  cumhaz <- if (cond) {
    surv_params$lambda * exp(y_part) * (exp(surv_params$alpha * slope * t) - 1) /
      (surv_params$alpha * slope)
  } else {

    # For interaction
    #alph <- surv_params$alpha
    alph <- surv_params$alpha + surv_params$alpha_dli

    # Careful about lhs alpha..
    lhs <- surv_params$lambda * exp(y_part) *
      (exp(surv_params$alpha * slope * t_dli) - 1) / (surv_params$alpha * slope)

    # rhs y_part is wrong?
    y_part_rhs <- alph * (intercept + long_params$beta_X * X) + surv_params$gamma_X * X
    rhs <- surv_params$lambda * exp(y_part_rhs + surv_params$gamma_dli) * (
      exp(alph * slope * t) - exp(alph * slope * t_dli)
    ) / (alph * slope)

    lhs + rhs
  }

  cumhaz + log(U)
}

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

# Tests
timep <- 0
individual_params[, lp_pre_dli := surv_params$alpha * (
  intercept + slope * timep + long_params$beta_X * X
)]
individual_params[, lp_post_dli := surv_params$gamma_dli + (surv_params$alpha + surv_params$alpha_dli) * (
  intercept + slope * timep + long_params$beta_X * X
)]

individual_params[, diff_lp := lp_post_dli - lp_pre_dli]
mean(individual_params$diff_lp)

individual_params[, "time_death" := ifelse(time_death <= 0, 0.0001, time_death)]
individual_params[, ':=' (
  ind_dli = as.numeric(t_dli < time_death),
  ind_death = 1
)]
individual_params[, "t_dli" := ifelse(ind_dli == 0, 1e5, t_dli)]
table(individual_params$ind_dli)
hist(individual_params$time_death)

measure_times <- individual_params[, .(
  t_meas = seq(0.1, 20, length.out = n_measurement_pp)
), by = "id"][order(id, t_meas)]

dat_long <- merge(measure_times, individual_params, by = "id")

# Simulate trajectories
dat_long[, ':=' (
  y = rnorm(
    n = .N,
    mean = intercept + slope * t_meas + long_params$beta_X * X,
    sd = sqrt(long_params$sigma)
  ),
  y_true = intercept + slope * t_meas + long_params$beta_X * X,
  id = factor(id)
)]

# Clean trajectories
dat_long <- dat_long[t_meas < time_death]
dat_long <- dat_long[, .SD[.N >= 2], by = "id"]
dat_long[, id := droplevels(id)]
dat_long[, .(av = .N), by = list(id, t_meas <= t_dli)][, .(mean(av)), by = t_meas]

# Try t merge here
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

# See https://cran.r-project.org/web/packages/survival/vignettes/timedep.pdf
df2 <- tmerge(
  data1 = df_surv,
  data2 = data.frame(dat_long),
  id = id,
  y = tdc(t_meas, y),
  y_true = tdc(t_meas, y_true)
)

coxph(
  Surv(tstart, tstop, status) ~ X + dli * y_true,
  data = df2,
  model = TRUE,
  x = TRUE,
  cluster = id
) |> coef()

#setDT(df2)
#df2[, y_var := -21, by = id] # Not sure if good practice

#coxFit_tdc <-
coxph(
  Surv(tstart, tstop, status) ~ X + dli + y,
  data = df2,
  model = TRUE,
  x = TRUE,
  cluster = id
)

coxph(
  Surv(tstart, tstop, status) ~ X + dli,
  data = df2,
  model = TRUE,
  x = TRUE,
  cluster = id
)


coxph(
  Surv(tstart, tstop, status) ~ X + dli + y_true,
  data = df2,
  model = TRUE,
  x = TRUE,
  cluster = id
)

coxph(
  Surv(tstart, tstop, status) ~ X + dli * y_true,
  data = df2,
  model = TRUE,
  x = TRUE,
  cluster = id
)


coxph(
  Surv(tstart, tstop, status) ~ X + y_true + y_true:dli,
  data = df2,
  model = TRUE,
  x = TRUE,
  cluster = id
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
  interFact = list(value = ~ dli), # and without
  iter.EM = 200
)

summary(JM_mod)

# Cannot estimate main effect DLI since it is in the alpha change too


# Bayesian way ------------------------------------------------------------



summary(JM_mod2)

# With jmbayes2
JMbayes2_mod <- jm(
  Mixed_objects = lmeFit,
  Surv_object = coxFit,
  time_var = "t_meas",
  functional_forms = list("y" = ~ value(y) * dli),
  data_Surv = df_surv
)

summary(JM_mod)
JMbayes2_mod


# Look at data
hist(individual_params$t_dli[individual_params$t_dli < 1e5], main = "DLI times") # DLI times

# Make better data
ggplot(dat_long, aes(t_meas, y, group = id, col = id)) +
  geom_line(alpha = 0.5) +
  facet_wrap(~ X) +
  theme(legend.position = "none")

dat_long[t_meas < t_dli, .("meano" = mean(y)), by = id][["meano"]] |> mean()
dat_long[t_dli < t_meas, .("meano" = mean(y)), by = id][["meano"]] |> mean()
summary(JM_mod)[["CoefTable-Event"]][1:4, ]

# dli coef is log(HR) for a increase of 0


# Try sim just normal Cox model, interaction time-varyiing covar with a normal variable
# (with Austin method), try interpreting main effect of time-varying covar



# Basic collapsability ----------------------------------------------------


X <- rbinom(n = n, size = 1, prob = 0.5)
Z <- rnorm(n)
t <- rexp(n, rate = 0.01 * exp(
  0.25 * Z + 0.5 * X - 0.25 * X * Z
))
simdat <- cbind.data.frame(t, delta = 1, X, Z)

coxph(Surv(t, delta) ~ Z + X, data = simdat)
coxph(Surv(t, delta) ~ Z * X, data = simdat)
