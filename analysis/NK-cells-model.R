##*********************##
## Modelling NK counts ##
##*********************##

full_dat <- readRDS("data/merged-data.rds")

JM_dat <- full_dat[, c(
  "IDAA",
  "NK_abs", # NK counts
  "intSCT2_2", # Time of measurement (-1 day if at endpoint)
  "EFS3_s", # 1 if REL/ non-relapse fail, zero other wise
  "endpoint2", # Time to endpoint
  "relation", # Donor relation
  "VCMVPAT_pre",
  "PATSEX",
  "age_SCT"
)]

# Limit to first 6 months..
JM_dat[endpoint2 >= 6, ':=' (
  endpoint2 = 6,
  EFS3_s = 0
)]

# No NA measurements, and no long. measures after endpoint
JM_dat <- JM_dat[!is.na(NK_abs) & intSCT2_2 < endpoint2, ]
JM_dat[NK_abs == 0, NK_abs := 0.5] # For undetectables
JM_dat[, NK_abs_log := log(NK_abs)]
JM_dat <- JM_dat[, .SD[.N > 1], by = "IDAA"] # At least two measurements per pat
JM_dat[, IDAA := factor(IDAA)]

# Pick a few patients to visualise some results
set.seed(85)
ref_pats <- sample(JM_dat$IDAA, size = 8, replace = FALSE)

# Plot measurements
JM_dat %>%
  melt.data.table(
    id.vars = c("IDAA", "intSCT2_2"),
    measure.vars = c("NK_abs", "NK_abs_log"),
    variable.name = "scale",
    value.name = "NK_counts"
  ) %>%
  ggplot(aes(intSCT2_2 , NK_counts, col = IDAA)) +
  geom_point(alpha = 0.5) +
  geom_line(aes(group = IDAA), alpha = 0.25) +
  facet_wrap(. ~ scale, scales = "free_y") +
  labs(x = "Time since alloHCT (months)", y = "NK cell counts") +
  theme(legend.position = "none")


# Continuous time (quadratic) ---------------------------------------------


# Estimate independent random effects
lmeFit <- lme(
  NK_abs_log ~ (intSCT2_2 + I(intSCT2_2^2)) * VCMVPAT_pre,
  random = list(IDAA = pdDiag(form = ~ intSCT2_2 + I(intSCT2_2^2))), #~ intSCT2_2 + I(intSCT2_2^2) | IDAA,
  control = lmeControl(opt = 'optim'),
  data = JM_dat
)

summary(lmeFit)

JM_dat_wide <- dcast(
  JM_dat,
  formula = IDAA + VCMVPAT_pre + age_SCT + relation + EFS3_s + endpoint2 ~ .,
  fun = length,
)

coxFit <- coxph(Surv(endpoint2, EFS3_s) ~ VCMVPAT_pre, data = JM_dat_wide, x = TRUE)
summary(coxFit)

jointFit <- jointModel(
  lmeObject = lmeFit,
  survObject = coxFit,
  timeVar = "intSCT2_2", #
  method = "piecewise-PH-aGH", # Piecewise constant baseline hazard
  iter.EM = 200
)

summary(jointFit)
summary(lmeFit)

# Make some individual predictions
preds_indiv_quad <- predict(
  object = jointFit,
  newdata = JM_dat,
  FtTimes = seq(0, 6, by = 0.1),
  type = "Subject",
  idVar = "IDAA",
  interval = "confidence",
  returnData = TRUE
)

preds_indiv_quad[IDAA %in% ref_pats & !is.na(se.fit) & pred > 0] %>%
  ggplot(aes(intSCT2_2, pred)) +
  geom_ribbon(aes(ymin = low, ymax = upp), fill = "lightgray", alpha = 0.75) +
  geom_point(
    data = JM_dat[IDAA %in% ref_pats],
    aes(intSCT2_2, NK_abs_log),
    col = "blue"
  ) +
  geom_line() +
  facet_wrap(~ IDAA, nrow = 2) +
  labs(x = "Time since alloHCT (months)", y = "NK cell counts (log)")



# Marginal prediction
DF <- with(
  JM_dat,
  expand.grid(
    VCMVPAT_pre = levels(VCMVPAT_pre),
    intSCT2_2 = seq(min(intSCT2_2), max(intSCT2_2), len = 100)
  )
)

preds_marg_quad <- predict(jointFit, DF, interval = "confidence", return = TRUE)
preds_marg_quad

preds_marg_quad %>%
  ggplot(aes(intSCT2_2, pred)) +
  geom_point(
    data = JM_dat,
    aes(y = NK_abs_log, col= IDAA), alpha = 0.5
  ) +
  geom_ribbon(aes(ymin = low, ymax = upp), fill = "lightgray", alpha = 0.75) +
  geom_line() +
  facet_wrap(~ VCMVPAT_pre) +
  labs(x = "Time since alloHCT (months)", y = "NK cell counts (log)") +
  theme(legend.position = "none")



# Continous - natural splines ---------------------------------------------


lmeFit_ns <- lme(
  NK_abs_log ~ ns(intSCT2_2, 3) * VCMVPAT_pre,
  random = list(IDAA = pdDiag(form = ~ ns(intSCT2_2, 3))),
  control = lmeControl(opt = 'optim'),
  data = JM_dat
)

summary(lmeFit_ns)

jointFit_ns <- jointModel(
  lmeObject = lmeFit_ns,
  survObject = coxFit,
  timeVar = "intSCT2_2", #
  method = "piecewise-PH-aGH", # Assume weibull basleine hazard, more efficient
  iter.EM = 200
)

# Make some individual predictions
preds_indiv_ns <- predict(
  object = jointFit_ns,
  newdata = JM_dat,
  FtTimes = seq(0, 6, by = 0.1),
  type = "Subject",
  idVar = "IDAA",
  interval = "confidence",
  returnData = TRUE
)


preds_indiv_ns[IDAA %in% ref_pats & !is.na(se.fit)] %>%
  ggplot(aes(intSCT2_2, pred)) +
  geom_ribbon(aes(ymin = low, ymax = upp), fill = "lightgray", alpha = 0.75) +
  geom_point(
    data = JM_dat[IDAA %in% ref_pats],
    aes(intSCT2_2, NK_abs_log),
    col = "blue"
  ) +
  geom_line() +
  facet_wrap(~ IDAA, nrow = 2, scales = "fixed") +
  labs(x = "Time since alloHCT (months)", y = "NK cell counts") +
  scale_y_continuous(
    breaks = log(c(5, 25, 100, 500, 1500)),
    labels = c(5, 25, 100, 500, 1500)
  )


# Marginal
preds_marg_ns <- predict(jointFit_ns, DF, interval = "confidence", return = TRUE)
preds_marg_ns

preds_marg_ns %>%
  ggplot(aes(intSCT2_2, pred)) +
  geom_point(
    data = JM_dat,
    aes(y = NK_abs_log, col= IDAA), alpha = 0.5
  ) +
  geom_ribbon(aes(ymin = low, ymax = upp), fill = "lightgray", alpha = 0.75) +
  geom_line() +
  facet_wrap(~ VCMVPAT_pre) +
  labs(x = "Time since alloHCT (months)", y = "NK cell counts") +
  theme(legend.position = "none") +
  labs(x = "Time since alloHCT (months)", y = "NK cell counts") +
  scale_y_continuous(
    breaks = log(c(5, 25, 100, 500, 1500)),
    labels = c(5, 25, 100, 500, 1500)
  )

# Make a survival prediction
#surv_JM <- survfitJM(jointFit_ns, newdata = JM_dat, idVar = "IDAA")


# Discrete - binned time --------------------------------------------------


JM_dat[, binned_time := cut(
  intSCT2_2, breaks = seq(0, 6, by = 0.75), right = TRUE, include.lowest = TRUE
)]

# Make measure to be mean cell count in period
JM_dat_binned <- JM_dat[, .(
  NK_abs_log = mean(NK_abs_log)
), by = c("binned_time", "IDAA", "VCMVPAT_pre")]


# Run models
lmeFit_binned <- lme(
  NK_abs_log ~ binned_time * VCMVPAT_pre, # skip interaction here
  random = list(IDAA = pdDiag(form = ~ binned_time)), # assume independent
  control = lmeControl(opt = 'optim'),
  data = JM_dat_binned
)

summary(lmeFit_binned)

df_binned <- data.frame(
  "binned_time" = levels(JM_dat_binned$binned_time),
  "PATSEX" = factor("Male", levels = c("Male", "Female"))
)

preds_binned <- predict(lmeFit_binned, JM_dat_binned)

library(gghighlight)

cbind.data.frame(JM_dat_binned, "preds" = preds_binned) %>%
  ggplot(aes(binned_time, preds)) +
  geom_line(aes(group = IDAA)) +
  geom_point(size = 1) +
  gghighlight(IDAA == "5453", use_group_by = FALSE,
              unhighlighted_params = list(alpha = 0.5, linetype = "dashed")) +
  geom_jitter(data = JM_dat_binned,
              aes(binned_time, NK_abs_log, col = IDAA),
              width = 0.2, alpha = 0.5) +
  theme(legend.position = "none") +
  facet_wrap(~ VCMVPAT_pre) +
  labs(x = "Time since alloHCT (months)", y = "NK cell counts (log)")


