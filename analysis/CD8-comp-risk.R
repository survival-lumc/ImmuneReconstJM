##************************************##
## CD8 modelling with competing risks ##
##************************************##

full_dat <- readRDS("data/merged-data.rds")

JM_dat <- full_dat[, c(
  "IDAA",
  "CD8_abs", # NK counts
  "intSCT2_2", # Time of measurement (-1 day if at endpoint)
  "endpoint3_s", # Competing risks
  "endpoint2", # Time to endpoint
  "endpoint_specify",
  "relation", # Donor relation
  "TCDmethod",
  "SCTyear" # edit to pre-2010 for DLI
)]


# Look at trajectories AFTER DLI ------------------------------------------

JM_dat <- JM_dat[!is.na(CD8_abs) & intSCT2_2 <= 24]
JM_dat[CD8_abs == 0, CD8_abs := 0.5] # For undetectables
JM_dat[, CD8_abs_log := log(CD8_abs)]
JM_dat <- JM_dat[, .SD[.N > 1], by = "IDAA"]

JM_dat_dli <- JM_dat[grepl(x = endpoint_specify, pattern = "DLI")] #DLI
JM_dat_dli <- JM_dat_dli[, post_DLI := factor(
  as.numeric(intSCT2_2 > endpoint2),
  levels = c(0, 1),
  labels = c("pre_DLI", "post_DLI")
)]


set.seed(86)
ref_pats <- sample(JM_dat_dli$IDAA, size = 8, replace = FALSE)


JM_dat_dli[IDAA %in% ref_pats] %>%
  ggplot(aes(intSCT2_2, CD8_abs_log)) +
  geom_point(
    aes(fill = post_DLI),
    size = 3.5, pch = 21,
    alpha = 0.8,
    col = "#359fda",
    #fill = colorspace::lighten("#359fda", amount = 0.3)
  ) +
  geom_path(alpha = 0.25, linetype = "dotted") +
  #geom_line(alpha = 0.25, aes(linetype = post_DLI))
  geom_vline(aes(xintercept = endpoint2), linetype = "dashed") +
  facet_wrap(~ IDAA, nrow = 2, scales = "free_x") +
  labs(x = "Time since alloHCT (months)", y = "CD8 cell counts") +
  scale_y_continuous(
    breaks = log(c(5, 25, 100, 500, 1500)),
    labels = c(5, 25, 100, 500, 1500)
  ) +
  #xlim(c(0, 24)) +
  theme_minimal()




# Resr --------------------------------------------------------------------



# Limit to first 6 months..
JM_dat[endpoint2 >= 6, ':=' (
  endpoint2 = 6,
  endpoint3_s = "censored"
)]

# No NA measurements, and no long. measures after endpoint
JM_dat <- JM_dat[!is.na(CD8_abs) & intSCT2_2 < endpoint2, ]
JM_dat[CD8_abs == 0, CD8_abs := 0.5] # For undetectables
JM_dat[, CD8_abs_log := log(CD8_abs)]
JM_dat <- JM_dat[, .SD[.N > 1], by = "IDAA"] # At least two measurements per pat
JM_dat[, IDAA := factor(IDAA)]
JM_dat[, SCTyear_2010 := factor(
  ifelse(SCTyear < 2010, "before", "after"),
  levels = c("before", "after")
)]

set.seed(86)
ref_pats <- sample(JM_dat$IDAA, size = 8, replace = FALSE)

# Plot measurements
JM_dat[IDAA %in% ref_pats] %>%
  ggplot(aes(intSCT2_2, CD8_abs_log)) +
  geom_point(
    data = JM_dat[IDAA %in% ref_pats],
    aes(intSCT2_2, CD8_abs_log),
    size = 3, pch = 21, alpha = 0.8,
    col = "#359fda",
    fill = colorspace::lighten("#359fda", amount = 0.3)
  ) +
  geom_line(alpha = 0.5, linetype = "dotted") +
  #geom_smooth(se = FALSE) +
  facet_wrap(~ IDAA, nrow = 2, scales = "fixed") +
  labs(x = "Time since alloHCT (months)", y = "CD8 cell counts") +
  scale_y_continuous(
    breaks = log(c(5, 25, 100, 500, 1500)),
    labels = c(5, 25, 100, 500, 1500)
  ) +
  theme_minimal()


# Use crLong.. - see analysis log too
JM_dat[, delta := as.numeric(endpoint3_s) - 1]

JM_dat_wide <- dcast(
  JM_dat,
  formula = IDAA + SCTyear_2010 + TCDmethod + delta + endpoint2 ~ .,
  fun = length
)

# Prepare for mstate
tmat <- trans.comprisk(2, names = c("rel", "non_rel_failure"))
tmat

JM_dat_wide[, ':=' (
  ev_rel = as.numeric(delta == 1),
  ev_nrf = as.numeric(delta == 2)
)]

covs <- c("TCDmethod", "SCTyear_2010")
JM_msdat <- msprep(
  time = c(NA, "endpoint2", "endpoint2"),
  status = c(NA, "ev_rel", "ev_nrf"),
  data = data.frame(JM_dat_wide),
  trans = tmat,
  keep = covs,
  id = "IDAA"
)

JM_msdat_expand <- expand.covs(JM_msdat, covs, append = TRUE, longnames = FALSE)

JM_data_crprep <- crLong(
   data = JM_dat_wide,
   statusVar = "delta",
   censLevel = 0,
   nameStatus = "status",
   nameStrata = "trans"
)

# Does not work
coxFit <- coxph(Surv(Tstart, Tstop, status) ~
                  TCDmethod1.1 + TCDmethod2.1 + SCTyear_2010.1 +
                  TCDmethod1.2 + TCDmethod2.2 + SCTyear_2010.2 +
                  strata(trans),
                data = JM_msdat_expand,
                x = TRUE,
                model = TRUE)

# Include inter
coxFit2 <- coxph(Surv(endpoint2, status) ~ (TCDmethod + SCTyear_2010) * trans +
                  strata(trans),
                data = JM_data_crprep,
                x = TRUE,
                model = TRUE)

lmeFit <- lme(
  CD8_abs_log ~ ns(intSCT2_2, 3),
  random = list(IDAA = pdDiag(form = ~ ns(intSCT2_2, 3))),
  control = lmeControl(opt = 'optim'),
  data = JM_dat
)

jointFit <- jointModel(
  lmeObject = lmeFit,
  survObject = coxFit2,
  timeVar = "intSCT2_2", #
  method = "spline-PH-aGH",# or piecewise
  CompRisk = TRUE,
  interFact = list(value = ~ trans, data = JM_msdat_expand),
  iter.EM = 200
)

summary(jointFit)

new_dat <- expand.grid(
  "delta" = c(0, 1, 2),
  "intSCT2_2" = seq(0, 6, by = 0.1)
)

preds_marg <- predict(
  object = jointFit,
  newdata = new_dat,
  #FtTimes = seq(0, 6, by = 0.1),
  #type = "Marginal",
  #idVar = "IDAA",
  interval = "confidence",
  returnData = TRUE
)

preds_marg %>%
  ggplot(aes(intSCT2_2, pred)) +
  geom_ribbon(aes(ymin = low, ymax = upp), fill = "lightgray", alpha = 0.75) +
  geom_line(alpha = 0.75) +
  #geom_smooth(se = FALSE) +
  facet_wrap(~ delta, ncol = 3) +
  labs(x = "Time since alloHCT (months)", y = "CD8 cell counts") +
  scale_y_continuous(
    breaks = log(c(5, 25, 100, 500, 1500)),
    labels = c(5, 25, 100, 500, 1500)
  ) +
  theme_minimal()



?predict.jointModel


# Or just JM bayes


# With JMbayes2 -----------------------------------------------------------



library(JMbayes2)

coxFit3 <- coxph(Surv(endpoint2, status) ~ (TCDmethod + SCTyear_2010) * strata(trans),
                 data = JM_data_crprep,
                 x = TRUE,
                 model = TRUE)


jFit_CR <- jm(coxFit3, lmeFit, time_var = "intSCT2_2",
             n_iter = 5000L, n_burnin = 1000L)

summary(jFit_CR)
ggtraceplot(jFit_CR)

predict(
  jFit_CR
)
