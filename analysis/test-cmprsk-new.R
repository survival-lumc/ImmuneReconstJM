tar_load(dat_merged)

# Prepare minimal dataset for analysis
dat <- dat_merged[!is.na(CD8_abs_log), c(
  "IDAA",
  "CD8_abs_log",
  "intSCT2_2",
  "endpoint2",
  "endpoint4_s",
  "endpoint4_s_dli",
  "endpoint_specify",
  "TCD2",
  "hirisk",
  "SCTyear_2010",
  "VCMVPAT_pre"
)]

# Limit to first 12 months after SCT, keep only cell measures pre-endpoint
dat[endpoint2 >= 12, ':=' (endpoint2 = 12, endpoint4_s_dli = "censored")]
dat <- dat[intSCT2_2 < endpoint2] # measures pre-endpoint
dat <- dat[, .SD[.N > 1], by = "IDAA"] # At least 2 measurements per subject

# Collapse RD: ALT and UD: ALT
dat[, TCD2 := factor(
  ifelse(TCD2 == "UD: ALT + ATG", "ATG", "noATG"),
  levels = c("ATG", "noATG")
)]

# No contrast available for ATG.. for now exclude the 6 pats with this failure
dat[endpoint4_s_dli == "non-relapse failure: GvHD", .SD[1], by = "IDAA"][["TCD2"]] %>%
  table()

dat <- dat[endpoint4_s_dli != "non-relapse failure: GvHD"]

# Get wide data
JM_dat_wide <- data.table::dcast(
  dat,
  formula = IDAA + SCTyear_2010 + hirisk + TCD2 + VCMVPAT_pre + endpoint4_s_dli + endpoint2 ~ .,
  fun = length
)

# Prep long data for competing risks
JM_dat_long <- JM::crLong(
  data = JM_dat_wide,
  statusVar = "endpoint4_s_dli",
  censLevel = "censored",
  nameStatus = "status",
  nameStrata = "trans"
)

JM_dat_long[trans == "non-relapse failure: other", trans := "NRF_other"]
JM_dat_long[, trans := droplevels(trans)]
levels(JM_dat_long$trans) <- c("censored_DLI", "relapse", "NRF_other")

# Fit longitudinal part
lmeFit <- nlme::lme(
  CD8_abs_log ~ ns(intSCT2_2, 3) * TCD2 + VCMVPAT_pre,
  random = list(IDAA = pdDiag(form = ~ ns(intSCT2_2, 3))),
  #random = ~ ns(intSCT2_2, 3) | IDAA,
  control = lmeControl(opt = 'optim'),
  data = dat
)

summary(lmeFit)

# JM_dat_long <- transform(
#   JM_dat_long,
#   hirisk_rel = interaction((trans == "relapse"), hirisk, drop = TRUE, sep = ":")
# )

# Try expanding model matrix
modmat_JM_long <- data.table::data.table(
  stats::model.matrix(
    ~ (SCTyear_2010 + hirisk + TCD2) * trans + trans,
    data = JM_dat_long
  )
)

# Keep only relevant contrasts
modmat_sub <- modmat_JM_long[, c(
  "transrelapse",
  "transNRF_other",
  "SCTyear_2010pre_2010",
  "hiriskyes",
  "TCD2noATG:transrelapse",
  "hiriskyes:transrelapse",
  "TCD2noATG:transNRF_other"
)]

colnames(modmat_sub) <- gsub(colnames(modmat_sub), pattern = ":", replacement = "_")

JM_long_new <- cbind(
  JM_dat_long[, c("IDAA", "endpoint2", "status", "trans")],
  modmat_sub
)

JM_long_new[trans != "censored_DLI", ':=' (
  hiriskyes = 0,
  SCTyear_2010pre_2010 = 0
)]

form <- stats::reformulate(
  termlabels = c(
    colnames(modmat_sub),
    "strata(trans)",
    "cluster(IDAA)"
  ),
  response = "Surv(endpoint2, status)"
)

# Fit cox part
coxFit <- coxph(
  formula = form,
  data = JM_long_new,
  x = TRUE,
  model = TRUE
)

broom::tidy(coxFit)

# Check DLI model and relapse
broom::tidy(
  coxph(
    formula = Surv(endpoint2, endpoint4_s_dli == "relapse") ~
      TCD2 + hirisk,
    data = JM_dat_wide
  )
)

# Try same with mstate!!
JMfit <- JM::jointModel(
  lmeObject = lmeFit,
  survObject = coxFit,
  timeVar = "intSCT2_2",
  method = "spline-PH-aGH",# or piecewise
  CompRisk = TRUE,
  interFact = list("value" = ~ strata(trans)), # since strata(trans) is
  # in the model frame of the surv object, not trans!!
  iter.EM = 200
)

round(summary(JMfit)$`CoefTable-Long`, 3)
round(summary(lmeFit)$tTable, 3)
round(summary(JMfit)$`CoefTable-Event`, 3)
round(summary(coxFit)$coefficients, 3)

# Make a plot..
newdat <- expand.grid(
  "TCD2" = levels(dat$TCD2),
  "VCMVPAT_pre" = levels(dat$VCMVPAT_pre),
  "intSCT2_2" = seq(0.1, 12, by = 0.1)
)

preds <- JM::predict.jointModel(
  JMfit,
  newdata = newdat,
  type = "Marginal",
  idVar = "IDAA",
  returnData = TRUE,
  interval = "confidence"
)

preds %>%
  ggplot(aes(intSCT2_2, pred, ymin = low, ymax = upp,
             group = TCD2, col = TCD2)) +
  geom_ribbon(fill = "gray", alpha = 0.5, col = NA) +
  geom_line(size = 1.5, alpha = 0.75) +
  facet_wrap(~ VCMVPAT_pre) +
  labs(x = "Time since alloHCT (months)", y = "CD8 cell counts") +
  scale_y_continuous(
    breaks = log(c(5, 25, 100, 500, 1500)),
    labels = c(5, 25, 100, 500, 1500)
  ) +
  theme_bw()

preds %>%
  ggplot(aes(intSCT2_2, exp(pred), ymin = low, ymax = upp,
             group = TCD2, col = TCD2)) +
  #geom_ribbon(fill = "gray", alpha = 0.5, col = NA) +
  geom_line(size = 1.5, alpha = 0.75) +
  facet_wrap(~ VCMVPAT_pre) +
  labs(x = "Time since alloHCT (months)", y = "CD8 cell counts") +
  # scale_y_continuous(
  #   breaks = log(c(5, 25, 100, 500, 1500)),
  #   labels = c(5, 25, 100, 500, 1500)
  # ) +
  theme_bw()



set.seed(80)
ref_pats <- sample(dat$IDAA, size = 8, replace = FALSE)


dat_merged[, post_endpoint := factor(
  ifelse(intSCT2_2 > endpoint2, 1, 0),
  levels = c(0, 1),
  labels = c("Pre-endpoint", "Post-endpoint")
)]

dat_merged[post_endpoint == "Post-endpoint", ev_label := ifelse(
  intSCT2_2 == min(intSCT2_2), endpoint4_s_dli, NA_character_
), by = "IDAA"]

dat_merged[IDAA %in% ref_pats] %>%
  ggplot(aes(intSCT2_2, CD8_abs_log)) +
  geom_point(
    aes(fill = post_endpoint),
    size = 3.5, pch = 21,
    alpha = 0.8,
    col = "#359fda",
    #fill = colorspace::lighten("#359fda", amount = 0.3)
  ) +
  geom_text(aes(label = ev_label, y = 1), hjust = 0, size = 3) +
  geom_path(alpha = 0.25, linetype = "dotted") +
  geom_vline(aes(xintercept = endpoint2), linetype = "dashed") +
  facet_wrap(~ IDAA, nrow = 2, scales = "free_x") +
  labs(x = "Time since alloHCT (months)", y = "CD8 cell counts") +
  scale_y_continuous(
    breaks = log(c(5, 25, 100, 500, 1500)),
    labels = c(5, 25, 100, 500, 1500)
  ) +
  xlim(c(0, 12)) +
  theme_minimal()

set.seed(1822)
sample_cens <- sample(dat[dat$endpoint4_s_dli == "censored"]$IDAA,
                      size = 8, replace = FALSE)

dat_merged[IDAA %in% sample_cens] %>%
  ggplot(aes(intSCT2_2, CD8_abs_log)) +
  geom_point(
    aes(fill = post_endpoint),
    size = 3.5, pch = 21,
    alpha = 0.8,
    col = "#359fda",
    #fill = colorspace::lighten("#359fda", amount = 0.3)
  ) +
  geom_text(aes(label = ev_label, y = 1), hjust = 0, size = 3) +
  geom_path(alpha = 0.25, linetype = "dotted") +
  geom_vline(aes(xintercept = endpoint2), linetype = "dashed") +
  facet_wrap(~ IDAA, nrow = 2, scales = "free_x") +
  labs(x = "Time since alloHCT (months)", y = "CD8 cell counts") +
  scale_y_continuous(
    breaks = log(c(5, 25, 100, 500, 1500)),
    labels = c(5, 25, 100, 500, 1500)
  ) +
  xlim(c(0, 12)) +
  theme_minimal()


# JMBayes -----------------------------------------------------------------


library(JMbayes)

tmat <- trans.comprisk(3, names = c("rel", "NRF_other", "cens_DLI"))
tmat

JM_dat_wide[, ':=' (
  ev_rel = as.numeric(endpoint4_s_dli == "relapse"),
  ev_nrf = as.numeric(endpoint4_s_dli == "non-relapse failure: other"),
  ev_dli = as.numeric(endpoint4_s_dli == "censored_DLI")
)]

covs <- c("TCD2", "SCTyear_2010", "hirisk")
JM_msdat <- msprep(
  time = c(NA, "endpoint2", "endpoint2", "endpoint2"),
  status = c(NA, "ev_rel", "ev_nrf", "ev_dli"),
  data = data.frame(JM_dat_wide),
  trans = tmat,
  keep = covs,
  id = "IDAA"
)

JM_msdat_expand <- expand.covs(JM_msdat, covs, append = TRUE, longnames = FALSE)

coxCRfit <- coxph(Surv(Tstart, Tstop, status) ~ SCTyear_2010.3 + hirisk.3 +
                    TCD2.1 + hirisk.1 + TCD2.2 +
                    strata(trans) + cluster(IDAA),
                  data = JM_msdat_expand, x = TRUE, model = TRUE)

## fit mixed-effects sub-model
mixfit <- mvglmer(list(CD8_abs_log ~ ns(intSCT2_2, 3) * TCD2 + (ns(intSCT2_2, 3) | IDAA)),
                  data = dat,
                  families = list(gaussian))

## Specify interactions in order to allow different effect of the outcome for each risk
interacts <- list("CD8_abs_log" = ~ strata(trans) - 1)

## fit multistate model
JMfit <- mvJointModelBayes(mvglmerObject = mixfit,
                           survObject = coxCRfit,
                           timeVar = "intSCT2_2",
                           Interactions = interacts,
                           multiState = TRUE,
                           data_MultiState = JM_msdat_expand,idVar_MultiState =  "IDAA",
                           control = list(equal.strata.knots = TRUE,
                                          equal.strata.bound.knots = TRUE))

# Cannot fit..
summary(JMfit)


# JMBayes2 ----------------------------------------------------------------





