library(mstate)
library(survival)
library(JM)


# Data preparation --------------------------------------------------------


tar_load(dat_merged)

dat <- dat_merged[!is.na(CD8_abs_log), c(
  "IDAA",
  "CD8_abs_log",
  "intSCT2_5",
  "endpoint5",
  "endpoint5_s",
  "endpoint_specify5",
  "TCDmethod",
  "hirisk",
  "SCTyear_2010",
  "VCMVPAT_pre"
)]

# Limit to first 12 months after SCT, keep only cell measures pre-endpoint
dat[endpoint5 >= 12, ':=' (endpoint5 = 12, endpoint5_s = "censored")]
dat <- dat[intSCT2_5 < endpoint5] # measures pre-endpoint
dat <- dat[, .SD[.N > 1], by = "IDAA"] # At least 2 measurements per subject

# Combined ATG
dat[, ATG := factor(
  ifelse(TCDmethod == "ALT", "noATG", "yesATG"),
  levels = c("noATG", "yesATG")
)]

dat[, endpoint5_s := factor(
  endpoint5_s,
  levels = c(
    "censored",
    "7 days after cellular intervention",
    "relapse",
    "non-relapse failure: other",
    "non-relapse failure: GvHD"
  ),
  labels = c("cens", "cell_interv", "REL", "NRF_other", "NRF_gvhd")
)]

# Get wide data
JM_dat_wide <- data.table::dcast(
  dat,
  formula = IDAA + SCTyear_2010 + hirisk + ATG + VCMVPAT_pre +
    endpoint5_s + endpoint5 ~ .,
  fun = length
)


# Full attempt with all events + mstate -----------------------------------

event_names <- levels(JM_dat_wide$endpoint5_s)
tmat <- trans.comprisk(
  K = 4,
  names = c("event_free", event_names[-1])
)

tmat

# Make numeric indicators
ind_cols <- paste0("ind_", event_names[-1])
JM_dat_wide[, (ind_cols) := lapply(
  event_names[-1], function(col) as.numeric(endpoint5_s == col)
)]

covs <- c("SCTyear_2010", "hirisk", "ATG", "VCMVPAT_pre")

JM_msdat <- msprep(
  time = c(NA, rep("endpoint5", times = length(ind_cols))),
  status = c(NA, ind_cols),
  data = data.frame(JM_dat_wide),
  trans = tmat,
  keep = covs,
  id = "IDAA"
)

JM_msdat_expand <- expand.covs(
  JM_msdat,
  covs,
  append = TRUE,
  longnames = FALSE
)

# Cannot do left truncation
coxCRfit <- coxph(
  Surv(Tstop, status) ~
    SCTyear_2010.1 + hirisk.1 + # cellular intervention
    ATG.2 + hirisk.2 + # relapse
    ATG.3 + # NRF other; then no covars for gvhd (only current val in JM)
    strata(trans) + cluster(IDAA),
  data = JM_msdat_expand,
  x = TRUE,
  model = TRUE
)

broom::tidy(coxCRfit)

model.matrix(~ strata(trans) - 1, data = coxCRfit$model)
# vs:
model.matrix(~ strata(trans), data = coxCRfit$model)

# See: https://github.com/drizopoulos/JMbayes2/blob/b6088dad58f6e4ad0fbd2b7353d4d1dc3fbfc07a/Development/Dev_Local_GP/MS_CR/Competing_Risks_Reproduce.R

# Fit mixed model (fit simple linear one for now)
lmeFit <- nlme::lme(
  CD8_abs_log ~ ns(intSCT2_5, 3) * ATG + VCMVPAT_pre,
  #random = list(IDAA = pdDiag(form = ~ ns(intSCT2_2, 3))),
  random = ~ ns(intSCT2_5, 3) | IDAA,
  control = lmeControl(opt = 'optim'),
  data = dat
)


# Try to fit joint model now
JMfit <- JM::jointModel(
  lmeObject = lmeFit,
  survObject = coxCRfit,
  timeVar = "intSCT2_5",
  method = "spline-PH-aGH",# or piecewise
  CompRisk = TRUE,
  interFact = list("value" = ~ strata(trans) - 1),
  iter.EM = 200
)

coxCRfit$model$trans2 <- as.numeric(coxCRfit$model$`strata(trans)` == "trans=2")
coxCRfit$model$trans3 <- as.numeric(coxCRfit$model$`strata(trans)` == "trans=3")
coxCRfit$model$trans4 <- as.numeric(coxCRfit$model$`strata(trans)` == "trans=4")


JMfit <- JM::jointModel(
  lmeObject = lmeFit,
  survObject = coxCRfit,
  timeVar = "intSCT2_5",
  method = "spline-PH-aGH",# or piecewise
  CompRisk = TRUE,
  interFact = list("value" = ~ trans2 + trans3 + trans4 - 1, data = coxCRfit$model),
  iter.EM = 200
)

# since strata(trans) is in the model frame of the surv object, not trans!

round(summary(JMfit)$`CoefTable-Long`, 3)
round(summary(lmeFit)$tTable, 3)
round(summary(JMfit)$`CoefTable-Event`, 3)[seq_len(9),]
round(summary(coxCRfit)$coefficients, 3)

library(kableExtra)
dt <- data.frame(round(summary(JMfit)$`CoefTable-Event`, 3)[seq_len(9),]) %>%
  rownames_to_column(var = "Coefficient") %>%
  mutate(
    event_num = str_extract(Coefficient, "[+1-9]$"),
    event = factor(event_num, levels = seq_len(4),
                   labels = c("Cell. intervention", "Relapse", "NRF: Other", "NRF: GVHD"))
  ) %>%
  mutate(
    Coefficient = gsub(x = Coefficient, pattern = "\\:.*", replacement = ""),
    Coefficient = gsub(x = Coefficient, pattern = "\\.[+1-9]$", replacement = ""),
    Coefficient = factor(
      Coefficient,
      levels = c("ATG", "hirisk", "SCTyear_2010", "Assoct"),
      labels = c("ATG", "Disease risk: high", "SCT pre-2010", "Assoct. (CD8)")
    )
  ) %>%
  arrange(event, Coefficient) %>%
  select(
    Event = event,
    Coefficient,
    `log(HR)` = Value,
    "SE" = `Std.Err`,
    "p-value" = `p.value`
  ) %>%
  mutate(
    "HR [95% CI]" = paste0(
      as.character(round(exp(`log(HR)`), 3)), " [",
      as.character(round(exp(`log(HR)` - pnorm(0.975) * SE), 3)), ";",
      as.character(round(exp(`log(HR)` + pnorm(0.975) * SE), 3)), "]"
    )
  ) %>%
  select(-Event)

dt %>%
  kableExtra::kbl(format = "html") %>%
  kable_paper("striped", full_width = F) %>%
  pack_rows("Cell. intervention", 1, 3) %>%
  pack_rows("Relapse", 4, 6) %>%
  pack_rows("NRF: Other", 7, 8) %>%
  pack_rows("NRF: GVHD", 9, 9)

# Plot
newdat <- expand.grid(
  "ATG" = levels(dat$ATG),
  "VCMVPAT_pre" = levels(dat$VCMVPAT_pre),
  "intSCT2_5" = seq(0.1, 12, by = 0.1)
)

preds <- JM::predict.jointModel(
  JMfit,
  newdata = newdat,
  type = "Marginal",
  idVar = "IDAA",
  returnData = TRUE,
  interval = "confidence"
)

dat %>%
  ggplot(
    aes(
      x = intSCT2_5,
      y = CD8_abs_log,
      col = ATG
    )
  ) +
  geom_point(size = 1.5, alpha = 0.75) +
  geom_line(aes(group = IDAA), alpha = 0.5) +
  facet_wrap(~ VCMVPAT_pre) +
  labs(x = "Time since alloHCT (months)", y = "CD8 cell counts") +
  scale_y_continuous(
    breaks = log(c(5, 25, 100, 500, 1500)),
    labels = c(5, 25, 100, 500, 1500)
  ) +
  theme_bw()

library(gghighlight)
dat %>%
  ggplot(
    aes(
      x = intSCT2_5,
      y = CD8_abs_log,
      col = SCTyear_2010
    )
  ) +
  geom_point(size = 1.5, alpha = 0.75) +
  geom_line(aes(group = IDAA), alpha = 0.5) +
  facet_wrap(
    ATG ~ VCMVPAT_pre,
    labeller = labeller(
      ATG = as_labeller(c("noATG" = "no ATG", "yesATG" = "ATG")),
      VCMVPAT_pre = as_labeller(
        x = c("Negative" = "Patient CMV: negative", "Positive" = "Patient CMV: positive")
      )
    )
  ) +
  gghighlight(
    SCTyear_2010 == "pre_2010",
    calculate_per_facet = TRUE,
    use_direct_label = FALSE,
    use_group_by = FALSE
  ) +
  labs(x = "Time since alloHCT (months)", y = "CD8 cell counts") +
  scale_y_continuous(
    breaks = log(c(5, 25, 100, 500, 1500)),
    labels = c(5, 25, 100, 500, 1500)
  ) +
  theme_bw() +
  scale_color_manual(values = "blue")
  #theme(legend.position = "none")


dat %>%
  ggplot(
    aes(
      x = intSCT2_5,
      y = CD8_abs_log,
      col = IDAA
    )
  ) +
  geom_point(size = 1.5, alpha = 0.75) +
  geom_line(aes(group = IDAA), alpha = 0.5) +
  facet_wrap(
    ATG ~ VCMVPAT_pre,
    labeller = labeller(
      ATG = as_labeller(c("noATG" = "no ATG", "yesATG" = "ATG")),
      VCMVPAT_pre = as_labeller(
        x = c("Negative" = "Patient CMV: negative", "Positive" = "Patient CMV: positive")
      )
    )
  ) +
  labs(x = "Time since alloHCT (months)", y = "CD8 cell counts") +
  scale_y_continuous(
    breaks = log(c(5, 25, 100, 500, 1500)),
    labels = c(5, 25, 100, 500, 1500)
  ) +
  theme_bw() +
  theme(legend.position = "none")

preds %>%
  ggplot(
    aes(
      x = intSCT2_5,
      y = pred,
      ymin = low,
      ymax = upp,
      group = ATG,
      col = ATG
    )
  ) +
  geom_ribbon(fill = "gray", alpha = 0.5, col = NA) +
  geom_line(size = 1.5, alpha = 0.75) +
  facet_wrap(~ VCMVPAT_pre, labeller = labeller(
    VCMVPAT_pre = as_labeller(
      x = c("Negative" = "Patient CMV: negative", "Positive" = "Patient CMV: positive")
    )
  )) +
  labs(x = "Time since alloHCT (months)", y = "CD8 cell counts") +
  scale_y_continuous(
    breaks = log(c(5, 25, 100, 500, 1500)),
    labels = c(5, 25, 100, 500, 1500)
  ) +
  theme_bw() +
  coord_cartesian(xlim = c(0, 6))


# Maybe try setting one association to zero?


# JMBayes -----------------------------------------------------------------


## fit mixed-effects sub-model
set.seed(1984)
mixfit <- JMbayes::mvglmer(
  list(CD8_abs_log ~ ns(intSCT2_5, 3) * ATG + VCMVPAT_pre +
         (ns(intSCT2_5, 3) | IDAA)),
  data = as.data.frame(dat), # Cannot be data.table!!
  families = list(gaussian)
)

# Now we need the truncation
coxCRfit_tstart <- coxph(
  Surv(Tstart, Tstop, status) ~
    SCTyear_2010.1 + hirisk.1 + # cellular intervention
    ATG.2 + hirisk.2 + # relapse
    ATG.3 + # NRF other; then no covars for gvhd (only current val in JM)
    strata(trans) + cluster(IDAA),
  data = JM_msdat_expand,
  x = TRUE,
  model = TRUE
)

## fit multistate model - will take ages
JMfit_jmbayes <- JMbayes::mvJointModelBayes(
  mixfit,
  coxCRfit_tstart,
  timeVar = "intSCT2_5",
  Interactions = list("CD8_abs_log" = ~ strata(trans)),
  multiState = TRUE,
  data_MultiState = JM_msdat_expand,
  idVar_MultiState = "IDAA",
  control = list(
    equal.strata.knots = TRUE,
    equal.strata.bound.knots = TRUE
  )
)

summary(JMfit_jmbayes)


# JMBayes2 ----------------------------------------------------------------



#fforms <- list("CD8_abs_log" = ~ value(CD8_abs_log) + value(CD8_abs_log):(strata(trans) - 1))
#fforms <- list("CD8_abs_log" = ~ value(CD8_abs_log):strata(trans)) # not implausible either
#fforms <- list("CD8_abs_log" = ~ value(CD8_abs_log):trans) #nope
#fforms <- list("CD8_abs_log" = ~ value(CD8_abs_log):(strata(trans) - 1)) # most similar output

fitjm2 <- JMbayes2::jm(
  Surv_object = coxCRfit,
  Mixed_objects = list(lmeFit),
  time_var = 'intSCT2_5',
  data_Surv = JM_msdat_expand,
  id_var = 'IDAA',
  functional_forms = fforms
)

summary(fitjm2)$Survival
round(summary(JMfit)$`CoefTable-Event`, 3)[seq_len(9),]










