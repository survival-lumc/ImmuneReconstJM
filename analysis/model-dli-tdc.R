library(data.table)
library(JMbayes2)
library(ggalluvial)
library(ggsankey)
library(tidyverse)
library(ggforce)

source("data-raw/prepare-raw-data.R")
source("R/individual-cell-models.R")

variables_raw <- data.table::data.table(readRDS("data-raw/2021-05-31_v6/variables.rds"))
lymphocytes_raw <- data.table::data.table(readRDS("data-raw/2021-05-31_v6/lymphocytes.rds"))
dat_merged <- prepare_raw_data(lymphocytes_raw, variables_raw)


# Prepare data ------------------------------------------------------------

vars <- c(
  "IDAA",
  "CD4_abs_log",
  "intSCT2_5",
  "endpoint6",
  "endpoint6_s",
  "endpoint_specify6",
  "uDLI",
  "uDLI_s",
  "sec_endpoint",
  "sec_endpoint_s",
  "sec_endpoint_specify",
  "TCDmethod",
  "hirisk",
  "SCT_May2010",
  "VCMVPAT_pre"
)

dat <- dat_merged[, ..vars]
dat <- dat[intSCT2_5 < endpoint6]
dat <- dat[, .SD[.N >= 2], by = "IDAA"]
dat <- dat[!is.na(CD4_abs_log)]
dat[, ATG := factor(
  ifelse(TCDmethod == "ALT", "noATG", "yesATG"),
  levels = c("noATG", "yesATG")
)]


# dat[, endpoint5_s := factor(
#   endpoint5_s,
#   levels = c(
#     "censored",
#     "7 days after cellular intervention",
#     "relapse",
#     "non-relapse failure: other",
#     "non-relapse failure: GvHD"
#   ),
#   labels = c("cens", "cell_interv", "REL", "NRF_other", "NRF_gvhd")
# )]

# Make the wide dataset
dat_wide <- data.table::dcast(
  data = dat,
  formula = IDAA + SCT_May2010 + hirisk + ATG + VCMVPAT_pre +
    endpoint6_s + endpoint6 + uDLI + uDLI_s + endpoint_specify6 ~ .,
  fun = length
)
data.table::setnames(dat_wide, old = ".", new = "n_measurements")



# Descriptives ------------------------------------------------------------

dat_wide_trans <- data.table::dcast(
  data = dat_merged,
  formula = IDAA + endpoint5_s + uDLI_s + sec_endpoint_s +
    sec_endpoint_specify + endpoint_specify5 + endpoint6_s ~ .,
  fun = length
)
data.table::setnames(dat_wide_trans, old = ".", new = "n_measurements")
trans_dat <- dat_wide_trans[, .(n = .N), by = c(
  "endpoint5_s", "endpoint_specify5", "sec_endpoint_specify"
)]
trans_dat_long <- gather_set_data(data = trans_dat, x = 1:3)

ggplot(trans_dat_long, aes(x, id = id, split = y, value = n)) +
  geom_parallel_sets(alpha = 0.3, axis.width = 0.1) +
  geom_parallel_sets_axes(axis.width = 0.1) +
  geom_parallel_sets_labels(colour = 'white')

trans_dat %>%
  ggplot(aes(y = n, axis1 = uDLI_s, axis2 = endpoint6_s, axis3 = endpoint_specify6)) +
  geom_alluvium(aes(fill = endpoint_specify6)) +
  guides(fill = FALSE) +
  geom_stratum() +
  geom_text(stat = "stratum", aes(label = after_stat(stratum)))

df <- mtcars %>%
  make_long(cyl, vs, am, gear, carb) %>%
  arrange(x)

df_test <- trans_dat %>%
  make_long(endpoint5_s, endpoint_specify5, sec_endpoint_specify, value = n)

df_test <- dat_wide_trans %>%
  make_long(endpoint5_s, endpoint_specify5, sec_endpoint_specify) %>%
  filter(!(is.na(node) & is.na(next_x) & is.na(next_node)))

ggplot(
  df_test,
  aes(
    x = x,
    next_x = next_x,
    node = node,
    next_node = next_node,
    fill = factor(node),
    label = node
  )
) +
  geom_sankey(flow.alpha = .6, node.color = "gray30") +
  geom_sankey_label(size = 3, color = "white", fill = "gray40") +
  scale_fill_viridis_d() +
  theme_sankey(base_size = 18) +
  labs(x = NULL) +
  theme(legend.position = "none")



# Alluvial version --------------------------------------------------------

trans_dat <- dat_wide_trans[, .(n = .N), by = c(
  "endpoint6_s", #"endpoint_specify5",
  "endpoint5_s"
)]

ggplot(trans_dat,
       aes(y = n, axis1 = endpoint5_s, axis2 = endpoint6_s)) +
  geom_alluvium(aes(fill = endpoint6_s), width = 1/12) +
  geom_stratum(width = 1/12, fill = "gray", color = "black") +
  geom_label(stat = "stratum", aes(label = after_stat(stratum))) +
  #geom_text(aes(label = n), stat = "flow", size = 3) +
  scale_x_discrete(limits = c("Pre-DLI", "Post-DLI"), expand = c(.05, .05)) +
  scale_fill_brewer(type = "qual", palette = "Set1") +
  theme(legend.position = "none")

# Attempt with JM ---------------------------------------------------------


# Cause-specific model of gvhd
dat_wide[, gvhd_s := ifelse(
  endpoint6_s == "non-relapse failure: GvHD", 1, 0
)]

dat_wide[, id := IDAA]

df_long <- survival::tmerge(
  data1 = data.frame(dat_wide),
  data2 = data.frame(dat_wide),
  id = id,
  status = event(endpoint6, gvhd_s),
  DLI = tdc(uDLI)
)

mod_cox <- coxph(
  Surv(tstart, tstop, status) ~ DLI,
  cluster = id,
  data = df_long,
  x = TRUE,
  model = TRUE
)

lmeFit <- lme(
  fixed = CD4_abs_log ~ ns(intSCT2_5, 3) * ATG + VCMVPAT_pre,
  random = ~ ns(intSCT2_5, 3) | IDAA,
  control = lmeControl(opt = "optim"),
  data = dat
)

jointFit_gvhd <- jointModel(
  lmeObject = lmeFit,
  survObject = mod_cox,
  timeVar = "intSCT2_5",
  method = "spline-PH-aGH",
  iter.EM = 200
)



# JMbayes -----------------------------------------------------------------

jointFit_gvhd <- jointModelBayes(lmeFit, mod_cox, timeVar = "intSCT2_5", )
preds <- predict(jointFit_gvhd, newdata = newdat, type = "Marginal",
                 returnData = TRUE,
                 interval = "confidence")

preds %>%
  ggplot(aes(intSCT2_5, pred, col = ATG)) +
  geom_line() +
  facet_grid(. ~ VCMVPAT_pre)


# JMbayes2 ----------------------------------------------------------------



library(JMbayes2)
jointFit_jmb2 <- jm(
  Surv_object = mod_cox,
  Mixed_objects = lmeFit,
  time_var = "intSCT2_5"
)

ggtraceplot(jointFit_jmb2, "betas") #long
ggtraceplot(jointFit_jmb2, "alphas") #association
ggtraceplot(jointFit_jmb2, "gammas") #survival


# Prediction
newdat <- expand.grid(
  "ATG" = levels(dat_wide$ATG),
  "VCMVPAT_pre" = levels(dat_wide$VCMVPAT_pre),
  "intSCT2_5" = seq(0.1, 12, by = 0.1)
)

preds <- predict(
  jointFit_jmb2,
  #newdata = newdat,
  return_newdata = TRUE
)



# Try with comp risks -----------------------------------------------------


dat_wide[, cr_ind := ifelse(
  endpoint6_s %in% c("non-relapse failure: GvHD", "relapse"),
  endpoint6_s, "censored"
)]
dat_wide[, ':=' (
  cr_ind = factor(
    x = cr_ind,
    levels = c("censored", "relapse", "non-relapse failure: GvHD"),
    labels = c("CENS", "REL", "GVHD")
  ),
  REL_ind = as.numeric(cr_ind == "relapse"),
  GVHD_ind = as.numeric(cr_ind == "non-relapse failure: GvHD"),
  event_ind = as.numeric(endpoint6_s %in% c("non-relapse failure: GvHD", "relapse")) # for TDC prep
)]

dat_wide[, id := IDAA]


# Now use mstate
library(mstate)
tmat <- transMat(
  list(c(2, 3, 4), c(3, 4), c(), c()),
  names = c("event-free", "DLI", "REL", "GVHD")
)

covs <- c("VCMVPAT_pre", "ATG", "SCTyear_2010", "hirisk")
JM_msdat <- mstate::msprep(
  time = c(NA, "uDLI", "endpoint6", "endpoint6"),
  status = c(NA, "uDLI_s", "REL_ind", "GVHD_ind"),
  data = data.frame(dat_wide),
  trans = tmat,
  keep = covs,
  id = "IDAA"
)

# Group transitions and set TDC
JM_msdat <- JM_msdat[JM_msdat$trans != 1, ]
JM_msdat$trans <- ifelse(JM_msdat$trans %in% c(2, 4), 1, JM_msdat$trans)
JM_msdat$trans <- ifelse(JM_msdat$trans %in% c(3, 5), 2, JM_msdat$trans)
JM_msdat$DLI <- ifelse(JM_msdat$Tstart > 0, 1, 0)
JM_msdat$from[JM_msdat$from == 2] <- 1
JM_msdat$to[JM_msdat$to == 4] <- 2
JM_msdat <- JM_msdat[order(JM_msdat$IDAA, JM_msdat$trans),]

tmat_new <- trans.comprisk(2, names = c("event-free", "REL", "GVHD"))
attr(JM_msdat, "trans") <- tmat_new
JM_msdat_expand <- mstate::expand.covs(
  JM_msdat,
  c(covs, "DLI"),
  append = TRUE,
  longnames = FALSE
)

mod_comp <- coxph(
  Surv(Tstart, Tstop, status) ~ hirisk.1 + DLI.1 + ATG.1 +
    DLI.2 + strata(trans),
  cluster = IDAA,
  model = TRUE,
  data = JM_msdat_expand
)

mod_comp

library(JMbayes)
mixfit <- mvglmer(
  list(CD4_abs_log ~ ns(intSCT2_5, 3) * ATG + VCMVPAT_pre + (ns(intSCT2_5, 3) | IDAA)),
  data = dat,
  families = list(gaussian)
)

mod_comp$strata <- mod_comp$model$`strata(trans)`

jointFit_tdc <- mvJointModelBayes(
  mixfit,
  mod_comp,
  timeVar = "intSCT2_5",
  Interactions = list("CD4_abs_log" = ~ strata(trans) - 1),
  multiState = TRUE,
  idVar_MultiState = "IDAA"
)



preds <- predict(jointFit_gvhd, newdata = newdat, type = "Marginal",
                 returnData = TRUE,
                 interval = "confidence")

preds %>%
  ggplot(aes(intSCT2_5, pred, col = ATG)) +
  geom_line() +
  facet_grid(. ~ VCMVPAT_pre)


# JMbayes2 ----------------------------------------------------------------

library(JMbayes2)
ff <- list(
  "CD4_abs_log" = ~ value(CD4_abs_log) + value(CD4_abs_log):(strata(trans) - 1)
)

jointFit_jmb2 <- jm(
  Surv_object = mod_comp,
  Mixed_objects = lmeFit,
  time_var = "intSCT2_5",
  functional_forms = ff,
  id_var = "IDAA"
)

summary(jointFit_jmb2)

newdat <- expand.grid(
  "ATG" = levels(dat_wide$ATG),
  "VCMVPAT_pre" = levels(dat_wide$VCMVPAT_pre)#,
  #"intSCT2_5" = seq(0.1, 12, by = 0.1)
)

predict(
  jointFit_jmb2,
  newdata = newdat,
  times = seq(0.1, 6, by = 0.1),
  return_newdata = TRUE
)


# try with cause-specific GVHD


# Notes:
# - JM simple code for exogenous and simple event no longer works (maybe due to survival)
# - JMbayes with multi-state does not work
# - JMbayes2 works, predictions functions not there yet

# To-do: model with cause-specific gvhd and time-dep DLI, and modulated long. event model
# (See code Greg in paper supplement)
# Some plots for Wed
