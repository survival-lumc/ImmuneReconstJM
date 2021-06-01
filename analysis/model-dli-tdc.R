##**********************************##
## Initial model(s) with DLI as TDC ##
##**********************************##


library(data.table)
library(JMbayes)
library(JMbayes2)
library(ggalluvial)
library(ggsankey)
library(tidyverse)
library(ggforce)
library(mstate)

source("data-raw/prepare-raw-data.R")
source("R/individual-cell-models.R")

variables_raw <- data.table::data.table(readRDS("data-raw/2021-05-31_v6/variables.rds"))
lymphocytes_raw <- data.table::data.table(readRDS("data-raw/2021-05-31_v6/lymphocytes.rds"))
dat_merged <- prepare_raw_data(lymphocytes_raw, variables_raw)


# Prepare data ------------------------------------------------------------

vars <- c(
  "IDAA",
  "CD4_abs_log",
  "CD8_abs_log",
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
dat <- dat[!(is.na(CD4_abs_log) | is.na(CD8_abs_log))]
dat[, ATG := factor(
  ifelse(TCDmethod == "ALT", "noATG", "yesATG"),
  levels = c("noATG", "yesATG")
)]

dat[, uDLI_s := factor(uDLI_s, c(0, 1), c("none", "uDLI"))]
dat[, post_DLI := factor(as.numeric(uDLI <= intSCT2_5), c(0, 1), c("pre_DLI", "post_DLI"))]
dat[, endpoint6_s := factor(
  endpoint6_s,
  levels = c(
    "censored",
    "7 days after cellular intervention",
    "relapse",
    "non-relapse failure: other",
    "non-relapse failure: GvHD"
  ),
  labels = c("cens", "cell_interv", "REL", "NRF_other", "NRF_gvhd")
)]


# Make the wide dataset
dat_wide <- data.table::dcast(
  data = dat,
  formula = IDAA + SCT_May2010 + hirisk + ATG + VCMVPAT_pre +
    endpoint6_s + endpoint6 + uDLI + uDLI_s + endpoint_specify6 ~ .,
  fun = length
)
data.table::setnames(dat_wide, old = ".", new = "n_measurements")


# Sankey plot (detailed) --------------------------------------------------


dat_wide_trans <- data.table::dcast(
  data = dat_merged,
  formula = IDAA + endpoint5_s + uDLI_s + sec_endpoint_s +
    sec_endpoint_specify + endpoint_specify5 + endpoint6_s + endpoint_specify6 ~ .,
  fun = length
)
# trans_dat <- dat_wide_trans[, .(n = .N), by = c(
#   "endpoint5_s", "endpoint_specify5", "sec_endpoint_specify"
# )]

dat_wide_trans %>%
  make_long(endpoint5_s, endpoint_specify5, sec_endpoint_specify) %>%
  filter(!(is.na(node) & is.na(next_x) & is.na(next_node))) %>%
  ggplot(
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
  #scale_x_discrete(limits = c("1st endpoint", "1st endpoint (detailed)", "Post-DLI")) +
  theme_sankey(base_size = 18) +
  labs(x = NULL) +
  theme(legend.position = "none")


# Sankey plot (mstate-like) -----------------------------------------------


dat_flow <- data.table::copy(dat_wide_trans)
dat_flow[, ':=' (
  HSCT = "HSCT",
  uDLI_s = factor(uDLI_s, c(0, 1), c("none", "uDLI")),
  endp_comp = endpoint6_s,
  endp_specify = ifelse(
    endpoint6_s == as.character(endpoint_specify6), NA, as.character(endpoint_specify6)
  )
)]

dat_flow %>%
  make_long(HSCT, uDLI_s, endp_comp, endp_specify) %>%
  filter(!(is.na(node) & is.na(next_x) & is.na(next_node))) %>%
  ggplot(
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
  scale_fill_viridis_d(direction = -1) +
  theme_sankey(base_size = 18) +
  labs(x = NULL) +
  theme(legend.position = "none")



# Longitudinal plots ------------------------------------------------------


dat[uDLI_s == "uDLI"] %>%
  ggplot(aes(intSCT2_5, CD4_abs_log, group = IDAA)) +
  geom_line(alpha = 0.5) +
  geom_point(aes(col = post_DLI)) +
  facet_grid(~ endpoint6_s, scales = "free_x") +
  theme(legend.position = "none")



# Event timings post-DLI --------------------------------------------------


# No admin censoring
table(dat_wide$endpoint6_s)
table(dat_wide$endpoint6_s, dat_wide$uDLI_s)

# Maybe at 2 years?
dat_wide %>%
  ggplot(aes(x = endpoint6)) +
  geom_histogram(bins = 30, fill = "lightblue", col = "black") +
  facet_grid(uDLI_s ~ endpoint6_s) +
  theme_bw()


# Number measurements post-DLI
dat[, .(n_measurements = .N), by = c("post_DLI", "IDAA", "endpoint6_s")] %>%
  .[, .(mu_postDLI = mean(n_measurements)), by = c("post_DLI", "endpoint6_s")]


# Cause-specific of GVHD --------------------------------------------------


# Cause-specific model of gvhd
dat_wide[, gvhd_s := ifelse(endpoint6_s == "NRF_gvhd", 1, 0)]
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

# Try relative time
dat[, t_from_dli := round((uDLI <= intSCT2_5) * (intSCT2_5 - uDLI), 4)]

# Try same with DAT (see supplement)
lmeFit <- lme(
  fixed = CD4_abs_log ~ ns(intSCT2_5, 3) * ATG + VCMVPAT_pre,
  random = ~ ns(intSCT2_5, 3) | IDAA,
  control = lmeControl(opt = "optim"),
  data = dat
)

# Does not converge with splines
lmeFit_tdep <- lme(
  fixed = CD4_abs_log ~ intSCT2_5 * ATG + VCMVPAT_pre + t_from_dli,
  random = ~ intSCT2_5 + t_from_dli | IDAA,
  control = lmeControl(opt = "optim", msMaxIter = 100),
  data = data.frame(dat)
)

jointFit <- jm(mod_cox, lmeFit_tdep, time_var = "intSCT2_5")
summary(jointFit)

newdat <- data.frame(dat[IDAA == "221"])

preds <- predict(
  jointFit, process = "longitudinal",
  type = "mean_subject",
  return_newdata = TRUE,
  newdata = newdat,
  times = seq(0.1, 6, by = 0.1)
)

ggtraceplot(jointFit, "betas") #long
ggtraceplot(jointFit, "alphas") #association
ggtraceplot(jointFit, "gammas") #survival


# JMbayes -----------------------------------------------------------------

mixed_fits <- mvglmer(
  list(
    CD4_abs_log ~ ns(intSCT2_5, 3) * ATG + VCMVPAT_pre + (ns(intSCT2_5, 3) | IDAA),
    CD8_abs_log ~ ns(intSCT2_5, 3) * ATG + VCMVPAT_pre + (ns(intSCT2_5, 3) | IDAA)
  ),
  data = dat,
  families = list(gaussian, gaussian)
)

jointFit_gvhd <- mvJointModelBayes(
  mixed_fits,
  mod_cox,
  timeVar = "intSCT2_5"#,
  #Interactions = list("CD4_abs_log" = ~ DLI, "CD8_abs_log" = ~ DLI)
)
summary(jointFit_gvhd)

newdat <- expand.grid(
  "ATG" = levels(dat$ATG),
  "VCMVPAT_pre" = levels(dat$VCMVPAT_pre),
  "intSCT2_5" = seq(0.1, 6, by = 0.1)
)

preds <- predict(
  object = jointFit_gvhd,
  newdata = newdat,
  type = "Marginal",
  returnData = TRUE,
  interval = "confidence"
)

preds %>%
  ggplot(aes(intSCT2_5, pred, col = ATG)) +
  geom_line() +
  facet_grid(. ~ VCMVPAT_pre)

# JMbayes - tdep ----------------------------------------------------------

jointFit_gvhd <- jointModelBayes(lmeFit_tdep, mod_cox, timeVar = "intSCT2_5")
summary(jointFit_gvhd)

newdat <- expand.grid(
  "ATG" = levels(dat$ATG),
  "VCMVPAT_pre" = levels(dat$VCMVPAT_pre),
  "intSCT2_5" = seq(0.1, 12, by = 0.1),
  "t_from_dli" = 0
)

t_dli <- 6
newdat$t_from_dli <- (newdat$intSCT2_5 >= t_dli) * (newdat$intSCT2_5 - t_dli)

preds <- predict(
  object = jointFit_gvhd,
  newdata = newdat,
  type = "Marginal",
  returnData = TRUE,
  interval = "confidence"
)

preds %>%
  ggplot(aes(intSCT2_5, pred, col = ATG)) +
  geom_line() +
  facet_grid(. ~ VCMVPAT_pre)


# Try with comp risks -----------------------------------------------------


dat_wide[, cr_ind := ifelse(
  endpoint6_s %in% c("REL", "NRF_gvhd"),
  as.character(endpoint6_s), "cens"
)]

dat_wide[, ':=' (
  cr_ind = factor(
    x = cr_ind,
    levels = c("cens", "REL", "NRF_gvhd"),
  ),
  REL_ind = as.numeric(cr_ind == "REL"),
  GVHD_ind = as.numeric(cr_ind == "NRF_gvhd"),
  DLI_ind = as.numeric(uDLI_s) - 1
)]


# Now use mstate
tmat <- transMat(
  list(c(2, 3, 4), c(3, 4), c(), c()),
  names = c("event-free", "DLI", "REL", "GVHD")
)

covs <- c("VCMVPAT_pre", "ATG", "SCT_May2010", "hirisk")
JM_msdat <- mstate::msprep(
  time = c(NA, "uDLI", "endpoint6", "endpoint6"),
  status = c(NA, "DLI_ind", "REL_ind", "GVHD_ind"),
  data = data.frame(dat_wide),
  trans = tmat,
  keep = covs,
  id = "IDAA"
)

# Group transitions and set TDC
JM_msdat <- JM_msdat[JM_msdat$trans != 1, ]
JM_msdat$trans <- ifelse(JM_msdat$trans %in% c(2, 4), 1, JM_msdat$trans)
JM_msdat$trans <- ifelse(JM_msdat$trans %in% c(3, 5), 2, JM_msdat$trans)
JM_msdat$DLI <- as.numeric(JM_msdat$Tstart > 0)
JM_msdat$from[JM_msdat$from == 2] <- 1
JM_msdat$to[JM_msdat$to == 4] <- 2
JM_msdat <- JM_msdat[order(JM_msdat$IDAA, JM_msdat$trans),]
View(JM_msdat)

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

# JMbayes cannot handle exogenous + competing risks (splitting error due to Ids)
# You would have to explictly model intermediate state


# With JMbayes2
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
