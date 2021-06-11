##**********************************##
## Initial model(s) with DLI as TDC ##
##**********************************##


library(data.table)
library(JM)
library(JMbayes)
library(JMbayes2)
library(ggsankey)
library(tidyverse)
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
  "CD3_abs_log",
  "CD19_abs_log",
  "NK_abs_log",
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
dat[endpoint6 >= 24, ':=' (
  endpoint6 = 24,
  endpoint6_s = "censored"
)]

dat <- dat[intSCT2_5 < endpoint6]

# Check this later: (more than one measurment at single time point)
#dat <- dat[, .SD[!duplicated(intSCT2_5)], by = IDAA]

dat <- dat[, .SD[.N >= 2], by = "IDAA"]
dat <- dat[!(is.na(CD4_abs_log) | is.na(CD8_abs_log))] # change this
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

# Drop IDAA levels
dat[, IDAA := droplevels(IDAA)]


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

ggalluvial::geom
ggforce::paral

# Longitudinal plots ------------------------------------------------------


set.seed(1083) # Heterog
set.seed(1084) # More homog
dat[uDLI_s == "uDLI"][IDAA %in% sample(IDAA, size = 6, replace = FALSE)] %>%
  melt.data.table(
    id.vars = c("post_DLI", "intSCT2_5", "IDAA"),
    measure.vars = c("CD4_abs_log", "CD8_abs_log", "NK_abs_log"),
    variable.name = "lymphocyte",
    value.name = "count"
  ) %>%
  ggplot(aes(intSCT2_5, count)) +
  geom_point(
    aes(fill = post_DLI),
    size = 3.5, pch = 21,
    alpha = 0.8,
    col = "#359fda",
    #fill = colorspace::lighten("#359fda", amount = 0.3)
  ) +
  geom_path(alpha = 0.25, linetype = "dotted") +
  facet_grid(IDAA ~ lymphocyte, scales = "free_x") +
  labs(x = "Time since alloHCT (months)", y = "CD8 cell counts") +
  scale_y_continuous(
    breaks = log(c(5, 25, 100, 500, 1500)),
    labels = c(5, 25, 100, 500, 1500)
  ) +
  theme_minimal() +
  guides(fill = guide_legend(title = "DLI"))



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
mod_cox

# Try relative time
dat[, t_from_dli := round((uDLI <= intSCT2_5) * (intSCT2_5 - uDLI), 4)]

# Try same with DAT (see supplement)
#dat <- dat[intSCT2_5 < 0]
lmeFit <- lme(
  fixed = CD4_abs_log ~ ns(intSCT2_5, 3) * ATG + VCMVPAT_pre,
  random = ~ ns(intSCT2_5, 3) | IDAA,
  control = lmeControl(opt = "optim"),#, msMaxIter = 100),
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
    CD8_abs_log ~ ns(intSCT2_5, 3) * ATG + VCMVPAT_pre + (ns(intSCT2_5, 3) | IDAA)
  ),
  data = dat,
  families = list(gaussian)
)

jointFit_gvhd <- jointModelBayes(
  lmeFit,
  mod_cox,
  timeVar = "intSCT2_5",
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
    DLI.2 + strata(trans), # add sct and hirisk.2
  cluster = IDAA,
  model = TRUE,
  x = TRUE,
  data = JM_msdat_expand
)

mod_comp

# JMbayes cannot handle exogenous + competing risks (splitting error due to Ids)
# You would have to explictly model intermediate state


# With JMbayes2
ff <- list(
  "CD4_abs_log" = ~ value(CD4_abs_log) + value(CD4_abs_log):(strata(trans) - 1)
) # Add interaction

ff2 <- list(
  "CD4_abs_log" = ~ DLI * (value(CD4_abs_log) + value(CD4_abs_log):(strata(trans) - 1))
)

jointFit_jmb2 <- jm(
  Surv_object = mod_comp,
  Mixed_objects = lmeFit,
  time_var = "intSCT2_5",
  functional_forms = ff,
  id_var = "IDAA"
)

summary(jointFit_jmb2)



# Try again with JM -------------------------------------------------------


lmeFit
mod_cox

mod_cox <- coxph(
  Surv(tstart, tstop, status) ~ DLI + cluster(IDAA),
  data = df_long,
  x = TRUE,
  model = TRUE
)

# Need to arrange data otherwise
dat[, "gvhd_s" := as.numeric(endpoint6_s == "NRF_gvhd")]
dat2 <- data.frame(dat[, c("IDAA", "CD4_abs_log", "uDLI_s", "uDLI", "intSCT2_5",
              "gvhd_s", "endpoint6")])

dat2$start <- dat2$intSCT2_5

splitID <- split(dat2[c("start", "endpoint6")], dat2$IDAA)
dat2$stop <- unlist(lapply(splitID,
                          function(d) c(d$start[-1], d$endpoint6[1]) ))

dat2$event <- with(dat2, ave(gvhd_s, IDAA,
                           FUN = function (x) c(rep(0, length(x)-1), x[1]) ))
dat2$dli <- as.numeric(dat2$stop > dat2$uDLI)
mod_cox2 <- coxph(
  Surv(start, stop, event) ~ dli + cluster(IDAA),
  data = dat2,
  x = TRUE,
  model = TRUE
)

mod_cox2

jm_mod <- jointModel(
  lmeFit,
  mod_cox,
  timeVar = "intSCT2_5",
  method = "spline-PH-aGH", # Christ almighty is this important
  interFact = list("value" = ~ DLI),
  iter.EM = 200
)
summary(jm_mod)

# With comprisk - wont work because of counting notation?
jm_mod2 <- jointModel(
  lmeFit,
  mod_comp,
  CompRisk = TRUE,
  timeVar = "intSCT2_5",
  method = "spline-PH-aGH", # Christ almighty is this important
  interFact = list("value" = ~ strata(trans) - 1),
  #interFact = list("value" = ~ `DLI.1` * (strata(trans) - 1)),
  iter.EM = 200 # see p.68 rizop book
)
summary(jm_mod2)

newdat <- expand.grid(
  "ATG" = levels(dat_wide$ATG),
  "VCMVPAT_pre" = levels(dat_wide$VCMVPAT_pre),
  "intSCT2_5" = seq(0.1, 12, by = 0.1)
)

JM::predict.jointModel(
  jm_mod2,
  newdata = newdat,
  type = "Marginal",
  idVar = "IDAA",
  returnData = TRUE,
  interval = "confidence"
) %>%
  ggplot(aes(x = intSCT2_5, y = pred, group = ATG, col = ATG)) +
  geom_ribbon(
    aes(ymin = low, ymax = upp),
    fill = "gray",
    alpha = 0.5,
    col = NA
  ) +
  geom_line(size = 1.5, alpha = 0.75) +
  #geom_text(aes(label = label), hjust = 0, na.rm = TRUE, fontface = "bold") +
  facet_wrap(facets = ~ VCMVPAT_pre) +
  labs(x = "Time since alloHCT (months)", y = "CD8 cell counts") +
  scale_y_continuous(
    breaks = log(c(5, 25, 100, 500, 1500)),
    labels = c(5, 25, 100, 500, 1500)
  ) +
  theme_minimal() +
  theme(legend.position = "none") +
  coord_cartesian(xlim = c(0, 14)) +
  scale_color_brewer(palette = "Paired")



# Test with inter ---------------------------------------------------------


mod_comp$model$trans1 <- as.numeric(mod_comp$model$`strata(trans)` == "trans=1") #rel
mod_comp$model$trans2 <- as.numeric(mod_comp$model$`strata(trans)` == "trans=2")
model.matrix(~ trans1 + trans1:DLI.1 + trans2 + trans2:DLI.2 - 1,
             data = mod_comp$model) %>%  View()

jm_mod3 <- jointModel(
  lmeFit,
  mod_comp,
  CompRisk = TRUE,
  timeVar = "intSCT2_5",
  method = "spline-PH-aGH", # Christ almighty is this important
  interFact = list("value" = ~ trans1 + trans1:DLI.1 + trans2 + trans2:DLI.2 - 1),
  iter.EM = 200 # see p.68 rizop book
)

summary(jm_mod3)
#summary(jm_mod2)

JM::predict.jointModel(
  jm_mod3,
  newdata = newdat,
  type = "Marginal",
  idVar = "IDAA",
  returnData = TRUE,
  interval = "confidence"
) %>%
  ggplot(aes(x = intSCT2_5, y = pred, group = ATG, col = ATG)) +
  geom_ribbon(
    aes(ymin = low, ymax = upp),
    fill = "gray",
    alpha = 0.5,
    col = NA
  ) +
  geom_line(size = 1.5, alpha = 0.75) +
  #geom_text(aes(label = label), hjust = 0, na.rm = TRUE, fontface = "bold") +
  facet_wrap(facets = ~ VCMVPAT_pre) +
  labs(x = "Time since alloHCT (months)", y = "CD8 cell counts") +
  scale_y_continuous(
    breaks = log(c(5, 25, 100, 500, 1500)),
    labels = c(5, 25, 100, 500, 1500)
  ) +
  theme_minimal() +
  #theme(legend.position = "none") +
  coord_cartesian(xlim = c(0, 14)) +
  scale_color_brewer(palette = "Paired")
