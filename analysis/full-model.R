##***************************************##
## (To be made into an rmd with targets) ##
##***************************************##


# Load packages
library(data.table)
library(JM)
library(JMbayes2)
#library(ggsankey)
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
  "TCDmethod",
  "hirisk",
  "SCT_May2010",
  "VCMVPAT_pre"
)

dat <- dat_merged[, ..vars]

# Admin censoring at 2y post-HSCT
dat[endpoint6 >= 24, ':=' (
  endpoint6 = 24,
  endpoint6_s = "censored"
)]

# Keep measurements prior to endpoint
dat <- dat[intSCT2_5 < endpoint6]

# Check this later: (more than one measurment at single time point)
#dat <- dat[, .SD[!duplicated(intSCT2_5)], by = IDAA]

# Keep patients with at least 2 measurements, and no missing in cell counts
cell_vars <- grep(x = colnames(dat), pattern = "_log$", value = TRUE)
dat <- na.omit(dat[, .SD[.N >= 2], by = "IDAA"], cols = cell_vars)

# Prepare a few variables
dat[, uDLI_s := factor(uDLI_s, c(0, 1), c("none", "uDLI"))]

dat[, ATG := factor(
  ifelse(TCDmethod == "ALT", "noATG", "yesATG"),
  levels = c("noATG", "yesATG")
)]

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


# Fit separate long submodels ---------------------------------------------


names(cell_vars) <- cell_vars
long_submodels <- lapply(cell_vars, function(cell) {

  form <- stats::reformulate(
    response = cell,
    termlabels = "splines::ns(intSCT2_5, 3) * ATG + VCMVPAT_pre"
  )

  lmeFit <- do.call(
    nlme::lme,
    list(
      "fixed" = form,
      "random" = ~ splines::ns(intSCT2_5, 3) | IDAA,
      "control" = nlme::lmeControl(opt = "optim"),
      "data" = dat
    )
  )

  return(lmeFit)
})


# Competing risks submodel ------------------------------------------------


# Make cell interv a nuisance
dat_wide[, cr_ind := ifelse(
  endpoint6_s != "cell_interv",
  as.character(endpoint6_s), "cens"
)]

# Prep factors for msprep
dat_wide[, ':=' (
  cr_ind = factor(
    x = cr_ind,
    levels = c("cens", "REL", "NRF_gvhd", "NRF_other"),
  ),
  REL_ind = as.numeric(cr_ind == "REL"),
  GVHD_ind = as.numeric(cr_ind == "NRF_gvhd"),
  NRF_ind = as.numeric(cr_ind == "NRF_other"),
  DLI_ind = as.numeric(uDLI_s) - 1
)]

# Prepare transition mat (DLI as intermediate state)
tmat <- transMat(
  list(c(2, 3, 4, 5), c(3, 4, 5), c(), c(), c()),
  names = c("event-free", "DLI", "REL", "GVHD", "NRF_other")
)
tmat

# Msprep data
covs <- c("VCMVPAT_pre", "ATG", "SCT_May2010", "hirisk")
JM_msdat <- mstate::msprep(
  time = c(NA, "uDLI", "endpoint6", "endpoint6", "endpoint6"),
  status = c(NA, "DLI_ind", "REL_ind", "GVHD_ind", "NRF_ind"),
  data = data.frame(dat_wide),
  trans = tmat,
  keep = covs,
  id = "IDAA"
)

# Group transitions and set TDC
JM_msdat <- JM_msdat[JM_msdat$trans != 1, ]
JM_msdat$trans <- ifelse(JM_msdat$trans %in% c(2, 5), 1, JM_msdat$trans)
JM_msdat$trans <- ifelse(JM_msdat$trans %in% c(3, 6), 2, JM_msdat$trans)
JM_msdat$trans <- ifelse(JM_msdat$trans %in% c(4, 7), 3, JM_msdat$trans)
JM_msdat$DLI <- as.numeric(JM_msdat$Tstart > 0)

# Edit to and from state (without DLI)
JM_msdat$from[JM_msdat$from == 2] <- 1
JM_msdat$to[JM_msdat$to == 3] <- 2
JM_msdat$to[JM_msdat$to == 4] <- 3
JM_msdat$to[JM_msdat$to == 5] <- 4
JM_msdat <- JM_msdat[order(JM_msdat$IDAA, JM_msdat$trans),]
View(JM_msdat)

# Edit tmat now
tmat_new <- trans.comprisk(3, names = c("event-free", "REL", "GVHD", "NRF_other"))
attr(JM_msdat, "trans") <- tmat_new
JM_msdat_expand <- mstate::expand.covs(
  JM_msdat,
  c(covs, "DLI"),
  append = TRUE,
  longnames = FALSE
)

mod_comp <- coxph(
  Surv(Tstart, Tstop, status) ~
    hirisk.1 + DLI.1 + ATG.1 + # Relapse submodel
    DLI.2 + ATG.2 + # GVHD submodel
    DLI.3 + ATG.3 + # NRF_other submodel
    strata(trans),
  cluster = IDAA,
  model = TRUE,
  x = TRUE,
  data = JM_msdat_expand
)

mod_comp


# Separate models with JM -------------------------------------------------


# Set up model matrices of assoc. parameters
mod_comp$model$trans1 <- as.numeric(mod_comp$model$`strata(trans)` == "trans=1")
mod_comp$model$trans2 <- as.numeric(mod_comp$model$`strata(trans)` == "trans=2")
mod_comp$model$trans3 <- as.numeric(mod_comp$model$`strata(trans)` == "trans=3")
model.matrix(
  ~ trans1 + trans1:DLI.1 + trans2 + trans2:DLI.2 + trans3 + trans3:DLI.3 - 1,
  data = mod_comp$model
) %>%  View()

# ---- Run for CD4
jm_mod_CD4 <- jointModel(
  lmeObject = long_submodels$CD4_abs_log,
  survObject = mod_comp,
  CompRisk = TRUE,
  timeVar = "intSCT2_5",
  method = "spline-PH-aGH",
  interFact = list(
    "value" = ~ trans1 + trans1:DLI.1 + trans2 + trans2:DLI.2 + trans3 + trans3:DLI.3 - 1
  ),
  iter.EM = 200 # see p.68 rizop book
)

summary(jm_mod_CD4)

newdat <- expand.grid(
  "ATG" = levels(dat_wide$ATG),
  "VCMVPAT_pre" = levels(dat_wide$VCMVPAT_pre),
  "intSCT2_5" = seq(0.1, 24, by = 0.1)
)


JM::predict.jointModel(
  jm_mod_CD4,
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


# Other JM ----------------------------------------------------------------

table(dat_wide$uDLI_s, dat_wide$endpoint6_s)

mod_comp_bis <- coxph(
  Surv(Tstart, Tstop, status) ~
    hirisk.1 + ATG.1 + # Relapse submodel
    DLI.2 + ATG.2 + # GVHD submodel
    ATG.3 + # NRF_other submodel
    strata(trans),
  cluster = IDAA,
  model = TRUE,
  x = TRUE,
  data = JM_msdat_expand
)

mod_comp_bis

# Set up model matrices of assoc. parameters
mod_comp_bis$model$trans1 <- as.numeric(mod_comp_bis$model$`strata(trans)` == "trans=1")
mod_comp_bis$model$trans2 <- as.numeric(mod_comp_bis$model$`strata(trans)` == "trans=2")
mod_comp_bis$model$trans3 <- as.numeric(mod_comp_bis$model$`strata(trans)` == "trans=3")
model.matrix(
  ~ trans1 + trans2 + trans2:DLI.2 + trans3 - 1,
  data = mod_comp_bis$model
) %>%  View()

# ---- Run for CD4
jm_mod_CD4 <- jointModel(
  lmeObject = long_submodels$CD4_abs_log,
  survObject = mod_comp_bis,
  CompRisk = TRUE,
  timeVar = "intSCT2_5",
  method = "spline-PH-aGH",
  interFact = list(
    "value" = ~  trans1 + trans2 + trans2:DLI.2 + trans3 - 1
  ),
  iter.EM = 200 # see p.68 rizop book
)

summary(jm_mod_CD4) # note zero implies succesful convergence!

jm_mod_CD4_2 <- jointModel(
  lmeObject = long_submodels$CD4_abs_log,
  survObject = mod_comp,
  CompRisk = TRUE,
  timeVar = "intSCT2_5",
  method = "spline-PH-aGH",
  interFact = list(
    "value" = ~ trans1 + trans2 + trans2:DLI.2 + trans3 - 1
  ),
  iter.EM = 200 # see p.68 rizop book
)

summary(jm_mod_CD4_2)
jm_mod_CD4_2$coefficients$gammas
mod_comp$coefficients
long_submodels$CD4_abs_log$coefficients$fixed
jm_mod_CD4_2$coefficients$betas

jm_mod_CD8 <- jointModel(
  lmeObject = long_submodels$CD8_abs_log,
  survObject = mod_comp_bis,
  CompRisk = TRUE,
  timeVar = "intSCT2_5",
  method = "spline-PH-aGH",
  interFact = list(
    "value" = ~  trans1 + trans2 + trans2:DLI.2 + trans3 - 1
  ),
  iter.EM = 200 # see p.68 rizop book
)

summary(jm_mod_CD8)


# Multivariate model with JMbayes2 ----------------------------------------

# With JMbayes2
ff <- list(
  "CD4_abs_log" = ~ value(CD4_abs_log) + value(CD4_abs_log):(strata(trans) - 1)
) # Add interaction


# Separate first
ff <- list(
  "CD4_abs_log" = ~ DLI * (value(CD4_abs_log) + value(CD4_abs_log):(strata(trans) - 1))
)

ff <- list(
  "CD4_abs_log" = ~ DLI * (value(CD4_abs_log):(strata(trans) - 1))
) #?? no clue which one is correct, not the jm one did not converge


jointFit_jmb2 <- jm(
  Surv_object = mod_comp,
  Mixed_objects = long_submodels$CD4_abs_log,
  time_var = "intSCT2_5",
  functional_forms = ff,
  id_var = "IDAA",
  control = list(
    n_burnin = 2000,
    n_iter = 10000
  )
)

summary(jointFit_jmb2)
jointFit_jmb2$statistics$Median$alphas
jm_mod_CD4$coefficients$alpha

ggdensityplot(jointFit_jmb2, "alphas")
ggtraceplot(jointFit_jmb2, "betas")
ggtraceplot(jointFit_jmb2, "gammas")
ggtraceplot(jointFit_jmb2, "alphas")

table(dat_wide$DLI_ind, dat_wide$endpoint6_s)

# Maybe later model just GVHD and Rel; or omit DLI in one?

# + marginal prediction functions..
