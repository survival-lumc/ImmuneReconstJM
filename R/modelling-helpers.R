# dli model needs to be connected after with model..

#' Prepare wide and long format datasets for models before/up to DLI
#'
#' The long format dataset is used to run the mixed models, while the wide
#' data is a) used mainly for descriptives, and b) used for survival submodel
#'
#' @param dat_merged Dataset that merges raw lympocyte and variable data
#' using `prepare_raw_data()`.
#' @param admin_cens Time of administrative censoring (months)
#'
#' @return A list containing long and wide datasets
get_preDLI_datasets <- function(dat_merged,
                                admin_cens) {

  # Variables to keep
  vars <- c(
    "IDAA",
    "CD4_abs_log",
    "CD8_abs_log",
    "CD3_abs_log",
    "intSCT2_5",
    "endpoint7",
    "endpoint7_s",
    "endpoint_specify7",
    "endpoint6",
    "endpoint6_s",
    "endpoint_specify6",
    "uDLI",
    "uDLI_s",
    "TCDmethod",
    "TCD",
    "TCD2",
    "hirisk",
    "SCT_May2010",
    "CMV_PD",
    "earlylow_DLI",
    "HLAmismatch_GvH",
    "relation"
  )

  dat <- dat_merged[, ..vars]

  # Define the endpoint
  dat[, endpoint7_s := factor(
    endpoint7_s,
    levels = c(
      "censored",
      "non-relapse failure: GvHD",
      "relapse",
      "non-relapse failure: other",
      "7 days after cellular intervention"
    ),
    labels = c("cens", "gvhd", "relapse", "other_nrf", "cens")
  )]

  # Apply admin censoring
  dat[endpoint7 >= admin_cens, ':=' (
    endpoint7 = admin_cens,
    endpoint7_s = "cens"
  )]

  # Measurements at endpoint taken as just prior
  dat[intSCT2_5 == endpoint7, intSCT2_5 := intSCT2_5 - 0.01]
  dat <- dat[intSCT2_5 < endpoint7]

  # Remove (couple) of missings from cell variables
  cell_vars <- grep(x = colnames(dat), pattern = "_log$", value = TRUE)
  dat <- na.omit(dat, cols = cell_vars)

  # Drop factor levels
  factors <- which(vapply(dat, FUN = is.factor, FUN.VALUE = logical(1L)))
  dat[, (factors) := lapply(.SD, droplevels), .SDcols = factors]

  # Add ATG variable (only relevant NMA)
  dat[, ATG := factor(
    ifelse(TCDmethod %in% c("ALT + 1mg ATG", "ALT + 2mg ATG"), "ALT+ATG", "ALT"),
    levels = c("ALT", "ALT+ATG")
  )]

  # Make the wide dataset
  dat_wide <- data.table::dcast(
    data = dat,
    formula = IDAA + SCT_May2010 + hirisk + ATG + CMV_PD +
      endpoint7_s + endpoint7 + endpoint_specify7 +
      HLAmismatch_GvH + relation +
      uDLI + uDLI_s + earlylow_DLI + TCD + TCD2 + TCDmethod ~ .,
    fun = length
  )
  data.table::setnames(dat_wide, old = ".", new = "n_measurements")

  # Return them in list
  res <- list("long" = dat, "wide" = dat_wide)

  return(res)
}


#' Run survival submodel
#'
#' @param form Formula to use for the Cox model. Transition-specific variables
#' are suffixed by the competing event number, e.g. DLI.1 for DLI predictor for
#' GVHD
#' @param dat_wide Dataset returned by `get_preDLI_datasets()`.
#' @param ... A list of additional arguments to be passed to `survival::coxph()`
#'
#' @return A `coxph()` model
run_preDLI_cox <- function(form, dat_wide, ...) {

  # Prepare data
  tmat <- trans.comprisk(K = 3, names = c("gvhd", "relapse", "other_nrf"))
  covs <- c("CMV_PD", "hirisk", "ATG", "earlylow_DLI",
            "HLAmismatch_GvH", "relation", "SCT_May2010")

  msdat <- msprep(
    time = c(NA, rep("endpoint7", 3)),
    status = with(
      dat_wide, cbind(
        NA,
        1 * (endpoint7_s == "gvhd"),
        1 * (endpoint7_s == "relapse"),
        1 * (endpoint7_s == "other_nrf")
      )
    ),
    data = data.frame(dat_wide),
    trans = tmat,
    keep = covs,
    id = "IDAA"
  )

  msdat_expand <- expand.covs(msdat, covs, append = TRUE, longnames = FALSE)

  extra_args <- list(...)
  main_args <- list(
    "formula" = form,
    "cluster" = msdat_expand$IDAA,
    "model" = TRUE,
    "x" = TRUE,
    "data" = msdat_expand
  )

  mod_comp <- do.call(what = survival::coxph, args = c(main_args, extra_args))

  return(mod_comp)
}


#' Run `lme()`-based longitudinal submodels
#'
#' These are the models to be fed into `jointModel()` and `jm()`.
#'
#' @param dat Long-format dataset as returned by `get_datasets()`
#' @param cell_line Lymphocytes name for which to run the models. E.g.
#' "CD4_abs_log"
#' @param form_fixed Character, rhs of fixed effects formula. Will be
#' appended to cell_line
#' @param form_random Specification of random effects (formula object)
#' @param ... A list of additional arguments to be passed to `nlme::lme()`
#'
#' @return A mixed model class `lme()`
run_preDLI_longitudinal <- function(cell_line,
                                    form_fixed,
                                    form_random,
                                    dat,
                                    ...) {

  # Combine formulas
  form <- stats::reformulate(response = cell_line, termlabels = form_fixed)

  # Run fit
  extra_args <- list(...)
  main_args <- list(
    "fixed" = form,
    "random" = form_random,
    "control" = nlme::lmeControl(opt = "optim", msMaxIter = 200),
    "data" = dat
  )

  lmeFit <- do.call(what = nlme::lme, args = c(main_args, extra_args))

  return(lmeFit)
}


# Continue here with long, JM, jmbayes2 helpers..
# + code to send to shark
# + clean up ALL other code



# Joint modeling helpers --------------------------------------------------



run_jointModel <- function(long_obj,
                           surv_obj,
                           ...) {

  # Set up model matrices of assoc. parameters
  surv_obj$model$trans1 <- as.numeric(surv_obj$model$`strata(trans)` == "trans=1")
  surv_obj$model$trans2 <- as.numeric(surv_obj$model$`strata(trans)` == "trans=2")
  surv_obj$model$trans3 <- as.numeric(surv_obj$model$`strata(trans)` == "trans=3")

  # Fit jointmodel
  mod <- jointModel(
    lmeObject = long_obj,
    survObject = surv_obj,
    CompRisk = TRUE,
    method = "spline-PH-aGH",
    ...
  )

  return(mod)
}


# jmbayes 2 here..


# Post-DLI helpers --------------------------------------------------------


get_postDLI_datasets <- function(dat_merged,
                                 admin_cens_dli = 12, # admin cens months post DLI
                                 preDLI_model,
                                 preDLI_datasets) {

  # Variables to keep
  vars <- c(
    "IDAA",
    "CD4_abs_log",
    "CD8_abs_log",
    "CD3_abs_log",
    "intSCT2_5",
    "endpoint7",
    "endpoint7_s",
    "endpoint_specify7",
    "sec_endpoint",
    "sec_endpoint_s",
    "sec_endpoint_specify",
    "uDLI",
    "uDLI_s",
    "TCDmethod",
    "TCD",
    "TCD2",
    "hirisk",
    "SCT_May2010",
    "CMV_PD",
    "earlylow_DLI",
    "HLAmismatch_GvH",
    "relation"
  )

  dat <- dat_merged[, ..vars]

  # Define the endpoint - rel and nrf are merged
  dat[, sec_endpoint_s := factor(
    sec_endpoint_s,
    levels = c(
      "censored",
      "non-relapse failure: GvHD",
      "relapse",
      "non-relapse failure: other"
    ),
    labels = c("cens", "gvhd", "rel_nrf", "rel_nrf")
  )]

  # Measurements at endpoint taken as just prior
  dat[intSCT2_5 == sec_endpoint, intSCT2_5 := intSCT2_5 - 0.01]

  # Add ATG variable (only relevant NMA)
  dat[, ATG := factor(
    ifelse(TCDmethod %in% c("ALT + 1mg ATG", "ALT + 2mg ATG"), "ALT+ATG", "ALT"),
    levels = c("ALT", "ALT+ATG")
  )]

  # First predict true current vals from prev model at time of DLI
  # (for now used fitted values at timepoint just prior)
  preDLI_dat_long <- preDLI_datasets$long
  preDLI_dat_long[, preds_subj := fitted(
    preDLI_model,
    process = "Longitudinal",
    type = "Subject"
  )]
  setorder(preDLI_dat_long, "IDAA", "intSCT2_5")
  df_at_dli <- preDLI_dat_long[
    uDLI_s == "uDLI" & intSCT2_5 <= uDLI, .SD[.N], by = "IDAA"
    ][, c("IDAA", "preds_subj")]

  # Selection happens
  dat <- dat[uDLI_s == "uDLI"] # Only those with DLI

  # Admin censoring
  dat[sec_endpoint >= uDLI + admin_cens_dli, ':=' (
    sec_endpoint = uDLI + admin_cens_dli,
    sec_endpoint_s = "cens"
  )]
  dat <- dat[intSCT2_5 < sec_endpoint]

  # Remember to clock reset measurements!!
  dat <- dat[uDLI < intSCT2_5]
  dat[, intSCT2_5_reset := intSCT2_5 - uDLI]

  # Remove (couple) of missings from cell variables
  cell_vars <- grep(x = colnames(dat), pattern = "_log$", value = TRUE)
  dat <- na.omit(dat, cols = cell_vars)

  # Drop factor levels
  factors <- which(vapply(dat, FUN = is.factor, FUN.VALUE = logical(1L)))
  dat[, (factors) := lapply(.SD, droplevels), .SDcols = factors]

  # Make the wide dataset
  dat_wide <- data.table::dcast(
    data = dat,
    formula = IDAA + SCT_May2010 + hirisk + ATG + CMV_PD +
      sec_endpoint_s + sec_endpoint + sec_endpoint_specify +
      HLAmismatch_GvH + relation +
      uDLI + uDLI_s + earlylow_DLI + TCD + TCD2 + TCDmethod ~ .,
    fun = length
  )
  data.table::setnames(dat_wide, old = ".", new = "n_measurements")

  # Add dli col for measure
  dat_wide <- merge(dat_wide, df_at_dli, by = "IDAA")

  # Return them in list
  res <- list("long" = dat, "wide" = dat_wide)

  return(res)
}


# Testing here.. remove after
tar_load(
  c(
    dat_merged,
    NMA_preDLI_datasets,
    preDLI_CD3__jointModel_both
  )
)

NMA_dats_postDLI <- get_postDLI_datasets(
  dat_merged = dat_merged[TCD2 %in% c("NMA RD: ALT", "UD: ALT + ATG")],
  admin_cens_dli = 12,
  preDLI_CD3__jointModel_both,
  NMA_preDLI_datasets
)

dat_wide <- NMA_dats_postDLI$wide

run_postDLI_cox <- function(form, dat_wide, ...) {

  # Prepare data
  tmat <- trans.comprisk(K = 2, names = c("DLI", "GVHD", "REL_NRF"))
  covs <- c("CMV_PD", "hirisk", "ATG", "earlylow_DLI",
            "HLAmismatch_GvH", "relation", "SCT_May2010", "preds_subj")

  dat_wide_prepped <- copy(dat_wide)
  msdat <- msprep(
    time = c(NA, "sec_endpoint", "sec_endpoint"),
    status = with(
      dat_wide, cbind(
        NA,
        1 * (sec_endpoint_s == "gvhd"),
        1 * (sec_endpoint_s == "rel_nrf")
      )
    ),
    start = list(state = rep(1, nrow(dat_wide_prepped)), time = dat_wide_prepped$uDLI),
    data = data.frame(dat_wide_prepped),
    trans = tmat,
    keep = covs,
    id = "IDAA"
  )

  msdat_expand <- expand.covs(msdat, covs, append = TRUE, longnames = FALSE)

  extra_args <- list(...)
  main_args <- list(
    "formula" = form,
    "cluster" = msdat_expand$IDAA,
    "model" = TRUE,
    "x" = TRUE,
    "data" = msdat_expand
  )

  mod_comp <- do.call(what = survival::coxph, args = c(main_args, extra_args))

  return(mod_comp)
}

mod_comp <- run_postDLI_cox(
  form = Surv(time, status) ~
    earlylow_DLI.1 + preds_subj.1  + # GVHD
    hirisk.2 + preds_subj.2 + # REL and NRF
    strata(trans),
  dat_wide = dat_wide
)

lmeFit <- run_preDLI_longitudinal(
  cell_line = "CD3_abs_log",
  form_fixed = "ns(intSCT2_5_reset, 2) * earlylow_DLI + ATG", # add three-way interaction later
  form_random = list(IDAA = pdDiag(~ ns(intSCT2_5_reset, 2))),
  dat = NMA_dats_postDLI$long
)


# Make raw plots
NMA_dats_postDLI$long |>
  ggplot(aes(intSCT2_5_reset, CD3_abs_log, col = IDAA, group = IDAA)) +
  geom_line() +
  facet_grid(earlylow_DLI ~ sec_endpoint_s) +
  theme(legend.position = "none")

jm_fit_both <- jointModel(
  lmeObject = lmeFit,
  survObject = mod_comp,
  CompRisk = TRUE,
  timeVar = "intSCT2_5_reset",
  method = "spline-PH-aGH",
  derivForm = list(
    fixed = ~ 0 + dns(intSCT2_5_reset, 2) + dns(intSCT2_5_reset, 2):as.numeric(earlylow_DLI == "yes"),
    random = ~ 0 + dns(intSCT2_5_reset, 2),
    indFixed = c(2:3, 6:7),
    indRandom = c(2:3)
  ),
  interFact = list(
    "value" = ~ strata(trans) - 1,
    "slope" = ~ strata(trans) - 1
  ),
  parameterization = "both",
  control = list("iter.EM" = 1000) #increasuuhhh
)

jm_fit_both |>  summary()
# Plot fitteeedddd
# histogram of slopes?
# just curr value

dat_long <- NMA_dats_postDLI$long
dat_long[, preds_subj := fitted(
  jm_fit_both,
  process = "Longitudinal",
  type = "Subject"
)]


set.seed(20220119)
IDAA_subs <- sample(levels(dat_long$IDAA), replace = FALSE, size = 32)#size = 56)

ggplot(
  dat_long, aes(intSCT2_5_reset, CD3_abs_log)
) +
  geom_point() +
  geom_line(aes(y = preds_subj, group = IDAA)) +
  facet_wrap(~ IDAA) +
  theme(legend.position = "none")
