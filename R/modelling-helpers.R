#' Prepare wide and long format datasets for models before/up to DLI
#'
#' The long format dataset is used to run the mixed models, while the wide
#' data is a) used mainly for descriptives, and b) used for survival submodel
#'
#' @param dat_merged Dataset that merges raw lymphocyte and variable data
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
    "intSCT2_7",
    "endpoint7",
    "endpoint7_s",
    "endpoint_specify7",
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
    x = endpoint7_s,
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

  # Measurements are terminated by endpoint
  dat <- dat[intSCT2_7 < endpoint7]

  # Remove missings from cell variables early on in follow-up
  cell_vars <- grep(x = colnames(dat), pattern = "_log$", value = TRUE)
  dat <- na.omit(dat, cols = cell_vars)

  # Drop factor levels
  factors <- which(vapply(dat, FUN = is.factor, FUN.VALUE = logical(1L)))
  dat[, (factors) := lapply(.SD, droplevels), .SDcols = factors]

  # Add ATG/UD variable (only relevant for NMA cohort)
  dat[, ATG := factor(
    ifelse(TCDmethod %in% c("ALT + 1mg ATG", "ALT + 2mg ATG"), "UD(+ATG)", "RD"),
    levels = c("RD", "UD(+ATG)")
  )]

  # Make the wide dataset
  dat_wide <- data.table::dcast(
    data = dat,
    formula = IDAA + SCT_May2010 + hirisk + ATG + CMV_PD + endpoint7_s + endpoint7 +
      endpoint_specify7 + HLAmismatch_GvH + relation + uDLI + uDLI_s + earlylow_DLI +
      TCD + TCD2 + TCDmethod ~ .,
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
  covs <- c(
    "CMV_PD",
    "hirisk",
    "ATG",
    "earlylow_DLI",
    "HLAmismatch_GvH",
    "relation",
    "SCT_May2010"
  )

  # Prepare "long" data to run cause-specific hazard models
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
    "control" = nlme::lmeControl(opt = "optim", msMaxIter = 200L),
    "data" = dat
  )

  lmeFit <- do.call(what = nlme::lme, args = c(main_args, extra_args))

  return(lmeFit)
}


# Continue here with long, JM, jmbayes2 helpers..
# + code to send to shark
# + clean up ALL other code



# Joint modeling helpers --------------------------------------------------


# Need this to specify different interactions/slopes in different transitions
tweak_preDLI_modmat <- function(cox_model) {
  cox_model$model$trans1 <- as.numeric(cox_model$model$`strata(trans)` == "trans=1")
  cox_model$model$trans2 <- as.numeric(cox_model$model$`strata(trans)` == "trans=2")
  cox_model$model$trans3 <- as.numeric(cox_model$model$`strata(trans)` == "trans=3")
  return(cox_model)
}


# Post-DLI helpers --------------------------------------------------------


get_postDLI_datasets <- function(dat_merged,
                                 admin_cens_dli = 3) { # Three months after early low-dose

  # Variables to keep
  vars <- c(
    "IDAA",
    "CD4_abs_log",
    "CD8_abs_log",
    "CD3_abs_log",
    "intDLI1",
    "intSCT2_7",
    "endpoint7",
    "endpoint7_s",
    "endpoint_specify7",
    "sec_endpoint2",
    "sec_endpoint2_s",
    "sec_endpoint2_specify",
    "uDLI_actual",
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
  dat[, sec_endpoint2_s := factor(
    sec_endpoint2_s,
    levels = c(
      "censored",
      "7 days after cellular intervention",
      "non-relapse failure: GvHD",
      "relapse",
      "non-relapse failure: other"
    ),
    labels = c("cens", "cens", "gvhd", "rel_nrf", "rel_nrf")
  )]

  # Add ATG variable (only relevant NMA)
  dat[, ATG := factor(
    ifelse(TCDmethod %in% c("ALT + 1mg ATG", "ALT + 2mg ATG"), "UD(+ATG)", "RD"),
    levels = c("RD", "UD(+ATG)")
  )]

  # Selection happens
  dat <- dat[earlylow_DLI == "yes"] # Only those with early-low DLI

  # Endpoint needs to be in t since DLI
  dat[, sec_endpoint2 := sec_endpoint2 - uDLI_actual]

  # Admin censoring
  dat[sec_endpoint2 >= admin_cens_dli, ':=' (
    sec_endpoint2 = admin_cens_dli,
    sec_endpoint2_s = "cens"
  )]

  # Keep measurements between DLI and endpoint
  dat <- dat[intDLI1 < sec_endpoint2 & intDLI1 >= 0]

  # Remove missings from cell variables
  cell_vars <- grep(x = colnames(dat), pattern = "_log$", value = TRUE)
  dat <- na.omit(dat, cols = cell_vars)

  # Drop factor levels
  factors <- which(vapply(dat, FUN = is.factor, FUN.VALUE = logical(1L)))
  dat[, (factors) := lapply(.SD, droplevels), .SDcols = factors]

  # Make the wide dataset
  dat_wide <- data.table::dcast(
    data = dat,
    formula = IDAA + SCT_May2010 + hirisk + ATG + CMV_PD + sec_endpoint2_s + sec_endpoint2 +
      sec_endpoint2_specify + HLAmismatch_GvH + relation + uDLI + uDLI_s + earlylow_DLI +
      TCD + TCD2 + TCDmethod ~ .,
    fun = length
  )
  data.table::setnames(dat_wide, old = ".", new = "n_measurements")

  # Return them in list
  res <- list("long" = dat, "wide" = dat_wide)

  return(res)
}


# Same principle as preDLI function
run_postDLI_cox <- function(form, dat_wide, ...) {

  # Prepare data
  tmat <- trans.comprisk(K = 2, names = c("DLI", "GVHD", "REL_NRF"))
  covs <- c("ATG")

  dat_wide_prepped <- copy(dat_wide)
  msdat <- msprep(
    time = c(NA, "sec_endpoint2", "sec_endpoint2"),
    status = with(
      dat_wide, cbind(
        NA,
        1 * (sec_endpoint2_s == "gvhd"),
        1 * (sec_endpoint2_s == "rel_nrf")
      )
    ),
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
