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
    ifelse(TCDmethod %in% c("ALT + 1mg ATG", "ALT + 2mg ATG"), "ALT", "ALT+ATG"),
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


tar_load(dat_merged)

get_postDLI_datasets <- function(dat_merged,
                                 admin_cens_dli, # admin cens post DLI
                                 preDLI_model) {

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
  dat <- dat[uDLI_s == "uDLI"] # Only those with DLI

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
    ifelse(TCDmethod %in% c("ALT + 1mg ATG", "ALT + 2mg ATG"), "ALT", "ALT+ATG"),
    levels = c("ALT", "ALT+ATG")
  )]

  # First predict true current vals from prev model at time of DLI
  dat_wide_temp <- data.table::dcast(
    data = dat,
    formula = IDAA + SCT_May2010 + hirisk + ATG + CMV_PD + uDLI ~ .,
    fun = length
  )

  # Try predli mer


  dat_wide_temp[, intSCT2_5 := uDLI]
  #dat_wide_temp[, CD3_abs_log := 1] # filler

  reffs_lp <- coef(preDLI_model)[, 1:4]
  mmat_newdat <- model.matrix(preDLI_model$termsYx, data = dat_wide_temp)[, 1:4]



  # Levels not yet dropped so can match with as.numeric
  id_new <- as.numeric(dat_wide_temp$IDAA)
  reffs_lp[id_new, ]
  preDLI_model$x$idT
  #preDLI_model$EB$Zb

  as.numeric(NMA_preDLI_datasets$long$IDAA)

  tar_load(NMA_preDLI_datasets)
  id_df <- unique(cbind.data.frame("IDAA" = NMA_preDLI_datasets$long$IDAA, "id" = preDLI_model$id))
  preDLI_model$EB$Zb[match(dat_wide_temp$IDAA, id_df$IDAA)] |>  unique() |>  length()
  preDLI_model
  preDLI_model$data.id

  # Just zb as offset? Try it out

  predict(
    preDLI_model,
    newdata = dat_wide_temp[, intSCT2_5 := uDLI][],
    idVar = "IDAA",
    type = "Marginal" # should be individual!!
  ) +
  cofs <- coe
  mmat <- model.matrix(preDLI_model$termsYx,
               data = dat_wide_temp[, ':=' (
                 intSCT2_5 = uDLI,
                 CD3_abs_log = 1
               )][])

  # Only measurements after DLI

  #... below is older

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
    ifelse(TCDmethod %in% c("ALT + 1mg ATG", "ALT + 2mg ATG"), "ALT", "ALT+ATG"),
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
