# Data preparation for submodels ------------------------------------------


#' Prepare wide and long format datasets
#'
#' The long format dataset is used to run the mixed models, while the wide
#' data is a) used mainly for descriptives, and b) used to later prepare
#' data in format to analyse intermediate DLI.
#'
#' @param dat_merged Dataset that merges raw lympocyte and variable data
#' using `prepare_raw_data()`.
#'
#' @return A list containing long and wide datasets
get_datasets <- function(dat_merged, admin_cens = 24) {

  vars <- c(
    "IDAA",
    "CD4_abs_log",
    "CD8_abs_log",
    "CD3_abs_log",
    "CD19_abs_log",
    "NK_abs_log",
    "intSCT2_5",
    "endpoint5",
    "endpoint5_s",
    "endpoint_specify5",
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
    "CMV_PD"
  )

  dat <- dat_merged[, ..vars]

  # Drop unused levels in subsetted cohort
  factors <- which(vapply(dat, FUN = is.factor, FUN.VALUE = logical(1)))
  dat[, (factors) := lapply(.SD, droplevels), .SDcols = factors]

  # Admin censoring at 2y post-HSCT
  dat[endpoint6 >= admin_cens, ':=' (
    endpoint6 = admin_cens,
    endpoint6_s = "censored"
  )]

  # Keep measurements prior to endpoint
  dat <- dat[intSCT2_5 < endpoint6]

  # Check this later: (more than one measurement at single time point)
  #dat <- dat[, .SD[!duplicated(intSCT2_5)], by = IDAA]

  # Keep patients with at least 2 measurements, and no missing in cell counts
  cell_vars <- grep(x = colnames(dat), pattern = "_log$", value = TRUE)
  dat <- na.omit(dat[, .SD[.N >= 2], by = "IDAA"], cols = cell_vars)

  # Prepare a few variables
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
    formula = IDAA + SCT_May2010 + hirisk + ATG + CMV_PD +
      endpoint6_s + endpoint6 + uDLI + uDLI_s + endpoint_specify6 +
      TCD + TCD2 + TCDmethod ~ .,
    fun = length
  )
  data.table::setnames(dat_wide, old = ".", new = "n_measurements")

  # Return them in list
  res <- list("long" = dat, "wide" = dat_wide)

  return(res)
}


#' Prepare wide data to accomodate intermediate DLI and competing risks
#'
#' Uses `mstate::msprep()` to prepare dataset for analysis of intermediate
#' DLI. This is done by preparing for a competing risks multi-state analysis
#' with DLI as intermediate state, and thereafter grouping transitions (e.g
#' event-free -> REL and DLI -> REL are grouped as a single transition).
#'
#' See the tutorial by Putter et al. for further explanations.
#'
#' @param dat_wide Wide (single row per patient) data as returned by
#' `get_datasets()`.
#'
#' @return A data.frame of class "msdata"
prepare_dli_msdata <- function(dat_wide) {

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
    x = list(c(2, 3, 4, 5), c(3, 4, 5), c(), c(), c()),
    names = c("event-free", "DLI", "REL", "GVHD", "NRF_other")
  )

  # Msprep data
  covs <- c("CMV_PD", "ATG", "SCT_May2010", "hirisk")
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
  JM_msdat <- JM_msdat[order(JM_msdat$IDAA, JM_msdat$trans), ]

  # Edit tmat now
  tmat_new <- trans.comprisk(K = 3, names = c("event-free", "REL", "GVHD", "NRF_other"))
  attr(JM_msdat, "trans") <- tmat_new
  JM_msdat_expand <- mstate::expand.covs(
    data = JM_msdat,
    covs = c(covs, "DLI"),
    append = TRUE,
    longnames = FALSE
  )

  # Add transition indicators (for functional forms)
  JM_msdat_expand$trans1 <- as.numeric(JM_msdat_expand$trans == 1)
  JM_msdat_expand$trans2 <- as.numeric(JM_msdat_expand$trans == 2)
  JM_msdat_expand$trans3 <- as.numeric(JM_msdat_expand$trans == 3)

  return(JM_msdat_expand)
}


# Wrappers for running submodels ------------------------------------------


#' Run `lme()`-based longitudinal submodels
#'
#' These are the models to be fed into `jointModel()` and `jm()`.
#'
#' @param dat Long-format dataset as returned by `get_datasets()`
#' @param which_cells Lymphocytes names for which to run the models. E.g.
#' c("CD8_abs_log", "CD4_abs_log").
#' @param df_splines Degrees of freedom to use for the splines of time. Default
#' is 4, which corresponds to 3 internal knots.
#' @param ranef_structure Structure of random effect vcov matrix, either
#' "diagonal" or "unstructured"
#'
#' @return A mixed model class `lme()`
run_longitudinal_submodels <- function(dat,
                                       which_cells,
                                       df_splines = 4,
                                       ranef_structure = "diagonal") {

  # Retrieve all (log) cell-line names
  cell_vars <- grep(x = colnames(dat), pattern = "_log$", value = TRUE)
  names(cell_vars) <- cell_vars

  if (!missing(which_cells)) {
    if (any(!(which_cells %in% cell_vars))) stop("One or more of cell names invalid.")
    cell_vars <- cell_vars[cell_vars %in% which_cells]
  }

  # Set up mixed model structure
  spline_string <- paste0("ns(intSCT2_5, df = ", df_splines, ")")
  mod_rhs <- paste0(spline_string, " * ATG + CMV_PD")
  ranef_expr <- if (ranef_structure == "diagonal") {
    list(IDAA = pdDiag(as.formula(paste0("~ ", spline_string))))
  } else as.formula(paste0("~ ", spline_string, "| IDAA"))

  long_submodels <- lapply(cell_vars, function(cell) {

    form <- stats::reformulate(
      response = cell,
      termlabels = mod_rhs
    )

    lmeFit <- do.call(
      what = nlme::lme,
      args = list(
        "fixed" = form,
        "random" = ranef_expr,
        "control" = nlme::lmeControl(opt = "optim"),
        "data" = dat
      )
    )

    return(lmeFit)
  })

  return(long_submodels)
}


#' Run survival submodel
#'
#' @param form Formula to use for the Cox model. Transition-specific variables
#' are suffixed by the competing event number, e.g. DLI.1 for DLI predictor for
#' Relapse.
#' @param dli_msdata Dataset returned by `prepare_dli_msdata()`.
#' @param ... A list of additional arguments to be passed to `survival::coxph()`
#'
#' @return A `coxph()` model
run_cox_submodel <- function(form, dli_msdata, ...) {

  extra_args <- list(...)
  main_args <- list(
    "formula" = form,
    "cluster" = dli_msdata$IDAA,
    "model" = TRUE,
    "x" = TRUE,
    "data" = dli_msdata
  )

  mod_comp <- do.call(
    what = survival::coxph,
    args = c(main_args, extra_args)
  )

  return(mod_comp)
}
