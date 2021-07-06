#' Returns long and wide dataset (not yet for DLI)
get_datasets <- function(dat_merged) {

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

  # Return them in list
  res <- list("long" = dat, "wide" = dat_wide)

  return(res)
}


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

  # Edit tmat now
  tmat_new <- trans.comprisk(3, names = c("event-free", "REL", "GVHD", "NRF_other"))
  attr(JM_msdat, "trans") <- tmat_new
  JM_msdat_expand <- mstate::expand.covs(
    JM_msdat,
    c(covs, "DLI"),
    append = TRUE,
    longnames = FALSE
  )

  return(JM_msdat_expand)
}


run_longitudinal_submodels <- function(dat) {

  cell_vars <- grep(x = colnames(dat), pattern = "_log$", value = TRUE)

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

  return(long_submodels)
}

run_cox_submodel <- function(form, dli_msdata) {

  mod_comp <- do.call(
    survival::coxph,
    list(
      "formula" = form,
      "cluster" = dli_msdata$IDAA,
      "model" = TRUE,
      "x" = TRUE,
      "data" = dli_msdata
    )
  )

  return(mod_comp)
}


run_JM <- function(long_obj,
                   surv_obj,
                   fform) {

  # Set up model matrices of assoc. parameters
  surv_obj$model$trans1 <- as.numeric(surv_obj$model$`strata(trans)` == "trans=1")
  surv_obj$model$trans2 <- as.numeric(surv_obj$model$`strata(trans)` == "trans=2")
  surv_obj$model$trans3 <- as.numeric(surv_obj$model$`strata(trans)` == "trans=3")

  # Fit jointmodel
  jm_mod <- jointModel(
    lmeObject = long_obj,
    survObject = surv_obj,
    CompRisk = TRUE,
    timeVar = "intSCT2_5",
    method = "spline-PH-aGH",
    interFact = list(
      "value" = fform
    ),
    iter.EM = 200 # see p.68 rizop book
  )

  return(jm_mod)
}
