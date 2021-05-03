#' Prepare Leiden data for joint modeling
#'
#' @param merged_data (Raw) dataset containing all variables
#' @param cells Character vector of cells to use - if more than two used,
#' it is in preparation for multivariate longitudinal model
#' @param admin_cens_time Time at which to apply administrative censoring (months)
#' @param control_n_measurements Keep only patients with at least this number
#' of cell count measurments (recommended at least 2 to avoid software errors)
prepare_JM_data <- function(merged_data,
                            cells,
                            admin_cens_time = 12,
                            control_n_measurements = 2) {

  # Subset relevant variables
  vars <- c(
    "IDAA",
    cells,
    "intSCT2_5",
    "endpoint5",
    "endpoint5_s",
    "endpoint_specify5",
    "TCDmethod",
    "hirisk",
    "SCTyear_2010",
    "VCMVPAT_pre"
  )

  dat <- merged_data[, .SD, .SDcols = vars]

  # Keep complete values
  expr_keep_complete <- paste0("!is.na(", cells, ")")
  if (length(cells) > 1) expr_keep_complete <- paste(expr_keep_complete, collapse = "&")
  dat <- dat[eval(parse(text = expr_keep_complete))]

  # Admin censoring
  dat[endpoint5 >= admin_cens_time, ':=' (
    endpoint5 = admin_cens_time,
    endpoint5_s = "censored"
  )]

  # What to do with patients where intSCT2_5 == endpoint5 (all with relapse)?
  # ... Probably take measurement prior
  dat <- dat[intSCT2_5 < endpoint5]
  dat <- dat[, .SD[.N >= control_n_measurements], by = "IDAA"]

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

  # Make the wide dataset
  dat_wide <- data.table::dcast(
    data = dat,
    formula = IDAA + SCTyear_2010 + hirisk + ATG + VCMVPAT_pre +
      endpoint5_s + endpoint5 ~ .,
    fun = length
  )
  data.table::setnames(dat_wide, old = ".", new = "n_measurements")

  # Returns list of wide and long data
  datasets <- list("long" = dat, "wide" = dat_wide)
  return(datasets)
}

#' Fit model for one cell count
fit_indiv_JM <- function(datasets,
                         fform = ~ trans2 + trans3 + trans4 - 1) {

  # - Prep mstate
  JM_dat_wide <- data.table::copy(datasets$wide)
  event_names <- levels(JM_dat_wide$endpoint5_s)
  tmat <- mstate::trans.comprisk(
    K = 4,
    names = c("event_free", event_names[-1])
  )

  # Make numeric indicators
  ind_cols <- paste0("ind_", event_names[-1])
  JM_dat_wide[, (ind_cols) := lapply(
    event_names[-1], function(col) as.numeric(endpoint5_s == col)
  )]

  covs <- c("SCTyear_2010", "hirisk", "ATG", "VCMVPAT_pre")

  JM_msdat <- mstate::msprep(
    time = c(NA, rep("endpoint5", times = length(ind_cols))),
    status = c(NA, ind_cols),
    data = data.frame(JM_dat_wide),
    trans = tmat,
    keep = covs,
    id = "IDAA"
  )

  JM_msdat_expand <- mstate::expand.covs(
    JM_msdat,
    covs,
    append = TRUE,
    longnames = FALSE
  )

  # Cannot do left truncation
  coxCRfit <- survival::coxph(
    Surv(Tstop, status) ~
      SCTyear_2010.1 + hirisk.1 + # cellular intervention
      ATG.2 + hirisk.2 + # relapse
      ATG.3 + # NRF other; then no covars for gvhd (only current val in JM)
      strata(trans) + cluster(IDAA),
    data = JM_msdat_expand,
    x = TRUE,
    model = TRUE
  )

  lmeFit <- nlme::lme(
    CD8_abs_log ~ splines::ns(intSCT2_5, 3) * ATG + VCMVPAT_pre,
    #random = list(IDAA = pdDiag(form = ~ ns(intSCT2_2, 3))),
    random = ~ splines::ns(intSCT2_5, 3) | IDAA,
    control = nlme::lmeControl(opt = "optim"),
    data = datasets$long
  )

  # Prep functional form
  coxCRfit$model$trans2 <- as.numeric(coxCRfit$model$`strata(trans)` == "trans=2")
  coxCRfit$model$trans3 <- as.numeric(coxCRfit$model$`strata(trans)` == "trans=3")
  coxCRfit$model$trans4 <- as.numeric(coxCRfit$model$`strata(trans)` == "trans=4")

  # Run JM fit
  JMfit <- JM::jointModel(
    lmeObject = lmeFit,
    survObject = coxCRfit,
    timeVar = "intSCT2_5",
    method = "spline-PH-aGH",# or piecewise
    CompRisk = TRUE,
    interFact = list("value" = fform, data = coxCRfit$model),
    iter.EM = 200
  )

  return(JMfit)
}
