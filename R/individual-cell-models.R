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
fit_indiv_JM <- function(long_outcome,
                         datasets,
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

  form <- stats::reformulate(
    response = long_outcome,
    termlabels = "splines::ns(intSCT2_5, 3) * ATG + VCMVPAT_pre"
  )

  # Need to use do.call to evaluate arguments, otherwise this fails
  lmeFit <- do.call(
    nlme::lme,
    list(
      "fixed" = form,
      #random = list(IDAA = pdDiag(form = ~ ns(intSCT2_2, 3))),
      "random" = ~ splines::ns(intSCT2_5, 3) | IDAA,
      "control" = nlme::lmeControl(opt = "optim"),
      "data" = datasets$long
    )
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

summarise_event_submodel <- function(JMfit,
                                     digits = 3) {

  summary_raw <- data.table::data.table(
    round(summary(JMfit)$`CoefTable-Event`, digits),
    keep.rownames = TRUE
  )

  # Exclude baseline hazard spline coefficients
  data.table::setnames(summary_raw, "rn", "Coefficient")
  dt <- summary_raw[grep("^bs", x = Coefficient, invert = TRUE)]

  # Set up
  dt[, event_num := regmatches(
    x = Coefficient,
    m = regexpr(pattern = "[+1-9]$", text = Coefficient)
  )]

  dt[, event := factor(
    x = event_num,
    levels = seq_len(4),
    labels = c("Cell. intervention", "Relapse", "NRF: Other", "NRF: GVHD")
  )]

  dt[, Coefficient := gsub(x = Coefficient, pattern = "\\:.*", replacement = "")]
  dt[, Coefficient := gsub(x = Coefficient, pattern = "\\.[+1-9]$", replacement = "")]
  dt[, Coefficient := factor(
    x = Coefficient,
    levels = c("ATG", "hirisk", "SCTyear_2010", "Assoct"),
    labels = c("ATG", "Disease risk: high", "SCT pre-2010", "Assoct.")
  )]

  data.table::setorder(dt, event, Coefficient)
  dt_edit <- dt[, c("event", "Coefficient", "Value", "Std.Err", "p-value")]
  data.table::setnames(dt_edit, new = c("Event", "Coefficient", "log(HR)", "SE", "p-value"))
  dt_edit[, "HR [95% CI]" := paste0(
    as.character(round(exp(`log(HR)`), 3)), " [",
    as.character(round(exp(`log(HR)` - pnorm(0.975) * SE), 3)), ";",
    as.character(round(exp(`log(HR)` + pnorm(0.975) * SE), 3)), "]"
  )]
  dt_edit[, Event := NULL]

  return(dt_edit)
  # dt_edit %>%
  #   kableExtra::kbl(format = "html") %>%
  #   kable_paper("striped", full_width = F) %>%
  #   pack_rows("Cell. intervention", 1, 3) %>%
  #   pack_rows("Relapse", 4, 6) %>%
  #   pack_rows("NRF: Other", 7, 7) %>%
  #   pack_rows("NRF: GVHD", 8, 8)
}
