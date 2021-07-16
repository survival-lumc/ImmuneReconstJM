

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
