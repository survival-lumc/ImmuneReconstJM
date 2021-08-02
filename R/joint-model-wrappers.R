# Helpers for JM::jointModel() --------------------------------------------


run_jointModel <- function(long_obj,
                           surv_obj,
                           fform,
                           ...) {

  # Set up model matrices of assoc. parameters
  surv_obj$model$trans1 <- as.numeric(surv_obj$model$`strata(trans)` == "trans=1")
  surv_obj$model$trans2 <- as.numeric(surv_obj$model$`strata(trans)` == "trans=2")
  surv_obj$model$trans3 <- as.numeric(surv_obj$model$`strata(trans)` == "trans=3")

  # Fit jointmodel
  JM_mod <- jointModel(
    lmeObject = long_obj,
    survObject = surv_obj,
    CompRisk = TRUE,
    timeVar = "intSCT2_5",
    method = "spline-PH-aGH",
    interFact = list(
      "value" = fform
    ),
    iter.EM = 200,
    ...
  )

  return(JM_mod)
}


# Helpers for JMbayes2::jm() ----------------------------------------------


