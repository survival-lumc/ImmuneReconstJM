fitted_slopes_long <- function(object) {

  # Editing code here https://github.com/drizopoulos/JM/blob/9d1fecf6c1fd4f6089d530929b5219b661134460/R/fitted.jointModel.R

  # Here we collect terms regarding derivative w.r.t time
  derivForm <- object$derivForm
  indFixed <- derivForm$indFixed
  indRandom <- derivForm$indRandom
  TermsX.deriv <- object$termsYx.deriv
  TermsZ.deriv <- object$termsYz.deriv
  mfX.deriv <- model.frame(TermsX.deriv, data = object$data)
  mfZ.deriv <- model.frame(TermsZ.deriv, data = object$data)
  X.deriv <- model.matrix(derivForm$fixed, mfX.deriv)
  Z.deriv <- model.matrix(derivForm$random, mfZ.deriv)

  # Fixed and random effects
  betas <- object$coefficients$betas
  b <- ranef(object)

  # This is pat-specific slopes
  id <- object$id
  ff <- c(X.deriv %*% betas[indFixed] + rowSums(Z.deriv * b[id, indRandom, drop = FALSE]))
  names(ff) <- names(object$y$y)
  return(ff)
}

get_int_tangent <- function(x, y, slope) {
  m <- slope * x
  int <- -x * slope + y
  return(int)
}
