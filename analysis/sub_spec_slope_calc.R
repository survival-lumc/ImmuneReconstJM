# Editing code here https://github.com/drizopoulos/JM/blob/9d1fecf6c1fd4f6089d530929b5219b661134460/R/fitted.jointModel.R
#object <- preDLI_CD3__jointModel_both
object <- preDLI_CD3_jointModel_corr

# Here we collect terms regarding derivative wrt time
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

# Here we need to edit: repeat each line by 3 (# of comp risks)
#b <- b[rep(seq_len(nrow(b)), each = 3), ]
#Bullshit!!

# Now can continue with rest of code
id <- object$id

# This is pat-specific slopes
ff <- c(X.deriv %*% betas[indFixed] +
          rowSums(Z.deriv * b[id, indRandom, drop = FALSE]))
names(ff) <- names(object$y$y)
Xtime.deriv <- object$x$Xtime.deriv
Ztime.deriv <- object$x$Ztime.deriv

# This has len
ffEvent <- c(Xtime.deriv %*% betas[indFixed] +
               rowSums(Ztime.deriv * b[, indRandom, drop = FALSE]))
id <- object$id
idT <- object$x$idT


# Split here --------------------------------------------------------------



ffEvent2 <- unlist(mapply("c", split(ff, id), split(ffEvent, idT)),
                  use.names = FALSE)
times <- object$times
Time <- exp(object$y$logT)
ind <- unlist(mapply("c", split(times, id), split(Time, idT)),
              use.names = FALSE)
ffEvent2 <- ffEvent2[ind != 0]
names(ffEvent2) <- seq_along(ffEvent2)
ffEvent2
