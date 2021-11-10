# Jm vs jmbayes 2
# https://github.com/graemeleehickey/comprisk/blob/7e558447541f8734d94135cd2e509528dd53d650/rizopoulos2012.R
# Check with https://github.com/drizopoulos/JMbayes2/blob/29a4f72e192d7847686432501f6116b381e32764/Development/Dev_Local_GP/MS_CR/Competing_Risks_Reproduce.R


epileptic <- read.table("epileptic.txt", header = TRUE)
head(epileptic)

epileptic$time <- epileptic$time / 365.25
epileptic$with.time <- epileptic$with.time / 365.25
epileptic$treat <- as.numeric(epileptic$treat == "LTG")
epileptic <- epileptic[epileptic$time < epileptic$with.time, ]
epileptic$id <- factor(epileptic$id)


# Time-to-event data (one row per subject)
survdat <- epileptic
survdat <- survdat[!duplicated(survdat$id), ]
survdat$status <- rep("alive", nrow(survdat))
survdat$status[survdat$with.status.uae == 1] <- "UAE"
survdat$status[survdat$with.status.isc == 1] <- "ISC"
survdat <- survdat[ , c(1, 4, 8, 12)]
survdat.CR <- crLong(survdat, statusVar = "status",
                     censLevel = "alive", nameStrata = "CR")


# Mods
lmeFit <- lme(dose ~ treat * time,
              random = list(id = pdDiag(~ time)),
              #random = ~ time | id,
              data = epileptic) # check na.action
summary(lmeFit)

# Time-to-event submodel
coxFit <- coxph(Surv(with.time, status2) ~ treat * strata(CR),
                data = survdat.CR,
                x = TRUE,
                model = TRUE,
                cluster = id)

summary(coxFit)

dform <- list(fixed = ~ 1 + treat, indFixed = 3:4, random = ~ 1, indRandom = 2)

jointFit2 <- jointModel(
  lmeObject = lmeFit,
  survObject = coxFit,
  timeVar = "time",
  method = "spline-PH-aGH",
  CompRisk = TRUE,
  interFact = list(slope = ~ strata(CR)),
  #GHk = 9,
  #GKk = 15,
  #iter.qN = 500,
  #numeriDeriv = "cd",
  parameterization = "slope",
  derivForm = dform,
  verbose = FALSE
)

summary(jointFit2)

jmfit2_alt <- jm(
  Surv_object = coxFit,
  Mixed_objects = lmeFit,
  #data_Surv = survdat.CR,
  time_var = "time",
  #n_iter = 200L,
  #n_burnin = 100L,
  id_var = "id",
  functional_forms = ~ slope(dose):strata(CR)
)

coefficients(jmfit2)$association
jointFit2_alt$coefficients$Dalpha

coefficients(jmfit2)$gammas
jointFit2_alt$coefficients$gammas

summary(jointFit2_alt)
jmfit2_alt


# Try with msprep ---------------------------------------------------------


tmat <- trans.comprisk(K = 2, c("alive", "ISC", "UAE"))
survdat$ISC_ind <- as.numeric(survdat$status == "ISC")
survdat$UAE_ind <- as.numeric(survdat$status == "UAE")

msdat <- msprep(
  time = c(NA, "with.time", "with.time"),
  status = c(NA, "ISC_ind", "UAE_ind"),
  id = "id",
  data = survdat,
  trans = tmat,
  keep = "treat"
)

msdat_exp <- expand.covs(
  msdat, covs = "treat", longnames = FALSE
)

coxFit_alt <- coxph(
  Surv(time, status) ~ treat.1 + treat.2 + strata(trans),
  data = msdat_exp,
  x = TRUE,
  model = TRUE,
  cluster = id
)

coef(coxFit_alt)
coef(coxFit_alt)


jointFit2_alt <- jointModel(
  lmeObject = lmeFit,
  survObject = coxFit_alt,
  timeVar = "time",
  method = "spline-PH-aGH",
  CompRisk = TRUE,
  interFact = list(slope = ~ strata(trans) - 1),
  #GHk = 9,
  #GKk = 15,
  #iter.qN = 500,
  #numeriDeriv = "cd",
  parameterization = "slope",
  derivForm = dform,
  verbose = FALSE
)

summary(jointFit2_alt)
jointFit2$coefficients$gammas
jointFit2_alt$coefficients$gammas

jointFit2$coefficients$Dalpha
jointFit2_alt$coefficients$Dalpha

msdat_exp$trans1 <- as.numeric(msdat_exp$trans == 1)
msdat_exp$trans2 <- as.numeric(msdat_exp$trans == 2)


jmfit2_alt <- jm(
  Surv_object = coxFit_alt,
  Mixed_objects = lmeFit,
  #data_Surv = msdat_exp,
  time_var = "time",
  n_iter = 2000L,
  n_burnin = 1000L,
  id_var = "id",
  # This is the contrast form:
  #functional_forms = ~ slope(dose) + slope(dose):(strata(trans) - 1)
)

coefficients(jmfit2_alt)$association
jointFit2_alt$coefficients$Dalpha

coefficients(jmfit2_alt)$gammas
jointFit2_alt$coefficients$gammas

summary(jointFit2)
jmfit2_alt

# Refit models
