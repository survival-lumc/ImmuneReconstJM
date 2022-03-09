tar_load(
  c(
    dat_merged,
    NMA_preDLI_datasets
  )
)

# Try now with cr_setup function (interaction parametrisation).
# - also bivariate CD4 and CD8?

jm_mod <- jm(
  Surv_object = preDLI_cox,
  Mixed_objects = list(preDLI_CD3_long_reffcorr),
  time_var = "intSCT2_5",
  functional_forms = list(
    "CD3_abs_log" = ~ value(CD3_abs_log) + value(CD3_abs_log):(strata(trans) - 1)
  ),
  id_var = "IDAA",
  n_iter = 13000L, n_burnin = 3000L
)

ggtraceplot(jm_mod, "alphas", gridrows = 1, gridcols = 3, grid = TRUE)
# Try with slope..?
ggtraceplot(jm_mod, "betas", gridrows = 6, gridcols = 3, grid = TRUE)


jm_slopes <- jm(
  Surv_object = preDLI_cox,
  Mixed_objects = list(preDLI_CD3_long_reffcorr),
  time_var = "intSCT2_5",
  functional_forms = list(
    #"CD3_abs_log" = ~ value(CD3_abs_log) + slope(CD3_abs_log) +
    #  (value(CD3_abs_log) + slope(CD3_abs_log)):(strata(trans) - 1)
    "CD3_abs_log" = ~ value(CD3_abs_log) + value(CD3_abs_log):(strata(trans) - 1)
  ),
  id_var = "IDAA",
  n_iter = 8000L,
  n_burnin = 3000L,
  #priors = list("penalty_alphas" = "horseshoe") #"horseshoe")
  # Numver of alphas
  priors = list(
    Tau_alphas = lapply(seq_len(3), function(alpha) matrix(data = 0.75))
  ),
  base_hazard_segments = 5
)

jm_slopes
gelman_diag(jm_slopes)
jm_slopes$control$diff
jm_slope
# Try les base_hazard_segments?? JM using has 7 bs gammas coefs
# try base_hazard segments = 5

# Try dynamic prediction

# Dyn pred ----------------------------------------------------------------

dat_long <- NMA_preDLI_datasets$long
dat_wide <- NMA_preDLI_datasets$wide

dat_cox <- {
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
  msdat_expand
}

ND_long <- as.data.frame(dat_long[dat_long$IDAA == "5707", ])#[1:6, ]
ND_long$status <- 0
ND_long$IDAA <- droplevels(ND_long$IDAA)
#preDLI_cox$model
#preDLI_cox
ND_event <- as.data.frame(dat_cox[dat_cox$IDAA == "5707", ])
ND_event$status <- 0
ND <- list(newdataL = ND_long, newdataE = ND_event)

predict(jm_slopes, newdata = ND, return_newdata = TRUE,
        times = seq(6.1, 7, length = 0.1))

predLong <- predict(jm_slopes, newdata = ND, return_newdata = TRUE,
                    times = seq(6.1, 7, length = 0.1))
predLong$newdata2

predEvent <- predict(jm_slopes, newdata = ND, return_newdata = TRUE,
                     process = "event")

# ... try using crisksetup

predEvent$pred_CIF
plot(predLong, predEvent)

ylim_long_outcome_range = FALSE,
     col_line_event = c("#03BF3D", "#FF0000"),
     fill_CI_event = c("#03BF3D4D", "#FF00004D"), pos_ylab_long = c(0.1, 20))
