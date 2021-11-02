tar_load(c(datasets, cox_all_dli))

long <- run_longitudinal_submodels(
  datasets$long,
  which_cells = c("CD3_abs_log"),
  df_splines = 5,
  ranef_structure = "diagonal"
)

source("data-raw/prepare-raw-data.R")
source("R/submodel-wrappers.R")
source("R/joint-model-wrappers.R")

mod_simplified_cd3 <- run_jointModel(
  long_obj = long$CD3_abs_log,
  surv_obj = cox_all_dli,
  fform = ~ trans1 + trans2 + trans3 - 1
)

summary(mod_simplified_cd3)


long_manual <- lme(
  fixed = CD3_abs_log ~ ns(intSCT2_5, 3) * ATG + VCMVPAT_pre,
  random = list(IDAA = pdDiag(~ ns(intSCT2_5, 3))),
  data = datasets$long,
  control = nlme::lmeControl(opt = "optim")
)


# Predict just with mixed model
cbind.data.frame(
  newdat,
  "preds" = predict(
    long_manual,
    newdata = newdat,
    level = 0L
  )
) |>
  ggplot(aes(intSCT2_5, preds)) +
  geom_line(
    data = dat_long,
    aes(y = CD3_abs_log, col = ATG, group = IDAA),
    show.legend = FALSE,
    alpha = 0.7
  ) +
  geom_line(aes(group = ATG, col = ATG), size = 2) +
  facet_wrap(~ VCMVPAT_pre) +
  theme_bw()



dat_wide <- datasets$wide

newdat <- expand.grid(
  "ATG" = levels(dat_wide$ATG),
  "VCMVPAT_pre" = levels(dat_wide$VCMVPAT_pre),
  "intSCT2_5" = seq(0.1, 24, by = 0.1)
)

predict(
  mod_simplified_cd3,
  newdat,
  type = "Marginal",
  interval = "confidence",
  returnData = TRUE
) |>
  ggplot(aes(intSCT2_5, pred, col = ATG, group = ATG)) +
  geom_ribbon(aes(ymin = low, ymax = upp),
              col = "lightgray", fill = "lightgray",
              alpha = 0.7) +
  geom_line(size = 1.25) +
  facet_wrap(~ VCMVPAT_pre) +
  theme_minimal()

dat_long <- datasets$long

dat_long[, cr_ind := factor(
  x = ifelse(endpoint6_s != "cell_interv", as.character(endpoint6_s), "cens"),
  levels = c("cens", "REL", "NRF_gvhd", "NRF_other")
)]

dat_long |>
  ggplot(aes(intSCT2_5, CD3_abs_log)) +
  geom_line(aes(group = IDAA, col = IDAA)) +
  theme(legend.position = "none")

dat_long |>
  ggplot(aes(intSCT2_5, CD3_abs_log)) +
  geom_line(aes(group = IDAA, col = IDAA)) +
  theme(legend.position = "none") +
  facet_grid(ATG * VCMVPAT_pre ~ cr_ind) +
  coord_cartesian(xlim = c(0, 12))



dat_long |>
  ggplot(aes(intSCT2_5, CD3_abs_log)) +
  geom_line(aes(group = IDAA, col = ATG)) +
  geom_point(aes(group = IDAA, col = ATG)) +
  facet_grid(cr_ind ~ VCMVPAT_pre) +
  coord_cartesian(xlim = c(0, 12))


# Individual traj ---------------------------------------------------------

# Over fitting on patients?
dat_long |>
  ggplot(aes(intSCT2_5, CD19_abs_log, col = IDAA)) +
  geom_line(aes(group = IDAA), alpha = 0.9)+#, col = ATG)) +
  geom_point(aes(group = IDAA))+#, col = ATG)) +
  facet_grid(cr_ind ~ VCMVPAT_pre * ATG) + #  ~ IDAA
  coord_cartesian(xlim = c(0, 12)) +
  theme_bw() +
  theme(legend.position = "none")


dat_long[cr_ind == "NRF_gvhd"] |>
  ggplot(aes(intSCT2_5, CD3_abs_log, col = interaction(ATG, VCMVPAT_pre))) +
  geom_line(aes(group = IDAA), alpha = 0.9)+#, col = ATG)) +
  geom_point(aes(group = IDAA))+#, col = ATG)) +
  facet_wrap(~ IDAA) + #  ~ IDAA
  coord_cartesian(xlim = c(0, 12)) +
  theme(legend.position = "top")


# Extra -------------------------------------------------------------------

tar_load(c(cox_all_dli, long_submodels, dli_msdata))

cox_all_dli$model$trans1 <- as.numeric(cox_all_dli$model$`strata(trans)` == "trans=1")
cox_all_dli$model$trans2 <- as.numeric(cox_all_dli$model$`strata(trans)` == "trans=2")
cox_all_dli$model$trans3 <- as.numeric(cox_all_dli$model$`strata(trans)` == "trans=3")

run_jointModel(
  long_obj = long_submodels$CD8_abs_log,
  surv_obj = cox_all_dli,
  fform = ~ trans1 + trans2 + trans3 - 1
)

jointModel(
  long_submodels$CD4_abs_log,
  cox_all_dli,
  timeVar = "intSCT2_5",
  interFact = list("value" = ~ trans1 + trans2 + trans3 - 1),
  method = "spline-PH-aGH"
)

# Check shtuff
tar_load(c(JM_CD8_allDLI_nointer, JM_CD4_allDLI_nointer, JM_CD3_allDLI_nointer))
JM_CD8_allDLI_nointer$coefficients$gammas
JM_CD4_allDLI_nointer$coefficients$gammas
JM_CD3_allDLI_nointer$coefficients$gammas

JM_CD3_allDLI_nointer$convergence

summary(JM_CD3_allDLI_nointer)
# Increasement EM iterations for CD3...


tar_load(datasets)
dat_wide <- datasets$wide

newdat <- expand.grid(
  "ATG" = levels(dat_wide$ATG),
  "VCMVPAT_pre" = levels(dat_wide$VCMVPAT_pre),
  "intSCT2_5" = seq(0.1, 24, by = 0.1)
)

predict(
  JM_CD3_allDLI_nointer,
  newdat,
  type = "Marginal",
  interval = "confidence",
  returnData = TRUE
) |>
  ggplot(aes(intSCT2_5, pred, col = ATG, group = ATG)) +
  geom_ribbon(aes(ymin = low, ymax = upp), col = "grey70", fill = "grey70") +
  geom_line(size = 1.25) +
  facet_wrap(~ VCMVPAT_pre) +
  theme_minimal()


# Check CD3 trajectories --------------------------------------------------



