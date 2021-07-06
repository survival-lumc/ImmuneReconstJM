# Load all objects
tar_load(
  c(
    datasets,
    dli_msdata,
    long_submodels,
    cox_submodel_all_dli,
    cox_submodel_no_dli,
    CD4_all_dli,
    CD4_all_dli_avfform,
    CD8_all_dli,
    CD8_all_dli_avfform,
    reference_values
  )
)

# Rename wide data
dat_wide <- datasets$wide

# For longitudinal predictions
newdat <- expand.grid(
  "ATG" = levels(dat_wide$ATG),
  "VCMVPAT_pre" = levels(dat_wide$VCMVPAT_pre),
  "intSCT2_5" = seq(0.1, 24, by = 0.1)
)

# Cox model ---------------------------------------------------------------

summary(cox_submodel_all_dli)


# CD4 JM: comparing cox ---------------------------------------------------


# Without interaction with DLI in association
summary(CD4_all_dli_avfform)
CD4_all_dli_avfform$coefficients$gammas
cox_submodel_all_dli$coefficients
CD4_all_dli_avfform$coefficients$alpha

# With interaction - DLI.2 explodes? Check with sim data
summary(CD4_all_dli)
CD4_all_dli$coefficients$gammas
cox_submodel_all_dli$coefficients
CD4_all_dli$coefficients$alpha

confint(CD4_all_dli)

# Joint model just for gvhd -----------------------------------------------



mod_gvhd <- coxph(Surv(Tstart, Tstop, status) ~
                    DLI + #hirisk + SCT_May2010 + VCMVPAT_pre +
                    ATG, data = dli_msdata,
      subset = (trans == 2), x = TRUE, model = TRUE, cluster = IDAA)

#mod_gvhd$model$DLI <- dli_msdata[dli_msdata$trans == 2, ]$DLI

jm_mod_gvhd <- jointModel(
  lmeObject = long_submodels$CD4_abs_log,
  survObject = mod_gvhd,
  CompRisk = FALSE,
  timeVar = "intSCT2_5",
  method = "spline-PH-aGH",
  interFact = list(
    "value" = ~ DLI
  ),
  iter.EM = 200 # see p.68 rizop book
)

summary(jm_mod_gvhd)
mod_gvhd$coefficients
jm_mod_gvhd$coefficients$gammas
jm_mod_gvhd$coefficients$alpha


jm_mod_gvhd_nointer <- jointModel(
  lmeObject = long_submodels$CD4_abs_log,
  survObject = mod_gvhd,
  CompRisk = FALSE,
  timeVar = "intSCT2_5",
  method = "spline-PH-aGH",
  interFact = NULL,
  iter.EM = 200 # see p.68 rizop book
)

mod_gvhd$coefficients
jm_mod_gvhd_nointer$coefficients$gammas
jm_mod_gvhd_nointer$coefficients$alpha


# CD4 JM: comparing long --------------------------------------------------

long_submodels$CD4_abs_log$coefficients$fixed
CD4_all_dli$coefficients$betas



# Concerns on event number pre-post DLI? ----------------------------------

table(dat_wide$uDLI_s, dat_wide$endpoint6_s)


# CD4 longitudinal --------------------------------------------------------


JM::predict.jointModel(
  CD4_all_dli,
  newdata = newdat,
  type = "Marginal",
  idVar = "IDAA",
  returnData = TRUE,
  interval = "confidence"
) %>%
  ggplot(aes(x = intSCT2_5, y = pred)) +
  geom_ribbon(
    aes(ymin = low, ymax = upp, group = ATG, col = ATG),
    fill = "gray",
    alpha = 0.5,
    col = NA
  ) +
  geom_line(aes(group = ATG, col = ATG), size = 1.5, alpha = 0.75) +
  #geom_text(aes(label = label), hjust = 0, na.rm = TRUE, fontface = "bold") +
  facet_wrap(facets = ~ VCMVPAT_pre) +
  labs(x = "Time since alloHCT (months)", y = "CD4 cell counts") +
  # scale_y_continuous(
  #   breaks = log(c(5, 25, 100, 500, 1500)),
  #   labels = c(5, 25, 100, 500, 1500)
  # ) +
  theme_minimal() +
  #theme(legend.position = "none") +
  #coord_cartesian(xlim = c(0, 14)) +
  scale_color_brewer(palette = "Paired") +
  geom_vline(xintercept = 6, linetype = "dotted")


JM::predict.jointModel(
  CD4_all_dli,
  newdata = newdat,
  type = "Marginal",
  idVar = "IDAA",
  returnData = TRUE,
  interval = "confidence"
) %>%
  ggplot(aes(x = intSCT2_5, y = pred)) +
  geom_ribbon(
    aes(ymin = low, ymax = upp, group = ATG, col = ATG),
    fill = "gray",
    alpha = 0.5,
    col = NA
  ) +
  geom_line(aes(group = ATG, col = ATG), size = 1.5, alpha = 0.75) +
  #geom_text(aes(label = label), hjust = 0, na.rm = TRUE, fontface = "bold") +
  facet_wrap(facets = ~ VCMVPAT_pre) +
  labs(x = "Time since alloHCT (months)", y = "CD4 cell counts") +
   scale_y_continuous(
     breaks = log(c(5, 25, 100, 500, 1500)),
     labels = c(5, 25, 100, 500, 1500)
   ) +
  theme_minimal() +
  #theme(legend.position = "none") +
  #coord_cartesian(xlim = c(0, 14)) +
  scale_color_brewer(palette = "Paired")


# CD8 longitudinal --------------------------------------------------------


JM::predict.jointModel(
  CD8_all_dli,
  newdata = newdat,
  type = "Marginal",
  idVar = "IDAA",
  returnData = TRUE,
  interval = "confidence"
) %>%
  ggplot(aes(x = intSCT2_5, y = pred)) +
  geom_ribbon(
    aes(ymin = low, ymax = upp, group = ATG, col = ATG),
    fill = "gray",
    alpha = 0.5,
    col = NA
  ) +
  geom_line(aes(group = ATG, col = ATG), size = 1.5, alpha = 0.75) +
  #geom_text(aes(label = label), hjust = 0, na.rm = TRUE, fontface = "bold") +
  facet_wrap(facets = ~ VCMVPAT_pre) +
  labs(x = "Time since alloHCT (months)", y = "CD8 cell counts") +
  scale_y_continuous(
    breaks = log(c(5, 25, 100, 500, 1500)),
    labels = c(5, 25, 100, 500, 1500)
  ) +
  theme_minimal() +
  #theme(legend.position = "none") +
  #coord_cartesian(xlim = c(0, 14)) +
  scale_color_brewer(palette = "Paired")


# CD8 cox -----------------------------------------------------------------




# Without interaction with DLI in association
summary(CD8_all_dli_avfform)
summary(CD8_all_dli)

CD4_all_dli_avfform$coefficients$gammas
cox_submodel_all_dli$coefficients
CD4_all_dli_avfform$coefficients$alpha

# With interaction - DLI.2 explodes? Check with sim data
summary(CD4_all_dli)
CD4_all_dli$coefficients$gammas
cox_submodel_all_dli$coefficients
CD4_all_dli$coefficients$alpha
