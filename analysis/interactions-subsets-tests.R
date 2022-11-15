# Add this to seprate file in R/ (UD sub analysis)
tar_load(c(NMA_preDLI_datasets, NMA_postDLI_datasets, preDLI_cox, preDLI_long_corr_CD4))

# Try preDLI subsets
UD_preDLI_long <- NMA_preDLI_datasets$long[ATG == "UD(+ATG)"]
UD_preDLI_wide <- NMA_preDLI_datasets$wide[ATG == "UD(+ATG)"]

# Also look at RD
RD_preDLI_long <- NMA_preDLI_datasets$long[ATG == "UD"]
RD_preDLI_wide <- NMA_preDLI_datasets$wide[ATG == "UD"]

table(UD_preDLI_wide$endpoint7_s)
table(UD_preDLI_wide$endpoint7_s, UD_preDLI_wide$hirisk) # All hirisk peope get gvhd

# Too few in RD
table(RD_preDLI_wide$endpoint7_s)
table(RD_preDLI_wide$endpoint7_s, RD_preDLI_wide$hirisk)

# Look more generally at endpoints with ATG
table(NMA_preDLI_datasets$wide$endpoint7_s, NMA_preDLI_datasets$wide$ATG)
table(NMA_preDLI_datasets$wide$endpoint7_s, NMA_preDLI_datasets$wide$hirisk)


# Start with CD4s
UD_cox <- run_preDLI_cox(
  form = Surv(time, status) ~
    #hirisk.1 + # GVHD
    hirisk.2 + # Relapse
    strata(trans), # NRF other
  dat_wide = UD_preDLI_wide
)

UD_CD4_long <- run_preDLI_longitudinal(
  cell_line = "CD4_abs_log",
  form_fixed = "ns(intSCT2_7, 3) * hirisk + CMV_PD",
  form_random = ~ 0 + ns(intSCT2_7, 3) | IDAA,
  #form_random = list(IDAA = pdDiag(~ 0 + ns(intSCT2_7, 3))),
  dat = UD_preDLI_long
)

VarCorr(UD_CD4_long)
summary(UD_CD4_long)

UD_CD4_JM_value <- jointModel(
  lmeObject = UD_CD4_long,
  survObject = UD_cox,
  CompRisk = TRUE,
  method = "spline-PH-aGH",
  timeVar = "intSCT2_7",
  parameterization = "value",
  interFact = list("value" = ~ strata(trans) - 1),
  iter.EM = 250,
  iter.qN = 1000,
  lng.in.kn = 3L,
  numeriDeriv = "cd",
  eps.Hes = 1e-04,
  verbose = TRUE
)

summary(tar_read(preDLI_JM_value_corr_CD4))$`CoefTable-Event`
summary(UD_CD4_JM_value)$`CoefTable-Event`

cbind(
  "Orig" = c(summary(tar_read(preDLI_JM_value_corr_CD4))$`CoefTable-Event`[1:8, "Value"], 0),
  "Inter" = summary(mod_interaction)$`CoefTable-Event`[1:9, "Value"]
)
#summary(tar_read(preDLI_JM_value_corr_CD4))$`CoefTable-Event`


summary(tar_read(preDLI_JM_value_corr_CD4_inter))$`CoefTable-Event`[1:11, ]

UD_CD4_JM_value$Hessian |> matrixcalc::is.positive.definite()
summary(UD_CD4_JM_value)

newdat_UD_preDLI <- expand.grid(
  "CMV_PD" = levels(UD_preDLI_long$CMV_PD),
  "hirisk" = levels(UD_preDLI_long$hirisk),
  "intSCT2_7" = seq(0, 6, by = 0.02)
)

predict(
  UD_CD4_JM_value,
  newdata = newdat_UD_preDLI,
  type = "Marginal",
  idVar = "IDAA",
  returnData = TRUE,
  interval = "confidence"
) |>
  ggplot(aes(x = intSCT2_7, y = pred, group = hirisk)) +
  geom_ribbon(aes(ymin = low, ymax = upp), fill = "gray", alpha = 0.5, col = NA) +
  geom_line(aes(linetype = hirisk, col = hirisk), size = 1.5, alpha = 0.75) +
  # Add repel labels?
  facet_wrap(facets = ~ CMV_PD) +
  labs(
    x = "Time since alloHCT (months)",
    y = "Cell counts"
  ) +
  scale_y_continuous(
    breaks = log(c(1, 5, 25, 100, 500, 1500)),
    labels = c(1, 5, 25, 100, 500, 1500)
  ) +
  coord_cartesian(xlim = c(0, 6), ylim = c(log(0.1), log(1500))) + # add a common ylim for all
  theme(legend.position = "bottom") +
  theme_bw(base_size = 14)



# Actually try interaction analysis ---------------------------------------

# Hard code transitions first
preDLI_cox$model$trans1 <- as.numeric(preDLI_cox$model$`strata(trans)` == "trans=1")
preDLI_cox$model$trans2 <- as.numeric(preDLI_cox$model$`strata(trans)` == "trans=2")
preDLI_cox$model$trans3 <- as.numeric(preDLI_cox$model$`strata(trans)` == "trans=3")

# Check the actual one that did run
summary(tar_read(preDLI_JM_value_corr_CD4))

mod_interaction <- jointModel(
  lmeObject = preDLI_long_corr_CD4,
  survObject = preDLI_cox,
  CompRisk = TRUE,
  method = "spline-PH-aGH",
  timeVar = "intSCT2_7",
  parameterization = "value",
  interFact = list(
    "value" = ~ trans1 + trans1:ATG.1 + trans2 + trans3 - 1
  ),
  iter.qN = 1000,
  lng.in.kn = 3L,
  numeriDeriv = "cd",
  eps.Hes = 1e-04,
  verbose = TRUE
)

cbind(
  "Orig" = c(summary(tar_read(preDLI_JM_value_corr_CD4))$`CoefTable-Event`[1:8, "Value"], 0),
  "Inter" = summary(mod_interaction)$`CoefTable-Event`[1:9, "Value"]
)

summary(mod_interaction)
table(NMA_preDLI_datasets$wide$endpoint7_s, NMA_preDLI_datasets$wide$ATG)

table(
  NMA_preDLI_datasets$wide$endpoint7_s,
  NMA_preDLI_datasets$wide$ATG
) #outcome vs UD

table(
  NMA_preDLI_datasets$wide$endpoint7_s,
  NMA_preDLI_datasets$wide$hirisk
)

summary(tar_read(preDLI_JM_value_corr_CD4))$`CoefTable-Event`[1:8, ]
summary(mod_interaction)$`CoefTable-Event`[1:10, ]

# Post-DLI ----------------------------------------------------------------


UD_postDLI_long <- NMA_postDLI_datasets$long[ATG == "UD(+ATG)"]
UD_postDLI_wide <- NMA_postDLI_datasets$wide[ATG == "UD(+ATG)"]


tar_load(c(postDLI_cox, postDLI_long_corr_CD4, postDLI_JM_corr_CD4))

postDLI_cox$model$trans1 <- as.numeric(postDLI_cox$model$`strata(trans)` == "trans=1")
postDLI_cox$model$trans2 <- as.numeric(postDLI_cox$model$`strata(trans)` == "trans=2")

postDLI_inter_CD4 <- jointModel(
  lmeObject = postDLI_long_corr_CD4,
  survObject = postDLI_cox,
  CompRisk = TRUE,
  parameterization = "value",
  interFact = list("value" = ~ trans1 + trans1:ATG.1 + trans2 + - 1),
  method = "spline-PH-aGH",
  timeVar = "intDLI1",
  lng.in.kn = 2L,
  iter.EM = 1000,
  #tol1 = 1e-6,
  #tol2 = 1e-6,
  #tol3 = .Machine$double.eps,
  #iter.qN = 1000,
  #numeriDeriv = "cd",
  #eps.Hes = 1e-04,
  verbose = TRUE
)

summary(postDLI_inter_CD4)
tar_load(NMA_postDLI_datasets)

# Post DLI not enough for contrast
table(NMA_postDLI_datasets$wide$ATG, NMA_postDLI_datasets$wide$sec_endpoint2_s)
