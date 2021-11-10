# Run code here, then incorporate into target pipeline..

lymphocytes <- data.table(readRDS("data-raw/2021-09-22_v7/lymphocytes.rds"))
variables <- data.table(readRDS("data-raw/2021-09-22_v7/variables.rds"))
dat_merged <- prepare_raw_data(lymphocytes, variables)


dat_wide <- data.table::dcast(
  data = dat_merged,
  formula = IDAA + hirisk + CMV_PD +
    endpoint6_s + endpoint6 + GvHD_prophylaxis + TCD + TCD2 + TCDmethod ~ .,
  fun = length
)

table(dat_wide$TCD)
table(dat_wide$endpoint6_s, dat_wide$TCD)

table(dat_wide$TCD, dat_wide$TCD2)
table(dat_wide$TCD)


# Testing  ----------------------------------------------------------------


tar_load(datasets)
dat_long <- datasets$long
dat_wide <- datasets$wide


table(dat_wide$endpoint6_s, dat_wide$uDLI_s)

tar_load(c(long_submodels, cox_all_dli, JM_CD4_allDLI_nointer))
long_submodels$CD4_abs_log |> summary() # try adding correlations
cox_all_dli

JM_CD4_allDLI_nointer |> summary()
JM_CD4_allDLI_nointer$coefficients$betas


newdat <- expand.grid(
  "ATG" = levels(dat_wide$ATG),
  "CMV_PD" = levels(dat_wide$CMV_PD),
  "intSCT2_5" = seq(0.1, 18, by = 0.1)
)

predict(
  JM_CD4_allDLI_nointer,
  newdat,
  type = "Marginal",
  interval = "confidence",
  returnData = TRUE
) |>
  ggplot(aes(intSCT2_5, pred, col = ATG, group = ATG)) +
  geom_ribbon(aes(ymin = low, ymax = upp), col = NA, fill = "grey70", alpha = 0.7) +
  geom_line(size = 1.25) +
  facet_wrap(~ CMV_PD) +
  theme_minimal()


# Try raw data plots
dat_long[, cr_ind := factor(
  x = ifelse(endpoint6_s != "cell_interv", as.character(endpoint6_s), "cens"),
  levels = c("cens", "REL", "NRF_gvhd", "NRF_other")
)]

dat_long |>
  ggplot(aes(intSCT2_5, CD4_abs_log, col = IDAA)) +
  geom_line(aes(group = IDAA), alpha = 0.9)+#, col = ATG)) +
  geom_point(aes(group = IDAA))+#, col = ATG)) +
  facet_grid(cr_ind ~ CMV_PD * ATG) + #  ~ IDAA
  coord_cartesian(xlim = c(0, 24)) +
  theme_bw() +
  theme(legend.position = "none")




dat_long |>
  ggplot(aes(intSCT2_5, CD4_abs_log)) +
  geom_line(aes(group = IDAA), alpha = 0.7, col = "lightgray") +#, col = ATG)) +
  geom_point(aes(group = IDAA), alpha = 0.7,  col = "lightgray") +#, col = ATG)) +
  geom_smooth(
    aes(group = ATG, col = ATG),
    method = "loess",
    formula = y ~ x,
    #col = "black",
    se = FALSE,
    span = 0.75,
    size = 1.25
  ) +
  facet_grid(cr_ind ~ CMV_PD) + #  ~ IDAA
  coord_cartesian(xlim = c(0, 24)) +
  theme_bw() #+
  #theme(legend.position = "none")



# Mixed models in the subsets ---------------------------------------------

endpoints <- c("cens", "REL", "NRF_gvhd", "NRF_other")
names(endpoints) <- endpoints

newdat <- expand.grid(
  "ATG" = levels(dat_wide$ATG),
  "CMV_PD" = levels(dat_wide$CMV_PD),
  "intSCT2_5" = seq(0.1, 18, by = 0.1)
)

long_subsets <- lapply(endpoints, function(ev) {

  long_manual <- lme(
    fixed = CD4_abs_log ~ ns(intSCT2_5, 3) * ATG + CMV_PD,
    random = list(IDAA = pdDiag(~ ns(intSCT2_5, 3))),
    data = dat_long[cr_ind == ev],
    control = nlme::lmeControl(opt = "optim")
  )

  cbind.data.frame(
    newdat,
    "preds" = predict(
      long_manual,
      newdata = newdat,
      level = 0L
    ),
    "subset" = ev
  )
})


ggplot(data = rbindlist(long_subsets), aes(intSCT2_5, preds)) +
  # geom_line(
  #   data = dat_long,
  #   aes(y = CD4_abs_log, col = ATG, group = IDAA),
  #   show.legend = FALSE,
  #   alpha = 0.5
  # ) +
  geom_line(aes(group = ATG, col = ATG), size = 2) +
  facet_grid(subset ~ CMV_PD) +
  theme_bw()



# Raw dat -----------------------------------------------------------------



dat_long[cr_ind == "REL"] |>
  ggplot(aes(intSCT2_5, CD4_abs_log, group = IDAA, col = ATG)) +
  geom_line() +
  geom_point() +
  facet_wrap(~ CMV_PD)



# Fit mixed model on whole data -------------------------------------------



long_all <- lme(
  fixed = CD4_abs_log ~ ns(intSCT2_5, 4) * ATG,# + CMV_PD,
  random = list(IDAA = pdDiag(~ ns(intSCT2_5, 4))),#~ ns(intSCT2_5, 4) | IDAA
  data = dat_long,
  control = nlme::lmeControl(opt = "optim")
)

newdat <- expand.grid(
  "ATG" = levels(dat_wide$ATG),
  "CMV_PD" = levels(dat_wide$CMV_PD),
  "intSCT2_5" = seq(0.1, 18, by = 0.1)
)

cbind.data.frame(
  newdat,
  "preds" = predict(
    long_all,
    newdata = newdat,
    level = 0L
  )
) |>
  ggplot(aes(intSCT2_5, preds, col = ATG)) +
  geom_line(size = 1.25) +
  facet_wrap(~ CMV_PD)



long_df2 <- lme(
  fixed = CD4_abs_log ~ ns(intSCT2_5, 2) * ATG + CMV_PD,
  random = list(IDAA = pdDiag(~ ns(intSCT2_5, 2))),#~ ns(intSCT2_5, 4) | IDAA
  data = dat_long,
  control = nlme::lmeControl(opt = "optim")
)

long_df3 <- lme(
  fixed = CD4_abs_log ~ ns(intSCT2_5, 3) * ATG + CMV_PD,
  random = list(IDAA = pdDiag(~ ns(intSCT2_5, 3))),#~ ns(intSCT2_5, 4) | IDAA
  data = dat_long,
  control = nlme::lmeControl(opt = "optim")
)

long_df4 <- lme(
  fixed = CD4_abs_log ~ ns(intSCT2_5, 4) * ATG + CMV_PD,
  random = list(IDAA = pdDiag(~ ns(intSCT2_5, 4))),#~ ns(intSCT2_5, 4) | IDAA
  data = dat_long,
  control = nlme::lmeControl(opt = "optim")
)

long_df5 <- lme(
  fixed = CD4_abs_log ~ ns(intSCT2_5, 5) * ATG + CMV_PD,
  random = list(IDAA = pdDiag(~ ns(intSCT2_5, 5))),#~ ns(intSCT2_5, 4) | IDAA
  data = dat_long,
  control = nlme::lmeControl(opt = "optim")
)

long_df6 <- lme(
  fixed = CD4_abs_log ~ ns(intSCT2_5, 6) * ATG + CMV_PD,
  random = list(IDAA = pdDiag(~ ns(intSCT2_5, 6))),#~ ns(intSCT2_5, 4) | IDAA
  data = dat_long,
  control = nlme::lmeControl(opt = "optim")
)

long_df7 <- lme(
  fixed = CD4_abs_log ~ ns(intSCT2_5, 7) * ATG + CMV_PD,
  random = list(IDAA = pdDiag(~ ns(intSCT2_5, 7))),#~ ns(intSCT2_5, 4) | IDAA
  data = dat_long,
  control = nlme::lmeControl(opt = "optim")
)

anova(
  update(long_df3, . ~ ., method = "ML"),
  update(long_df4, . ~ ., method = "ML")
)

anova(
  update(long_df3, . ~ ., method = "ML"),
  update(long_df4, . ~ ., method = "ML"),
  update(long_df5, . ~ ., method = "ML"),
  update(long_df6, . ~ ., method = "ML"),
  update(long_df7, . ~ ., method = "ML")
)



# Smoothing of whole thing ------------------------------------------------

melt.data.table(
  data = dat_long,
  id.vars = c("IDAA", "ATG", "CMV_PD", "intSCT2_5", "endpoint6_s"),
  measure.vars = paste0(c("CD4", "CD8", "CD19", "CD3", "NK"), "_abs_log"),
  variable.name = "cell_line",
  value.name = "cell_count"
) |>
  ggplot(aes(intSCT2_5, cell_count)) +
  geom_line(aes(group = IDAA), alpha = 0.5, col = "lightgray") +
  geom_smooth(
    method = "loess",
    se = FALSE,
    aes(group = ATG, col = ATG),
    size = 1.25,
    formula = y ~ x
  ) +
  facet_grid(cell_line ~ endpoint6_s) +
  scale_y_continuous(
    breaks = log(c(5, 25, 100, 500, 1500)),
    labels = c(5, 25, 100, 500, 1500)
  ) +
  theme_bw()



