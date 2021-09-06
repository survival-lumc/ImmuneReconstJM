# Meeting 5/08/2021

tar_load(
  c(
    datasets,
    dli_msdata,
    long_submodels,
    reference_values,
    multivar_allDLI_nointer,
    multivar_allDLI_inter_penal,
    JM_CD4_allDLI_nointer,
    JM_CD4_allDLI_inter,
    JM_CD8_allDLI_nointer
  )
)



# Longitudinal trajectories -----------------------------------------------



dat_wide <- datasets$wide
newdat <- expand.grid(
  "ATG" = levels(dat_wide$ATG),
  "VCMVPAT_pre" = levels(dat_wide$VCMVPAT_pre),
  "intSCT2_5" = seq(0.1, 24, by = 0.1)
)

# Get model frames
mod <- multivar_allDLI_nointer
terms_FE <- mod$model_info$terms$terms_FE_noResp
model_frames_FE <- lapply(terms_FE, model.frame.default, data = newdat)
modmats_FE <- mapply(
  model.matrix.default,
  object = terms_FE,
  data = model_frames_FE,
  SIMPLIFY = FALSE
)

mcmc_elements <- names(mod$mcmc)
betas_mcmc <- as.list(mod$mcmc[grep("^betas", mcmc_elements, value = TRUE)])
n_samp_posterior <- 2000
betas <- lapply(betas_mcmc, function(b) {
  df_samps <- runjags::combine.mcmc(b)
  df_samps[sample(nrow(df_samps), size = n_samp_posterior, replace = FALSE), ]
})

preds_ls <- mapply(
  function(modmat, beta) {
    preds_df <- data.table(modmat %*% t(beta))

    res <- preds_df[, ':=' (
      pred = rowMeans(.SD),
      lower = apply(.SD, 1, quantile, probs = 0.025),
      upper = apply(.SD, 1, quantile, probs = 0.975)
    )][, c("pred", "lower", "upper")]

    cbind(newdat, res)
  },
  modmat = modmats_FE,
  beta = betas,
  SIMPLIFY = FALSE
)

df_preds <- rbindlist(preds_ls, idcol = "cell_type")

ggplot(df_preds, aes(intSCT2_5, pred, group = ATG)) +
  geom_ribbon(aes(ymin = lower, ymax = upper), fill = "lightgray", alpha = 0.7) +
  geom_line(aes(col = ATG), size = 1.25) +
  facet_grid(cell_type ~ VCMVPAT_pre) +
  labs(x = "Time since alloHCT (months)", y = "Cell counts") +
  scale_y_continuous(
    breaks = log(c(5, 25, 100, 500, 1500)),
    labels = c(5, 25, 100, 500, 1500)
  ) +
  theme_minimal()



# Exploring coefficients --------------------------------------------------


multivar_allDLI_nointer
coef(multivar_allDLI_nointer)
multivar_allDLI_inter_penal

tar_load(multivar_allDLI_slopeGVHD)
coef(multivar_allDLI_slopeGVHD)
coef(multivar_allDLI_nointer)
