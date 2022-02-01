tar_load(
  c(
    dat_merged,
    NMA_preDLI_datasets,
    NMA_postDLI_datasets_CD3,
    #NMA_postDLI_datasets_CD4,
    #NMA_postDLI_datasets_CD8,
    preDLI_CD3__jointModel_both,
    #preDLI_CD4_jointModel_both,
    #preDLI_CD8_jointModel_both,
    postDLI_CD3_jointModel_both,#
    #postDLI_CD4_jointModel_both,
    #postDLI_CD8_jointModel_both
  )
)

dat_long <- NMA_postDLI_datasets_CD3$long

dat_long[sec_endpoint_s == "rel_nrf"] |>
  ggplot(aes(intSCT2_5_reset, CD3_abs_log, group = IDAA)) +
  geom_point() +
  facet_wrap(~ IDAA)


mod <- run_preDLI_longitudinal(
  cell_line = "CD3_abs_log",
  form_fixed = "ns(intSCT2_5_reset, 3) * earlylow_DLI + ATG", #"ns(intSCT2_5_reset, 3)", #"ns(intSCT2_5_reset, 3) * earlylow_DLI + ATG",
  form_random = ~ ns(intSCT2_5_reset, 3) | IDAA,
  #form_random = list(IDAA = pdDiag(~ ns(intSCT2_5_reset, 3))),
  dat = dat_long
)

summary(mod)

dat_long[, "preds" := predict(mod)]

dat_long[sec_endpoint_s == "cens"] |>
  ggplot(aes(intSCT2_5_reset, CD3_abs_log, group = IDAA)) +
  geom_point() +
  #geom_point(aes(y = preds), col = "blue") +
  geom_line(aes(y = preds), col = "blue") +
  facet_wrap(~ IDAA) +
  coord_cartesian(ylim = c(0, 8))

dat_long[sec_endpoint_s == "gvhd"] |>
  ggplot(aes(intSCT2_5_reset, CD3_abs_log, group = IDAA)) +
  geom_point() +
  #geom_point(aes(y = preds), col = "blue") +
  geom_line(aes(y = preds), col = "blue") +
  facet_wrap(~ IDAA) +
  coord_cartesian(ylim = c(0, 8))

dat_long[sec_endpoint_s == "rel_nrf"] |>
  ggplot(aes(intSCT2_5_reset, CD3_abs_log, group = IDAA)) +
  geom_point() +
  #geom_point(aes(y = preds), col = "blue") +
  geom_line(aes(y = preds), col = "blue") +
  facet_wrap(~ IDAA) +
  coord_cartesian(ylim = c(0, 8))

# Marginals
newdat_jm <- expand.grid(
  "ATG" = levels(dat_long$ATG),
  #"CMV_PD" = levels(dat_long$CMV_PD),
  "earlylow_DLI" = levels(dat_long$earlylow_DLI),
  "intSCT2_5_reset" = seq(0, 12, by = 0.05)
)

newdat_jm$preds_marg <- predict(mod, newdata = newdat_jm, level = 0L)

newdat_jm |>
  ggplot(aes(intSCT2_5_reset, preds_marg)) +
  geom_line(aes(group = interaction(ATG, earlylow_DLI),
                col = interaction(ATG, earlylow_DLI),
                linetype = interaction(ATG, earlylow_DLI)))  +
  coord_cartesian(ylim = c(0, 8))
