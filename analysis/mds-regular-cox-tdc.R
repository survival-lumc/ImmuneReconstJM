tar_load(
  c(
    datasets,
    dli_msdata,
    long_submodels,
    JM_CD4_allDLI_nointer,
    multivar_allDLI_nointer,
    multivar_allDLI_inter_penal
  )
)

summary(JM_CD4_allDLI_nointer)
dat_wide <- datasets$wide
dat_long <- datasets$long
dat_wide[, gvhd_ind := as.numeric(endpoint6_s == "NRF_gvhd")]

df1 <- tmerge(
  data1 = data.frame(dat_wide),
  data2 = data.frame(dat_wide),
  id = IDAA,
  status = event(endpoint6, gvhd_ind),
  dli = tdc(uDLI)
)

df2 <- tmerge(
  data1 = df1,
  data2 = data.frame(dat_long),
  id = IDAA,
  cd4 = tdc(intSCT2_5, CD4_abs_log),
  cd8 = tdc(intSCT2_5, CD8_abs_log)
)

coxph(
  Surv(tstart, tstop, status) ~ dli + ATG + cd4,
  data = data.frame(df2),
  model = TRUE,
  x = TRUE,
  cluster = IDAA
)

JM_CD4_allDLI_nointer$coefficients$gammas
JM_CD4_allDLI_nointer$coefficients$alpha


coxph(
  Surv(tstart, tstop, status) ~ dli * cd4 + ATG,
  data = data.frame(df2),
  model = TRUE,
  x = TRUE,
  cluster = IDAA
)

tar_load(JM_CD4_allDLI_inter)
JM_CD4_allDLI_inter$coefficients$gammas
JM_CD4_allDLI_inter$coefficients$alpha


# Look at multivar ones ---------------------------------------------------


coef(multivar_allDLI_nointer)
coxph(
  Surv(tstart, tstop, status) ~ dli + cd4 + ATG + cd8,
  data = data.frame(df2),
  model = TRUE,
  x = TRUE,
  cluster = IDAA
)
coef(multivar_allDLI_inter_penal)

coxph(
  Surv(tstart, tstop, status) ~ dli * cd4 + ATG +  dli * cd8,
  data = data.frame(df2),
  model = TRUE,
  x = TRUE,
  cluster = IDAA
)

# For tomorrow
tar_load(multivar_DLIgvhd_inter)
multivar_DLIgvhd_inter
multivar_allDLI_nointer

# Make plot here?
