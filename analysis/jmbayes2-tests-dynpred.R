# ONLY IF TIME (for fun):
# Try dynamic prediction with new jmbayes2 (also to double check results)
# - need to use cr_setup() (interaction parametrisation)
# - Try bivariate CD4 and CD8, or just CD3

tar_load(c(dat_merged, NMA_preDLI_datasets))

# This was syntax with mstate parametrisation (just for guidance)
jm_mod <- jm(
  Surv_object = preDLI_cox,
  Mixed_objects = list(preDLI_CD3_long_reffcorr),
  time_var = "intSCT2_5",
  #"CD3_abs_log" = ~ value(CD3_abs_log) + slope(CD3_abs_log) +
  #  (value(CD3_abs_log) + slope(CD3_abs_log)):(strata(trans) - 1)
  functional_forms = list(
    "CD3_abs_log" = ~ value(CD3_abs_log) + value(CD3_abs_log):(strata(trans) - 1)
  ),
  id_var = "IDAA",
  n_iter = 13000L,
  n_burnin = 3000L
  #priors = list("penalty_alphas" = "horseshoe") #"horseshoe")
  # Numver of alphas
  #priors = list(
  #  Tau_alphas = lapply(seq_len(3), function(alpha) matrix(data = 0.75))
  #),
  #base_hazard_segments = 5
)

ggtraceplot(jm_mod, "alphas", gridrows = 1, gridcols = 3, grid = TRUE)
ggtraceplot(jm_mod, "betas", gridrows = 6, gridcols = 3, grid = TRUE)

gelman_diag(jm_slopes)



# Try dynamic pred --------------------------------------------------------

dat_long <- NMA_preDLI_datasets$long
dat_wide <- NMA_preDLI_datasets$wide

# crisk_setup(...)
# ...
