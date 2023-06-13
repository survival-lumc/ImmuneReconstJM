library("targets")
library("tarchetypes")
library("future")
library("future.callr")
library("here")

source(here::here("packages.R"))
source(here::here("data-raw/prepare-raw-data.R"))
functions <- list.files(here("R"), full.names = TRUE)
invisible(lapply(functions, source))
rm("functions")

# Global plots theme + settings
colrs <- Manu::get_pal("Hoiho")
confint_alpha <- 0.5
confint_col <- "lightgray"
global_font <- "Roboto Condensed"
log_axis_scales <- scale_y_continuous(
  breaks = log(c(0.1, 1, 5, 25, 100, 500, 1500)),
  labels = c(0.1, 1, 5, 25, 100, 500, 1500)
)

base_size_png <- 14

theme_set(
  theme_light(base_size = base_size_png, base_family = global_font) +
    theme(
      strip.background = element_rect(fill = colrs[[2]], colour = "transparent"),
      strip.text = element_text(colour = 'white')
    )
)



# Start -------------------------------------------------------------------


tar_load(NMA_preDLI_datasets)
dat_long_NK <- NMA_preDLI_datasets$long

dat_long_NK[, "endpoint_lab" := factor(
  fcase(
    endpoint7_s == "cens", "Censored",
    endpoint7_s == "gvhd", "GvHD",
    endpoint7_s == "relapse", "Relapse",
    endpoint7_s == "other_nrf", "Other failure"
  ),
  levels = c("Censored", "GvHD", "Relapse", "Other failure")
)]


ggplot(dat_long_NK, aes(intSCT2_7, NK_abs_log, group = IDAA)) +
  geom_line(show.legend = FALSE, linewidth = 1, alpha = 0.5, col = colrs[[6]]) +
  labs(x = "Time since alloSCT (months)", y = expression(paste("NK cell count (x10"^"6","/l)"))) +
  log_axis_scales +
  facet_wrap(. ~ endpoint_lab) +
  coord_cartesian(ylim = log(c(0.1, 3000)))

ggsave(
  here("analysis/figures/preDLI_NK_perEndpoint.png"),
  dpi = 300,
  width = 8,
  height = 8
)

mod_long <- run_preDLI_longitudinal(
  cell_line = "NK_abs_log",
  form_fixed = "ns(intSCT2_7, 3) * hirisk * ATG + CMV_PD",
  form_random = ~ 0 + ns(intSCT2_7, 3) | IDAA,
  dat = dat_long_NK
)

newdat_preDLI <- expand.grid(
  "ATG" = levels(dat_long_NK$ATG),
  "CMV_PD" = levels(dat_long_NK$CMV_PD),
  "hirisk" = levels(dat_long_NK$hirisk),
  "intSCT2_7" = seq(0, 6, by = 0.05)
)

cbind(
  newdat_preDLI,
  preds = predict(mod_long, level = 0L, newdata = newdat_preDLI)
) |>
  ggplot(aes(intSCT2_7, preds, col = ATG)) +
  log_axis_scales +
  geom_line(aes(linetype = hirisk), linewidth = 1.5) +
  facet_wrap(~ CMV_PD)

mod_cox <- run_preDLI_cox(
  form = Surv(time, status) ~
    ATG.1 + hirisk.1 + # GVHD
    ATG.2 + hirisk.2 + # Relapse
    ATG.3 + # NRF other
    strata(trans),
  dat_wide = NMA_preDLI_datasets$wide
)

NK_JM <- jointModel(
  lmeObject = mod_long,
  survObject = mod_cox,
  CompRisk = TRUE,
  method = "spline-PH-aGH",
  timeVar = "intSCT2_7",
  parameterization = "value",
  interFact = list("value" = ~ strata(trans) - 1),
  iter.qN = 1000L,
  lng.in.kn = tar_read(preDLI_basehaz_knots),
  numeriDeriv = "cd",
  eps.Hes = 1e-04,
  verbose = TRUE
)

summary(NK_JM)

predict(
  NK_JM,
  newdata = newdat_preDLI,
  type = "Marginal",
  idVar = "IDAA",
  returnData = TRUE,
  interval = "confidence"
) |>
  ggplot(aes(x = intSCT2_7, y = pred, group = interaction(ATG, hirisk))) +
  geom_ribbon(
    aes(ymin = low, ymax = upp),
    fill = confint_col,
    alpha = confint_alpha,
    col = NA
  ) +
  log_axis_scales +
  geom_line(aes(linetype = hirisk, col = ATG), size = 1.5) +
  labs(
    x = "Time since alloSCT (months)",
    y = expression(paste("NK cell count (x10"^"6","/l)")),
    linetype = "ITT",
    col = "Donor type"
  ) +
  scale_color_manual(
    labels = c("RD", "UD(+ATG)"),
    values = c(colrs[[6]], colrs[[1]])
  ) +
  scale_linetype_manual(
    labels = c("No", "Yes"),
    values = c("solid", "dotdash")
  ) +
  facet_wrap(~ CMV_PD)

ggsave(
  here("analysis/figures/preDLI_NK_marginals-test.png"),
  dpi = 300,
  width = 9,
  height = 5
)
