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


tar_load(
  c(
    NMA_preDLI_datasets,
    preDLI_tdc_dataset,
    preDLI_JM_value_corr_NK
  )
)
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

newdat_preDLI <- expand.grid(
  "ATG" = levels(dat_long_NK$ATG),
  "CMV_PD" = levels(dat_long_NK$CMV_PD),
  "hirisk" = levels(dat_long_NK$hirisk),
  "intSCT2_7" = seq(0, 6, by = 0.05)
)

NK_JM <- preDLI_JM_value_corr_NK
summary(NK_JM)
round(exp(confint(NK_JM, parm = "Event")), 3)

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
    linetype = "Risk group",
    col = "Donor type"
  ) +
  scale_color_manual(
    labels = c("RD", "UD(+ATG)"),
    values = c(colrs[[6]], colrs[[1]])
  ) +
  scale_linetype_manual(
    labels = c("Non-high risk", "High risk"),
    values = c("solid", "dotdash")
  ) +
  facet_wrap(
    ~ CMV_PD,
    labeller = as_labeller(
      c("-/-" = "-/-", "other P/D" = "Other")
    )
  )

#confin

ggsave(
  here("analysis/figures/preDLI_NK_marginals-test.png"),
  dpi = 300,
  width = 9,
  height = 5
)

ggsave(
  here("analysis/figures/figure07.pdf"), # also to eps/tiff?
  dpi = 300,
  #scale = 1.8,
  #units = "mm",
  width = 9,
  height = 5,
  device = cairo_pdf
)

ggsave(
  here("analysis/figures/figure07.eps"),
  dpi = 300,
  #scale = 1.8,
  #units = "mm",
  width = 9,
  height = 5,
  device = cairo_ps
)



# Try new figure ----------------------------------------------------------

predict(
  NK_JM,
  newdata = newdat_preDLI,
  type = "Marginal",
  idVar = "IDAA",
  returnData = TRUE,
  interval = "confidence"
) |>
  ggplot(aes(x = intSCT2_7, y = pred, group = interaction(CMV_PD, hirisk))) +
  geom_ribbon(
    aes(ymin = low, ymax = upp),
    fill = confint_col,
    alpha = confint_alpha,
    col = NA
  ) +
  log_axis_scales +
  geom_line(aes(linetype = hirisk, col = CMV_PD), size = 1.5) +
  labs(
    x = "Time since alloSCT (months)",
    y = expression(paste("NK cell count (x10"^"6","/l)")),
    linetype = "ITT",
    col = "CMV"
  ) +
  scale_color_manual(
    labels = c("-/-", "other P/D"),
    values = c(colrs[[6]], colrs[[1]])
  ) +
  scale_linetype_manual(
    labels = c("No", "Yes"),
    values = c("solid", "dotdash")
  ) +
  facet_wrap(~ ATG)

# or
predict(
  NK_JM,
  newdata = newdat_preDLI,
  type = "Marginal",
  idVar = "IDAA",
  returnData = TRUE,
  interval = "confidence"
) |>
  ggplot(aes(x = intSCT2_7, y = pred, group = CMV_PD)) +
  geom_ribbon(
    aes(ymin = low, ymax = upp),
    fill = confint_col,
    alpha = confint_alpha,
    col = NA
  ) +
  log_axis_scales +
  geom_line(aes(linetype = CMV_PD, col = CMV_PD), size = 1.5) +
  labs(
    x = "Time since alloSCT (months)",
    y = expression(paste("NK cell count (x10"^"6","/l)")),
    linetype = "CMV",
    col = "CMV"
  ) +
  scale_color_manual(
    labels = c("-/-", "other P/D"),
    values = c(colrs[[6]], colrs[[1]])
  ) +
  scale_linetype_manual(
    labels = c("-/-", "other P/D"),
    values = c("solid", "dotdash")
  ) +
  facet_grid(
    hirisk ~ ATG,
    labeller = as_labeller(
      c(
        "no" = "ITT: no",
        "yes" = "ITT: yes",
        "RD" = "RD",
        "UD(+ATG)" = "UD(+ATG)"
      )
    )
  )



# TDC bit -----------------------------------------------------------------


# Compare results with preDLI_Cox - it's the same
mod_standard <- coxph(
  form = Surv(tstart, tstop, status) ~
    ATG.1 + hirisk.1 + # GVHD
    ATG.2 + hirisk.2 + # Relapse
    ATG.3 + # NRF other
    strata(trans),
  data = preDLI_tdc_dataset
)
coef(tar_read(preDLI_cox))
coef(mod_standard)

# Mod full tings - num events = 20, 22, 33
endpoints <- c("gvhd", "relapse", "other_nrf")
cell_vars <- as.character(
  outer(
    X = paste0("log_", c("CD4", "NK")), #"CD8"
    Y = seq_along(endpoints),
    FUN = paste,
    sep = "."
  )
)
cell_vars_rhs <- reformulate(
  #termlabels = c(cell_vars, "strata(trans)"),
  #termlabels = c(paste0("ATG.", 1:3), cell_vars, "strata(trans)"),
  termlabels = c(".", cell_vars),
  response = "."
)
cell_vars_rhs

mod_full <- update(mod_standard, cell_vars_rhs)
forestmodel::forest_model(
  mod_full,
  exponentiate = TRUE,
  format_options = forestmodel::forest_model_format_options(point_size = 2, text_size = 4)
)

