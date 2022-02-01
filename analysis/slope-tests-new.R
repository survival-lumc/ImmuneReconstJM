tar_load(
  c(
    preDLI_CD3_long,
    preDLI_cox,
    NMA_preDLI_datasets
  )
)

preDLI_CD3_long$data <- as.data.frame(preDLI_CD3_long$data)
dat_long <- NMA_preDLI_datasets$long
dat_wide <- NMA_preDLI_datasets$wide

jm_slope <- jointModel(
  lmeObject = preDLI_CD3_long,
  survObject = preDLI_cox,
  timeVar = "intSCT2_5",
  CompRisk = TRUE,
  method = "spline-PH-aGH",
  parameterization = "both",
  iter.EM = 1000,
  interFact = list(
    "value" = ~ strata(trans) - 1,
    "slope" = ~ strata(trans) - 1,
    "data" = preDLI_cox$model
  ),
  derivForm = list(
    fixed = ~ 0 + dns(intSCT2_5, 3) + dns(intSCT2_5, 3):as.numeric(hirisk == "yes") +
      dns(intSCT2_5, 3):as.numeric(ATG == "ALT+ATG") +
      dns(intSCT2_5, 3):as.numeric(hirisk == "yes"):as.numeric(ATG == "ALT+ATG"),
    random = ~ 0 + dns(intSCT2_5, 3),
    indFixed = c(2:4, 8:13, 15:17),
    indRandom = c(2:4)
  )
)

# Try the tangentsss
summary(jm_slope)
object <- jm_slope
dat_long[, "currval" := fitted(
  jm_slope,
  process = "Longitudinal",
  type = "Subject"
)]
dat_long[, "slope" := ff]

get_int_tangent <- function(x, y, slope) {
  m <- slope * x
  int <- -x * slope + y
  return(int)
}

dat_long[, int_tang := get_int_tangent(
  x = intSCT2_5, y = currval, slope = slope
)]

delta_tan <- 0.5
dat_long[, ':=' (
  start_x = intSCT2_5 - delta_tan,
  start_y = int_tang + slope * (intSCT2_5 - delta_tan),
  end_x = intSCT2_5 + delta_tan,
  end_y = int_tang + slope * (intSCT2_5 + delta_tan)
)]

set.seed(08908)
dat_long[IDAA %in% sample(IDAA, replace = FALSE, size = 56)] |>
  ggplot(aes(intSCT2_5, CD3_abs_log, group = IDAA)) +
  geom_point(size = 0.75) +
  geom_line(aes(y = currval), size = 0.5) +
  # geom_segment(
  #   aes(x = start_x, y = start_y, xend = end_x, yend = end_y),
  #   size = 1,
  #   arrow = arrow(length = unit(0.03, "npc"))
  # ) +
  facet_wrap(~ IDAA)


dat_long[IDAA %in% sample(IDAA, replace = FALSE, size = 9)] |>
  ggplot(aes(intSCT2_5, CD3_abs_log, group = IDAA)) +
  geom_point(size = 0.75) +
  geom_line(aes(y = currval), size = 0.5) +
  geom_segment(
    aes(x = start_x, y = start_y, xend = end_x, yend = end_y),
    size = 1,
    arrow = arrow(length = unit(0.03, "npc"))
  ) +
  facet_wrap(~ IDAA)

## ANSWER: include dataset by interfact; and make function for slope...
