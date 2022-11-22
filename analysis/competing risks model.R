library(targets)
library(dplyr)
library(survival)
library(prodlim)
source(here::here("packages.R"))

tar_load(NMA_preDLI_datasets)
dat_wide <- NMA_preDLI_datasets$wide
colrs <- Manu::get_pal("Hoiho")

table(dat_wide$endpoint_specify7)

fit_compEvents_ITT <- prodlim(
  Hist(endpoint7, endpoint7_s, cens.code = "cens") ~ hirisk,
  data = dat_wide
)

summary(fit_compEvents_ITT, times = c(0, 3, 6))

causes <- list(
  "gvhd" = "Clinically significant GvHD",
  "relapse" = "Relapse",
  "other_nrf" = "Other failure"
)

png(
  here("analysis/figures/compEvents_ITT_all.png"), width = 18, height = 9,
  units = "cm", pointsize = 11, bg = "white", res = 300
)
par(mar = c(6, 4, 2, 0) + 0.1, lend = 1, cex = 0.5, ljoin = 1, mfrow = c(1, 3),
    family = "Roboto Condensed")
for (i in seq_along(causes)) {
  plot(
    fit_compEvents_ITT,
    cause = names(causes)[i],
    atrisk.at = c(0:6),
    xlim = c(0, 6),
    col = c(colrs[[6]], colrs[[5]]),#c("blue", "red"),
    marktime = TRUE,
    atrisk.line = c(3.5, 4.5),
    atrisk.cex = 0.5,
    axes = TRUE,
    percent = FALSE,
    ylab = "Cumulative incidence",
    xlab = "Months since alloSCT\n ",
    axis2.las = 1,
    confint = TRUE,
    lwd = 2,
    automar = FALSE,
    ylim = c(0, 1.01),
    axis2.at = seq(0, 1, 0.1),
    background.border = "transparent",
    background = FALSE, # background.horizontal=(0:4)*0.25,
    plot.main = causes[[i]],
    axis1.pos = 0,
    axis2.pos = 0,
    legend = FALSE,
    atrisk.title = "No. at risk  ",
    atrisk.labels = c("Non-high risk  ", "High risk  ")
  )
}
legend(
  "topright",
  bty = "n",
  legend = c("Non-high risk", "High risk"),
  col = c(colrs[[6]], colrs[[5]]),
  title = "",
  lty = 1,
  lwd = 2,
  inset = c(0.05, 0)
)
dev.off()
