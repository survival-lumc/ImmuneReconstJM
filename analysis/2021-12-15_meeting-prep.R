library(targets)
library(data.table)
library(ggplot2)
library(magrittr)
library(kableExtra)
library(JM)

tar_load(
  c(NMA_preDLI_CD3_jointModel_both, NMA_preDLI_datasets)
)

dat_long <- NMA_preDLI_datasets$long
dat_wide <- NMA_preDLI_datasets$wide

summary(NMA_preDLI_CD3_jointModel_both)


newdat_jm <- expand.grid(
  "ATG" = levels(dat_long$ATG),
  "CMV_PD" = levels(dat_long$CMV_PD),
  "hirisk" = levels(dat_long$hirisk),
  "intSCT2_5" = seq(0, 6, by = 0.05)
)

predict(
  NMA_preDLI_CD3_jointModel_both,
  newdata = newdat_jm,
  type = "Subject",
  idVar = "IDAA",
  returnData = TRUE,
  interval = "confidence"
) |>
  ggplot(
    aes(
      x = intSCT2_5,
      y = pred,
      group = interaction(ATG, hirisk),
      col = interaction(ATG, hirisk)
    )
  ) +
  geom_ribbon(
    aes(ymin = low, ymax = upp),
    fill = "gray",
    alpha = 0.3,
    col = NA
  ) +
  geom_line(size = 1.5, alpha = 0.75) +
  #geom_text(aes(label = label), hjust = 0, na.rm = TRUE, fontface = "bold") +
  facet_wrap(facets = ~ CMV_PD) +
  labs(x = "Time since alloHCT (months)", y = "CD3 cell counts") +
  #scale_y_continuous(
  #  breaks = log(c(5, 25, 100, 500, 1500)),
  #  labels = c(5, 25, 100, 500, 1500)
  #) +
  xlim(c(0, 6)) +
  theme_minimal(base_size = 14)

# Event summary table
summ_mod <- summary(NMA_preDLI_CD3_jointModel_both)
surv_tab <- summ_mod$`CoefTable-Event`
surv_tab


# Make table
summary_raw <- data.table::data.table(surv_tab, keep.rownames = TRUE)
data.table::setnames(summary_raw, "rn", "Coefficient")
dt <- summary_raw[grep("^bs", x = Coefficient, invert = TRUE)]

dt[, event_num := regmatches(
  x = Coefficient,
  m = regexpr(pattern = "[+1-9]$", text = Coefficient)
)]

dt[, event := factor(
  x = event_num,
  levels = seq_len(3),
  labels = c("GVHD", "Relapse", "NRF: Other")
)]

dt[, Coefficient := gsub(x = Coefficient, pattern = "\\:.*", replacement = "")]
dt[, Coefficient := gsub(x = Coefficient, pattern = "\\.[+1-9]$", replacement = "")]
dt[, Coefficient := factor(
  x = Coefficient,
  levels = c("ATG", "hirisk", "Assoct", "Assoct.s"),
  labels = c("ATG", "High risk (ITT)", "Current value", "Slope")
)]
data.table::setorder(dt, event, Coefficient)
dt_edit <- dt[, c("event", "Coefficient", "Value", "Std.Err", "p-value")]
data.table::setnames(dt_edit, new = c("Event", "Coefficient", "log(HR)", "SE", "p-value"))
dt_edit[, "HR [95% CI]" := paste0(
  as.character(round(exp(`log(HR)`), 3)), " [",
  as.character(round(exp(`log(HR)` - qnorm(0.975) * SE), 3)), ";",
  as.character(round(exp(`log(HR)` + qnorm(0.975) * SE), 3)), "]"
)]
dt_edit[, Event := NULL]
dt_edit[, c("log(HR)", "SE", "p-value") := lapply(.SD, function(x) round(x, 3)),
        .SDcols = c("log(HR)", "SE", "p-value")]
dt_edit

dt_edit %>%
  kableExtra::kbl(format = "html") %>%
  kable_paper("striped", full_width = F) %>%
  pack_rows("GVHD", 1, 3) %>%
  pack_rows("Relapse", 4, 7) %>%
  pack_rows("NRF: Other", 8, 10)


# Plot current slopes?
nrow(NMA_preDLI_CD3_jointModel_both$x$X)
plot(NMA_preDLI_CD3_jointModel_both)


#plot() basis derivative
plot(dat_long$intSCT2_5,
     NMA_preDLI_CD3_jointModel_both$x$Xtime[, 2])


# Model just with cause-specific Relapse (or GVHD on its own)
