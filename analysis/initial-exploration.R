# Add parts of this to data-raw and process to data
devtools::load_all()

# Load data
lymphocytes <- data.table::data.table(readRDS("data-raw/2021-02-05_v2/lymphocytes.rds"))
variables <- data.table::data.table(readRDS("data-raw/2021-02-05_v2/variables.rds"))


# Add to description file later
library(ggplot2)
library(magrittr)
library(lme4)
library(nlme)
library(survival)
library(JM)
library(data.table)
library(mstate)


theme_set(theme_bw(base_size = 14))

p <- data.table::melt.data.table(
  data = lymphocytes,
  id.vars = c("IDAA", "intSCT"),
  value.name = "cell_counts",
  measure.vars = c("NK_abs", paste0("CD", c(3, 4, 8 , 19), "_abs")),
  variable.name = "lymph_subset"
) %>%
  ggplot(aes(intSCT, cell_counts, col = factor(IDAA))) +
  geom_line(alpha = 0.5) +
  facet_wrap(. ~ lymph_subset, scales = "free") +
  labs(x = "Time since alloHCT (days)", y = "Cell counts") +
  theme(legend.position = "none")

p

# Same but on log scale (some zero values causing weird things)
p + scale_y_continuous(trans = "log")

# Try plotly
#plotly::ggplotly(p)

# Test merge
full_dat <- data.table::merge.data.table(
  lymphocytes, variables,
  by = c("IDAA", "endpoint", "endpoint2", "endpoint_s", "endpoint2_s")
)

factors <- which(sapply(full_dat, is.factor))
full_dat[, (factors) := lapply(.SD, droplevels), .SDcols = factors]

#
saveRDS(full_dat, file = "data/merged-data.rds")


# Fit an initial JM -------------------------------------------------------


# Check prevalence of zero measurments
full_dat[, lapply(.SD, function(col) {
  sum(is.na(col))
}), .SDcols = c("CD3_abs", "CD4_abs", "CD8_abs", "CD19_abs", "NK_abs")]


# Subset the necessary for simplified JM - use CD4, age, patsex
JM_dat <- full_dat[, c(
  "IDAA",
  "CD4_abs", # CD4 counts
  "intSCT2_2", # Time of measurement (-1 day if at endpoint)
  "EFS3_s", # 1 if REL/ non-relapse fail, zero other wise
  "endpoint2", # Time to endpoint
  "PATSEX",
  "age_SCT"
)]

JM_dat

# Also remove keep non-NA measurements
JM_dat <- JM_dat[!is.na(CD4_abs), ]

# If any CD4 counts = zero, set to 0.5
JM_dat[CD4_abs == 0, CD4_abs := 0.5]

# Make second var for CD4 counts on log scale
JM_dat[, CD4_log := log(CD4_abs)]

# Make sure no measurements taken AFTER endpoint
JM_dat <- JM_dat[intSCT2_2 < endpoint2, ]

# Check number of measurements per patient
JM_dat[, .(n_measures = .N), by = "IDAA"] %>% View() # Two pats have only 1 measurement.. remove
JM_dat <- JM_dat[, .SD[.N > 1], by = "IDAA"]

# Look at remaining measurements
JM_dat[, .(n_measures = .N), by = "IDAA"][["n_measures"]] %>%
  hist(breaks = 15, main = "Number of measurements per pat (pre-endpoint)")

# One patient with very long follow-up (time is in months)
hist(JM_dat$endpoint2, breaks = 20)

# Check trajectories on log and normal scale
p <- JM_dat %>%
  melt.data.table(
    id.vars = c("IDAA", "intSCT2_2"),
    measure.vars = c("CD4_abs", "CD4_log"),
    variable.name = "scale",
    value.name = "CD4_counts"
  ) %>%
  ggplot(aes(intSCT2_2 , CD4_counts, col = factor(IDAA))) +
  geom_point(alpha = 0.5) +
  geom_line(aes(group = factor(IDAA)), alpha = 0.25) +
  facet_wrap(. ~ scale, scales = "free_y") +
  labs(x = "Time since alloHCT (months)", y = "Cell counts") +
  theme(legend.position = "none") +
  coord_cartesian(xlim = c(0, 6))

ggsave(filename = "cd4_counts.png", plot = p, dpi = 300, width = 10, height = 6)

# JM dat is already prepped for a mixed model
lmeFit <- lme(
  CD4_log ~ intSCT2_2 + PATSEX,
  random = ~ intSCT2_2 | IDAA,
  data = JM_dat
)
summary(lmeFit)

# Lets plot the (linear) mixed model fit already
p_lmm <- cbind(JM_dat, "pred" = predict(lmeFit, JM_dat)) %>%
  .[, IDAA := as.factor(IDAA)] %>%
  ggplot(aes(intSCT2_2, CD4_log)) +
  geom_point(aes(col = IDAA), alpha = 0.25) +
  geom_line(aes(y = pred, group = IDAA, col = IDAA), alpha = 0.5) +
  labs(x = "Time from alloHCT", y = "CD4 counts (log)") +
  facet_wrap(~ PATSEX) +
  theme(legend.position = "none") +
  coord_cartesian(xlim = c(0, 6)) # Limit to first 2 years

p_lmm



# Lets put out data back into wide from for the Cox model
JM_dat_wide <- dcast(
  JM_dat,
  formula = IDAA + PATSEX + age_SCT + EFS3_s + endpoint2 ~ .,
  fun = length,
)

# Fit cox model with just sex
coxFit <- coxph(Surv(endpoint2, EFS3_s) ~ PATSEX, data = JM_dat_wide, x = TRUE)
summary(coxFit)


# Now we have everything to fit joint model (linear assoc)
jointFit <- jointModel(
  lmeObject = lmeFit,
  survObject = coxFit,
  timeVar = "intSCT2_2" #
  #method = "weibull-PH-GH" # Assume weibull basleine hazard, more efficient
)

summary(jointFit)
