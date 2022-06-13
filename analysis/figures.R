# figures for manuscript. Code based on/copied from meeting prep scripts

library(targets)
library(dplyr)
library(data.table)
library(ggplot2)
library(JM)

theme_set(theme_bw(base_size = 14))

tar_load(
  c(
    NMA_preDLI_datasets,
    NMA_postDLI_datasets,
    #preDLI_cox,
    #preDLI_long_corr_CD4,

    # preDLI: value
    preDLI_JM_value_corr_CD3,
    preDLI_JM_value_corr_CD4,
    preDLI_JM_value_corr_CD8,

    # preDLI: value + slope
    preDLI_JM_both_corr_CD3,
    preDLI_JM_both_corr_CD4,
    preDLI_JM_both_corr_CD8,

    # postDLI: value
    postDLI_JM_corr_CD3,
    postDLI_JM_corr_CD4,
    postDLI_JM_corr_CD8
  )
)

source("R/modelling-helpers.R")
source("R/plotting-helpers.R")

#### preDLI ####
dat_wide_pre <- NMA_preDLI_datasets$wide
dat_long_pre <- NMA_preDLI_datasets$long

dat_long_pre[, "endpoint_lab" := fcase(
  endpoint7_s == "cens", "Censored",
  endpoint7_s == "gvhd", "GvHD", # instead of GvHD
  endpoint7_s == "relapse", "Relapse",
  endpoint7_s == "other_nrf", "Other failure" # instead of Other NRF
)]
dat_long_pre[, endpoint_lab := factor(
  endpoint_lab, levels = c("Censored", "GvHD", "Relapse", "Other failure")
)]

# Added fitted vals and slope
dat_long_pre[, ':=' (
  curr_val_value = fitted(preDLI_JM_value_corr_CD3, type = "Subject"),
  curr_val_both = fitted(preDLI_JM_both_corr_CD3, type = "Subject"),
  slope = fitted_slopes_long(preDLI_JM_both_corr_CD3)
)]

# number of measurements
ggplot(dat_wide_pre, aes(n_measurements)) +
  geom_histogram(bins = 13, col = "black", fill = "lightblue") +
  labs(x = "# Repeated measurements per patient") +
  theme_minimal(base_size = 14)

# For plotting: get tangent intercept, and then get coordinates
dat_long_pre[, int_tang := get_int_tangent(x = intSCT2_7, y = curr_val_value, slope = slope)]
delta_tan <- 0.15 # Width of tangent arrow
dat_long_pre[, ':=' (
  start_x = intSCT2_7 - delta_tan,
  start_y = int_tang + slope * (intSCT2_7 - delta_tan),
  end_x = intSCT2_7 + delta_tan,
  end_y = int_tang + slope * (intSCT2_7 + delta_tan)
)]

# Make a data with last measurement
dat_long_pre_last <- dat_long_pre[, .SD[.N], by = "IDAA"]

# Sample 16 patients
set.seed(649846561)
IDAA_subs_pre <- sample(levels(dat_long_pre$IDAA), replace = FALSE, size = 16)

# pre-DLI: raw + indiv fits -----------------------------------------------

# First raw-plot
p_raw_indiv_pre <- ggplot(dat_long_pre[IDAA %in% IDAA_subs_pre], aes(intSCT2_7, CD3_abs_log)) +
  #ggplot(dat_long_pre, aes(intSCT2_7, CD3_abs_log)) +
  geom_point(
    size = 3.5, pch = 21,
    alpha = 0.8,
    col = "#359fda",
    fill = colorspace::lighten("#359fda", amount = 0.3)
  ) +
  geom_vline(aes(xintercept = endpoint7), linetype = "dashed") +
  facet_wrap(~ as.numeric(factor(IDAA))) +  # instead of IDAA to get 1:16
  labs(x = "Time since alloSCT (months)", y = expression(paste("CD3 (x10"^"6","/l)"))) +
  scale_y_continuous(
    breaks = log(c(5, 25, 100, 500, 1500)),
    labels = c(5, 25, 100, 500, 1500)
  ) +

  # Label endpoint
  geom_label(
    data = dat_long_pre_last[IDAA %in% IDAA_subs_pre],
    aes(x = endpoint7 + 0.25, y = log(2.5), label = endpoint_lab),
    #size = 5,
    hjust = 0,
    lineheight = .8,
    inherit.aes = FALSE,
    label.size = NA
  ) +
  geom_curve(
    data = dat_long_pre_last[IDAA %in% IDAA_subs_pre],
    #data = data.frame(x = 2, y = 29, xend = 2.5, yend = 20),
    mapping = aes(x = endpoint7 + 0.25, y = log(2.5), xend = endpoint7 + 0.025, yend = log(5)),
    colour = "black",
    size = 0.75,
    curvature = 0.2,
    arrow = arrow(length = unit(0.05, "npc"), type = "open"),
    inherit.aes = FALSE
  ) +
  coord_cartesian(xlim = c(0, 7.5))

p_raw_indiv_pre
ggsave("analysis/figures/preDLI_raw_indiv.png",dpi=300,width=12,height=8)

# Add individual fits
p_raw_indiv_pre +
  geom_line(aes(y = curr_val_value), size = 1, col = "darkblue")
ggsave("analysis/figures/preDLI_lines_indiv.png",dpi=300,width=12,height=8) # this is Figure 1


# all pre-DLI trajectories per endpoint -----------------------------------
ggplot(dat_long_pre, aes(intSCT2_7, CD3_abs_log, color=IDAA)) + # CD3_abs_log instead of curr_val to get raw values
  geom_line(show.legend = FALSE)+
  labs(x = "Time since alloSCT (months)", y = expression(paste("CD3 (x10"^"6","/l)"))) +
  scale_y_continuous(
    breaks = log(c(5, 25, 100, 500, 1500)),
    labels = c(5, 25, 100, 500, 1500)
  ) +
  facet_wrap(~endpoint_lab, labeller=as_labeller(c("Censored"="Censored", "GvHD"="GvHD", "Relapse"="Relapse", "Other failure"="Other failure")))
ggsave("analysis/figures/preDLI_perEndpoint.png",dpi=300,width=8,height=8) # this is Supplemental Figure 1


# pre-DLI: marginal fit ---------------------------------------------------

mod_names <- c(
  paste0("preDLI_JM_value_corr_CD", c(3, 4, 8)),
  paste0("preDLI_JM_both_corr_CD", c(3, 4, 8)),
  paste0("postDLI_JM_corr_CD", c(3, 4, 8))
)
names(mod_names) <- mod_names
mods <- lapply(mod_names, tar_read_raw)

table_mod_summary <- function(model) {

  # Event summary extract
  summ_mod <- summary(model)
  surv_tab <- summ_mod$`CoefTable-Event`

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
    labels = c("GvHD", "Relapse", "NRF: Other")
  )]

  dt[, Coefficient := gsub(x = Coefficient, pattern = "\\:.*", replacement = "")]
  dt[, Coefficient := gsub(x = Coefficient, pattern = "\\.[+1-9]$", replacement = "")]
  dt[, Coefficient := factor(
    x = Coefficient,
    levels = c("ATG", "hirisk", "Assoct", "Assoct.s"),
    labels = c("Unrel. donor (+ATG)", "High risk (ITT)", "Current value", "Slope")
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
  # Change order?
  dt_edit <- dt_edit[, c("Coefficient", "HR [95% CI]", "log(HR)", "SE", "p-value")]
  dt_edit[, c("log(HR)", "SE", "p-value") := lapply(.SD, function(x) round(x, 3)),
          .SDcols = c("log(HR)", "SE", "p-value")]

  return(dt_edit)
}

table(dat_wide_pre$endpoint7_s)
table(dat_wide_pre$endpoint7_s, dat_wide_pre$ATG)

# Marginals
mods_preDLI_value <- list(
  "CD3" = mods$preDLI_JM_value_corr_CD3,
  "CD4" = mods$preDLI_JM_value_corr_CD4,
  "CD8" = mods$preDLI_JM_value_corr_CD8
)

# For marginal fit
newdat_preDLI <- expand.grid(
  "ATG" = levels(dat_long_pre$ATG),
  "CMV_PD" = levels(dat_long_pre$CMV_PD),
  "hirisk" = levels(dat_long_pre$hirisk),
  "intSCT2_7" = seq(0, 6, by = 0.05)
)

marg_preds_preDLI <- lapply(mods_preDLI_value, function(mod) {
  predict(
    mod,
    newdata = newdat_preDLI,
    type = "Marginal",
    idVar = "IDAA",
    returnData = TRUE,
    interval = "confidence"
  )
})

# Good plot with one
rbindlist(marg_preds_preDLI, idcol = "cell_line") |>
  ggplot(aes(x = intSCT2_7, y = pred, group = interaction(ATG, hirisk))) +
  geom_ribbon(aes(ymin = low, ymax = upp), fill = "gray", alpha = 0.5, col = NA) +
  geom_line(aes(linetype = hirisk, col = ATG), size = 1.5, alpha = 0.75) +
  # Add repel labels?
  facet_grid(facets = cell_line ~ CMV_PD,
             labeller = as_labeller(c("-/-"="CMV: -/-", "other P/D"="CMV: other",
                                      "CD3"="CD3", "CD4"="CD4", "CD8"="CD8"
                                      #"CD3"="total T-cells", "CD4"="CD4+ T-cells", "CD8"="CD8+ T-cells"
                                      ))) +
  labs(
    x = "Time since alloSCT (months)",
    y = expression(paste("cell count (x10"^"6","/l)")),
    linetype = "ITT",
    col = "Donor type"
  ) +
  scale_y_continuous(
    breaks = log(c(5, 25, 100, 500, 1500)),
    labels = c(5, 25, 100, 500, 1500)
  ) +
  coord_cartesian(xlim = c(0, 6), ylim = c(log(0.1), log(1500))) + # add a common ylim for all
  scale_color_manual(
    labels = c("RD", "UD+ATG"),
    values = c("brown", "darkblue")
  ) +
  scale_linetype_manual(
    labels = c("No", "Yes (high risk)"),
    values = c("solid", "dotdash")
  ) +
  theme(legend.position = "bottom")
ggsave("analysis/figures/preDLI_trajectories.png",dpi=300,width=6,height=8) # this is Figure 2


# Summary of post DLI models ----------------------------------------------

dat_wide_postDLI <- NMA_postDLI_datasets$wide
dat_long_postDLI <- NMA_postDLI_datasets$long

table(dat_wide_postDLI$sec_endpoint2_s)
table(dat_wide_postDLI$sec_endpoint2_s, dat_wide_postDLI$ATG) # Silly to estimate effect of ATG?

newdat_postDLI <- expand.grid(
  "ATG" = levels(dat_long_postDLI$ATG),
  "CMV_PD" = levels(dat_long_postDLI$CMV_PD),
  "intDLI1" = seq(0, 3, by = 0.01)
)

mods_postDLI_value <- list(
  "CD3" = mods$postDLI_JM_corr_CD3,
  "CD4" = mods$postDLI_JM_corr_CD4,
  "CD8" = mods$postDLI_JM_corr_CD8
)

marg_preds_postDLI <- lapply(mods_postDLI_value, function(mod) {
  predict(
    mod,
    newdata = newdat_postDLI,
    type = "Marginal",
    idVar = "IDAA",
    returnData = TRUE,
    interval = "confidence"
  )
})


# Good plot with one
rbindlist(marg_preds_postDLI, idcol = "cell_line") |>
  ggplot(aes(x = intDLI1, y = pred, group = ATG)) + # was intDLI1 + 3
  geom_ribbon(aes(ymin = low, ymax = upp), fill = "gray", alpha = 0.5, col = NA) +
  geom_line(aes(col = ATG), size = 1.5, alpha = 0.75) +
  # Add repel labels?
  facet_grid(facets = cell_line ~ CMV_PD,
             labeller = as_labeller(c("-/-"="CMV: -/-", "other P/D"="CMV: other",
                                      "CD3"="CD3", "CD4"="CD4", "CD8"="CD8"
                                      #"CD3"="total T-cells", "CD4"="CD4+ T-cells", "CD8"="CD8+ T-cells"
                                      ))) +
  labs(
    x = "Time since DLI (months)", # was alloSCT
    y = expression(paste("cell count (x10"^"6","/l)")),
    #linetype = "ITT",
    col = "Donor type"
  ) +
  scale_y_continuous(
    breaks = log(c(5, 25, 100, 500, 1500)),
    labels = c(5, 25, 100, 500, 1500)
  ) +
  coord_cartesian(xlim = c(0, 3), ylim = c(log(0.1), log(1500))) + # xlim was c(0,6)
  scale_color_manual(
    labels = c("RD", "UD+ATG"),
    values = c("brown", "darkblue")
  ) +
  theme(legend.position = "bottom") #+
  #geom_vline(xintercept = 3, linetype = "dotted")
ggsave("analysis/figures/postDLI_trajectories.png",dpi=300,width=6,height=8) # this is FIgure 4


# -- Model summaries
ggplot(dat_wide_postDLI, aes(n_measurements)) +
  geom_histogram(bins = 13, col = "black", fill = "lightblue") +
  labs(x = "# Repeated measurements per patient") +
  theme_minimal(base_size = 14)

dat_long_postDLI[, "endpoint_lab" := fcase(
  sec_endpoint2_s == "cens", "Censored",
  sec_endpoint2_s == "gvhd", "GvHD",
  sec_endpoint2_s == "rel_nrf", "Relapse or other failure"
)]

# Added fitted vals and slope
dat_long_postDLI[, ':=' (
  curr_val = fitted(postDLI_JM_corr_CD3, type = "Subject")
)]

# For plotting: get tangent intercept, and then get coordinates
dat_long_postDLI$intSCT2_7_reset=dat_long_postDLI$intSCT2_7-dat_long_postDLI$uDLI_actual
dat_long_postDLI[, int_tang := get_int_tangent(x = intSCT2_7_reset, y = curr_val, slope = 1)]#slope
delta_tan <- 0.3 # Width of tangent arrow
dat_long_postDLI[, ':=' (
  start_x = intSCT2_7_reset - delta_tan,
  start_y = int_tang  * (intSCT2_7_reset - delta_tan), # + slope
  end_x = intSCT2_7_reset + delta_tan,
  end_y = int_tang  * (intSCT2_7_reset + delta_tan) # + slope
)]

# Make a data with last measurement
dat_long_post_last <- dat_long_postDLI[, .SD[.N], by = "IDAA"]

# Sample 16 patients
set.seed(649846561)
IDAA_subs_post <- sample(levels(dat_long_postDLI$IDAA), replace = FALSE, size = 16)

# compare pre and post
any(duplicated(c(IDAA_subs_pre,IDAA_subs_post))) # no -> use numbers 17:32

ggplot(
  dat_long_postDLI[IDAA %in% IDAA_subs_post],
  aes(intSCT2_7_reset, CD3_abs_log)
) +
  geom_point(
    size = 3.5,
    pch = 21,
    alpha = 0.8,
    col = "#359fda",
    fill = colorspace::lighten("#359fda", amount = 0.3)
  ) +
  labs(x = "Time since DLI (months)", y = expression(paste("CD3 (x10"^"6","/l)"))) +
  scale_y_continuous(
    breaks = log(c(5, 25, 100, 500, 1500)),
    labels = c(5, 25, 100, 500, 1500)
  ) +
  facet_wrap(~ (as.numeric(factor(IDAA))+16)) +  # instead of IDAA to get 17:32
  coord_cartesian(ylim = c(0, 8), xlim = c(0, 4.5)) +
  geom_line(aes(y = curr_val), size = 1, col = "darkblue")+
  # endpoint
  geom_vline(aes(xintercept = sec_endpoint2), linetype = "dashed") +
  # Label endpoint
  geom_label(
    data = dat_long_postDLI[IDAA %in% IDAA_subs_post],
    aes(x = sec_endpoint2 + 0.25, y = log(2.5), label = gsub("Event-free", "Censored",endpoint_lab)),
    #size = 5,
    hjust = 0,
    lineheight = .8,
    inherit.aes = FALSE,
    label.size = NA
  ) +
  geom_curve(
    data = dat_long_postDLI[IDAA %in% IDAA_subs_post],
    #data = data.frame(x = 2, y = 29, xend = 2.5, yend = 20),
    mapping = aes(x = sec_endpoint2 + 0.25, y = log(2.5), xend = sec_endpoint2 + 0.025, yend = log(5)),
    colour = "black",
    size = 0.75,
    curvature = 0.2,
    arrow = arrow(length = unit(0.05, "npc"), type = "open"),
    inherit.aes = FALSE
  )
ggsave("analysis/figures/postDLI_lines_indiv.png",dpi=300,width=12,height=8) # this is Supplemental Figure 2

# all pre-DLI trajectories per endpoint -----------------------------------
ggplot(dat_long_postDLI, aes(intSCT2_7_reset, CD3_abs_log, color=IDAA)) + # CD3_abs_log instead of curr_val to get raw values
  geom_line(show.legend = FALSE)+
  geom_point(data=dat_long_postDLI[!dat_long_postDLI$IDAA%in%dat_long_postDLI$IDAA[duplicated(dat_long_postDLI$IDAA)],],
             show.legend = FALSE, size = 0.5)+ # add dots for those with only one measurement
  labs(x = "Time since DLI (months)", y = expression(paste("CD3 (x10"^"6","/l)"))) +
  scale_y_continuous(
    breaks = log(c(5, 25, 100, 500, 1500)),
    labels = c(5, 25, 100, 500, 1500)
  ) +
  facet_wrap(~endpoint_lab)
ggsave("analysis/figures/postDLI_perEndpoint.png",dpi=300,width=8,height=8) # this is Supplemental Figure 3


# forest plots ----------------------------------------------------------
## to be updated

res <- rbind(
  data.table(
    summary(preDLI_JM_both_corr_CD3)$`CoefTable-Event`, "cell_line" = "CD3",
    keep.rownames = TRUE
  ),
  data.table(
    summary(preDLI_JM_both_corr_CD4)$`CoefTable-Event`, "cell_line" = "CD4",
    keep.rownames = TRUE
  ),
  data.table(
    summary(preDLI_JM_both_corr_CD8)$`CoefTable-Event`, "cell_line" = "CD8",
    keep.rownames = TRUE
  )
)

res <- res[grep(x = rn, pattern = "^bs", invert = TRUE)]
setnames(res, old = "rn", new = "param")
res[, param := factor(
  param,
  levels = c(
    "ATG.1", "hirisk.1", "Assoct:strata(trans)trans=1", "Assoct.s:strata(trans)trans=1",
    "ATG.2", "hirisk.2", "Assoct:strata(trans)trans=2", "Assoct.s:strata(trans)trans=2",
    "ATG.3", "Assoct:strata(trans)trans=3", "Assoct.s:strata(trans)trans=3"
  )
)]
res[, param_label := factor(
  param, labels = c(
    "ATG", "hirisk", "Value", "Slope",
    "ATG", "hirisk", "Value", "Slope",
    "ATG", "Value", "Slope"
  )
)]
setorder(res, cell_line, param)
res
res[, ':=' (
  HR = exp(Value),
  low = exp(Value - qnorm(0.975) * Std.Err),
  upp = exp(Value + qnorm(0.975) * Std.Err)
)]


# Try basic one
ggplot(data = res, aes(x = param, y = HR)) +
  geom_linerange(
    aes(
      ymin = low,
      ymax = upp,
      xmin = param,
      xmax = param,
      col = cell_line
    ),
    position = position_dodge(width = 0.75),
    size = 0.5,
    na.rm = TRUE
  ) +
  ggplot2::scale_y_continuous(
    name = "Hazard ratio (95% CI)",
    trans = "log",
    breaks = c(0.05, 0.1, 0.25, 0.5, 1, 2, 4, 8, 15, 25),
  ) +
  ggplot2::scale_x_discrete(limits = rev) +
  geom_point(
    aes(col = cell_line, shape = cell_line),
    position = position_dodge(width = 0.75),
    size = 1.25,
    na.rm = TRUE
  ) +
  guides(col = guide_legend(reverse = TRUE), shape = guide_legend(reverse = TRUE)) +
  coord_flip() +
  geom_hline(yintercept = 1, linetype = "dotted") +
  theme_bw()


####
res2 <- rbind(
  data.table(
    summary(preDLI_JM_value_corr_CD3)$`CoefTable-Event`, "cell_line" = "CD3",
    keep.rownames = TRUE
  ),
  data.table(
    summary(preDLI_JM_value_corr_CD4)$`CoefTable-Event`, "cell_line" = "CD4",
    keep.rownames = TRUE
  ),
  data.table(
    summary(preDLI_JM_value_corr_CD8)$`CoefTable-Event`, "cell_line" = "CD8",
    keep.rownames = TRUE
  )
)

res2 <- res2[grep(x = rn, pattern = "^bs", invert = TRUE)]
setnames(res2, old = "rn", new = "param")
res2[, param := factor(
  param,
  levels = c(
    "ATG.1", "hirisk.1", "Assoct:strata(trans)trans=1",
    "ATG.2", "hirisk.2", "Assoct:strata(trans)trans=2",
    "ATG.3", "Assoct:strata(trans)trans=3"
  )
)]
res2[, param_label := factor(
  param, labels = c(
    "ATG", "hirisk", "Value",
    "ATG", "hirisk", "Value",
    "ATG", "Value"
  )
)]
setorder(res2, cell_line, param)
res2
res2[, ':=' (
  HR = exp(Value),
  low = exp(Value - qnorm(0.975) * Std.Err),
  upp = exp(Value + qnorm(0.975) * Std.Err)
)]


# Try basic one
ggplot(data = res2, aes(x = param, y = HR)) +
  geom_linerange(
    aes(
      ymin = low,
      ymax = upp,
      xmin = param,
      xmax = param,
      col = cell_line
    ),
    position = position_dodge(width = 0.75),
    size = 0.5,
    na.rm = TRUE
  ) +
  ggplot2::scale_y_continuous(
    name = "Hazard ratio (95% CI)",
    trans = "log",
    breaks = c(0.05, 0.1, 0.25, 0.5, 1, 2, 4, 8, 15, 25),
  ) +
  ggplot2::scale_x_discrete(limits = rev) +
  geom_point(
    aes(col = cell_line, shape = cell_line),
    position = position_dodge(width = 0.75),
    size = 1.25,
    na.rm = TRUE
  ) +
  guides(col = guide_legend(reverse = TRUE), shape = guide_legend(reverse = TRUE)) +
  coord_flip() +
  geom_hline(yintercept = 1, linetype = "dotted") +
  theme_bw()
ggsave("analysis/figures/forest_preDLI.png",dpi=300,width=8,height=8) # this is Figure 3

### tryout by Eva
library(dplyr)
res2.2<-res2%>%
  mutate(endpoint=factor(ifelse(grepl("1",param),"GvHD",ifelse(grepl("2",param),"Relapse", "Other failure")),
                         levels=c("GvHD","Relapse","Other failure")),
         param2=factor(ifelse(grepl("ATG",param),"UD+ATG",ifelse(grepl("hirisk",param),"indication for\nearly low-dose DLI", "current count")),
                       levels=c("UD+ATG","indication for\nearly low-dose DLI","current count")),
         subset=factor(cell_line,levels=c("CD8","CD4","CD3")),
         label=paste0(ifelse(HR>=1,"  ", ""), round(HR,digits=2), " (", round(low,digits=2), "-", round(upp,digits=2), ")",
                      ifelse(HR<1, "   \n ", "\n ")))%>%
  arrange(endpoint,param2,subset)%>%
  mutate(label=factor(label,levels=label))
res2.2

ggplot(data = res2.2, aes(x = param2, y = HR)) +
  geom_vline(xintercept=c(1.5,2.5),color="grey92",size=0.5)+
  geom_linerange(
    aes(
      ymin = low,
      ymax = upp,
      xmin = param2,
      xmax = param2,
      col = subset
    ),
    position = position_dodge(width = 0.8),
    size = 1.1,
    na.rm = TRUE,
    show.legend = F
  ) +
  geom_hline(alpha=0,aes(col=subset,yintercept=1),size=1.1)+
  geom_text(aes(x=,param2,label=label,group=subset,hjust=ifelse(HR>=1,0,1)),
            position = position_dodge(width = 0.8),size=3,color="black")+ #grey40
  ggplot2::scale_y_continuous(
    name = "Hazard ratio (95% CI)",
    trans = "log",
    breaks = c(0.05, 0.1, 0.25, 0.5, 1, 2, 4, 8, 15, 25),
  ) +
  ggplot2::scale_x_discrete(limits = rev, expand=expansion(add=0.5)) +
  geom_point(
    aes(col = subset, shape = subset),
    position = position_dodge(width = 0.8),
    size = 2,
    na.rm = TRUE
  ) +
  guides(col = guide_legend(reverse = TRUE,override.aes = list(alpha=1)),
         shape = guide_legend(reverse = TRUE)) +
  coord_flip() +
  geom_hline(yintercept = 1, linetype = "dotted") +
  theme_bw()+
  #facet_wrap(~endpoint,scales = "free_y",ncol=1, strip.position = "right")+
  facet_grid(rows="endpoint",scales="free",space="free")+
  theme(axis.title.y = element_blank(), panel.spacing = unit(0,"cm"),
        axis.text=element_text(color="black",size=11),
        legend.text = element_text(size=11),legend.title=element_text(size=11),
        axis.title = element_text(size=12), strip.text = element_text(size=11),
        panel.grid.major.y=element_blank())
ggsave("analysis/figures/preDLI_forest.png",dpi=300,width=10,height=8)


#### postDLI
res3 <- rbind(
  data.table(
    summary(postDLI_JM_corr_CD3)$`CoefTable-Event`, "cell_line" = "CD3",
    keep.rownames = TRUE
  ),
  data.table(
    summary(postDLI_JM_corr_CD4)$`CoefTable-Event`, "cell_line" = "CD4",
    keep.rownames = TRUE
  ),
  data.table(
    summary(postDLI_JM_corr_CD8)$`CoefTable-Event`, "cell_line" = "CD8",
    keep.rownames = TRUE
  )
)

res3 <- res3[grep(x = rn, pattern = "^bs", invert = TRUE)]
setnames(res3, old = "rn", new = "param")
res3[, param := factor(
  param,
  levels = c(
    "ATG.1", "Assoct:strata(trans)trans=1",
    "Assoct:strata(trans)trans=2"
  )
)]
res3[, param_label := factor(
  param, labels = c(
    "ATG",  "Value",
    "Value"
  )
)]
setorder(res3, cell_line, param)
res3
res3[, ':=' (
  HR = exp(Value),
  low = exp(Value - qnorm(0.975) * Std.Err),
  upp = exp(Value + qnorm(0.975) * Std.Err)
)]


# Try basic one
ggplot(data = res3, aes(x = param, y = HR)) +
  geom_linerange(
    aes(
      ymin = low,
      ymax = upp,
      xmin = param,
      xmax = param,
      col = cell_line
    ),
    position = position_dodge(width = 0.75),
    size = 0.5,
    na.rm = TRUE
  ) +
  ggplot2::scale_y_continuous(
    name = "Hazard ratio (95% CI)",
    trans = "log",
    breaks = c(0.05, 0.1, 0.25, 0.5, 1, 2, 4, 8, 15, 25,50,100,150),
  ) +
  ggplot2::scale_x_discrete(limits = rev) +
  geom_point(
    aes(col = cell_line, shape = cell_line),
    position = position_dodge(width = 0.75),
    size = 1.25,
    na.rm = TRUE
  ) +
  guides(col = guide_legend(reverse = TRUE), shape = guide_legend(reverse = TRUE)) +
  coord_flip() +
  geom_hline(yintercept = 1, linetype = "dotted") +
  theme_bw()
ggsave("analysis/figures/forest_postDLI.png",dpi=300,width=8,height=8) # this is Figure 5


### tryout by Eva
res3.2<-res3%>%
  mutate(endpoint=factor(ifelse(grepl("1",param),"GvHD","Relapse or other failure"),
                         levels=c("GvHD","Relapse or other failure")),
         param2=factor(ifelse(grepl("ATG",param),"UD+ATG", "current count"),
                       levels=c("UD+ATG", "current count")),
         subset=factor(cell_line,levels=c("CD8","CD4","CD3")),
         label=paste0(ifelse(HR>=1,"  ", ""), format(round(HR,digits=2)), " (", round(low,digits=2), "-", round(upp,digits=2), ")",
                      ifelse(HR<1, "   \n \n ", "\n \n ")))%>%
  arrange(endpoint,param2,subset)%>%
  mutate(label=factor(label,levels=label))
res3.2

ggplot(data = res3.2, aes(x = param2, y = HR)) +
  geom_vline(xintercept=c(1.5,2.5),color="grey92",size=0.5)+
  geom_linerange(
    aes(
      ymin = low,
      ymax = upp,
      xmin = param2,
      xmax = param2,
      col = subset
    ),
    position = position_dodge(width = 0.8),
    size = 1.1,
    na.rm = TRUE,
    show.legend = F
  ) +
  geom_hline(alpha=0,aes(col=subset,yintercept=1),size=1.1)+
  geom_text(aes(x=,param2,label=label,group=subset,hjust=ifelse(HR>=1,0,1)),
            position = position_dodge(width = 0.8),size=3,color="black")+ #grey40
  ggplot2::scale_y_continuous(
    name = "Hazard ratio (95% CI)",
    trans = "log",
    breaks = c(0.05, 0.1, 0.25, 0.5, 1, 2, 4, 8, 15, 25,50,100,150),
    labels=c(0.05, "0.10", 0.25, "0.50", "1.00", "2.00", "4.00", "8.00", 15, 25,50,100,150)
  ) +
  ggplot2::scale_x_discrete(limits = rev, expand=expansion(add=0.5)) +
  geom_point(
    aes(col = subset, shape = subset),
    position = position_dodge(width = 0.8),
    size = 2,
    na.rm = TRUE
  ) +
  guides(col = guide_legend(reverse = TRUE,override.aes = list(alpha=1)),
         shape = guide_legend(reverse = TRUE)) +
  coord_flip() +
  geom_hline(yintercept = 1, linetype = "dotted") +
  theme_bw()+
  #facet_wrap(~endpoint,scales = "free_y",ncol=1, strip.position = "right")+
  facet_grid(rows="endpoint",scales="free",space="free")+
  theme(axis.title.y = element_blank(), panel.spacing = unit(0,"cm"),
        axis.text=element_text(color="black",size=11),
        legend.text = element_text(size=11),legend.title=element_text(size=11),
        axis.title = element_text(size=12), strip.text = element_text(size=11),
        panel.grid.major.y=element_blank(),plot.margin=margin(l=0))
ggsave("analysis/figures/postDLI_forest.png",dpi=300,height=6,width=11)
