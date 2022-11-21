# Reproducing manuscript figures ------------------------------------------

source(here::here("packages.R"))

# Global plots theme + settings
colrs <- Manu::get_pal("Hoiho")
confint_alpha <- 0.5
confint_col <- "lightgray"
global_font <- "Roboto Condensed"
log_axis_scales <- scale_y_continuous(
  breaks = log(c(0.1, 1, 5, 25, 100, 500, 1500)),
  labels = c(0.1, 1, 5, 25, 100, 500, 1500)
)

theme_set(
  theme_light(base_size = 14, base_family = global_font) +
    theme(
      strip.background = element_rect(fill = colrs[[2]], colour = "transparent"),
      strip.text = element_text(colour = 'white')
    )
)


# Load models/data
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

# Helper functions (won't need these once part of the pipeline)
source("R/facet_zoom2.R")
source("R/modelling-helpers.R")
source("R/plotting-helpers.R")
source("R/summarising-helpers.R")

# Load data
dat_wide_pre <- NMA_preDLI_datasets$wide
dat_long_pre <- NMA_preDLI_datasets$long

# Label endpoints
dat_long_pre[, "endpoint_lab" := factor(
  fcase(
    endpoint7_s == "cens", "Censored",
    endpoint7_s == "gvhd", "GvHD",
    endpoint7_s == "relapse", "Relapse",
    endpoint7_s == "other_nrf", "Other failure"
  ),
  levels = c("Censored", "GvHD", "Relapse", "Other failure")
)]

# Added fitted values and slope for CD3 model
dat_long_pre[, ':=' (
  curr_val_value = fitted(preDLI_JM_value_corr_CD3, type = "Subject"),
  curr_val_both = fitted(preDLI_JM_both_corr_CD3, type = "Subject"),
  slope = fitted_slopes_long(preDLI_JM_both_corr_CD3)
)]

# Make a data with last measurement
dat_long_pre_last <- dat_long_pre[, .SD[.N], by = "IDAA"]

# Sample 16 patients
set.seed(649846561)
IDAA_subs_pre <- sample(levels(dat_long_pre$IDAA), replace = FALSE, size = 16)


# pre-DLI: raw + indiv fits -----------------------------------------------


# First raw-plot
ggplot(dat_long_pre[IDAA %in% IDAA_subs_pre], aes(intSCT2_7, CD3_abs_log)) +
  geom_point(size = 3.5, alpha = 0.8, col = colrs[[1]]) +
  geom_vline(aes(xintercept = endpoint7), linetype = "dashed") +
  facet_wrap(~ as.numeric(factor(IDAA))) +  # instead of IDAA to get 1:16
  labs(x = "Time since alloSCT (months)", y = expression(paste("CD3 (x10"^"6","/l)"))) +
  log_axis_scales +
  geom_label(
    data = dat_long_pre_last[IDAA %in% IDAA_subs_pre],
    aes(x = endpoint7 + 0.05, y = log(2.5), label = endpoint_lab),
    hjust = 0,
    lineheight = .8,
    family = "Roboto Condensed",
    label.size = NA
  ) +
  geom_curve(
    data = dat_long_pre_last[IDAA %in% IDAA_subs_pre],
    mapping = aes(
      x = endpoint7 + 0.6,
      y = log(5),
      xend = endpoint7 + 0.05,
      yend = log(10)
    ),
    colour = "black",
    linewidth = 0.75,
    curvature = 0.3,
    arrow = arrow(length = unit(0.05, "npc"), type = "open"),
    inherit.aes = FALSE
  ) +
  coord_cartesian(xlim = c(0, 7.5)) +
  geom_line(aes(y = curr_val_value), linewidth = 1, col = colrs[[6]]) # indivdual fits

# Figure 1 (later to save as part of pipeline)
ggsave(
  here("analysis/figures/preDLI_lines_indiv.png"),
  dpi = 300,
  width = 12,
  height = 8
)

# all pre-DLI trajectories per endpoint -----------------------------------


# Do we exclude the censoring here?
ggplot(dat_long_pre, aes(intSCT2_7, CD3_abs_log, group = IDAA)) +
  geom_line(show.legend = FALSE, size = 1, alpha = 0.5, col = colrs[[1]])+
  labs(x = "Time since alloSCT (months)", y = expression(paste("CD3 (x10"^"6","/l)"))) +
  log_axis_scales +
  geom_smooth(
    linewidth = 2,
    aes(group = endpoint_lab),
    se = FALSE,
    method = "loess",
    formula = y ~ x,
    col = colrs[[6]]
  ) +
  facet_wrap(~ endpoint_lab) +
  coord_cartesian(ylim = log(c(0.1, 3000)))

# This is Supplemental Figure 1 (exclude or include the smooths?)
ggsave(
  here("analysis/figures/preDLI_perEndpoint.png"),
  dpi = 300,
  width = 8,
  height = 8
)


# pre-DLI: marginal fit ---------------------------------------------------


mod_names <- c(
  paste0("preDLI_JM_value_corr_CD", c(3, 4, 8)),
  paste0("preDLI_JM_both_corr_CD", c(3, 4, 8)),
  paste0("postDLI_JM_corr_CD", c(3, 4, 8))
)
names(mod_names) <- mod_names
mods <- lapply(mod_names, tar_read_raw)

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

# Good plot with one -> use this one to make cmv plot..
rbindlist(marg_preds_preDLI, idcol = "cell_line") |>
  ggplot(aes(x = intSCT2_7, y = pred, group = interaction(ATG, hirisk))) +
  geom_ribbon(aes(ymin = low, ymax = upp), fill = confint_col, alpha = confint_alpha, col = NA) +
  geom_line(aes(linetype = hirisk, col = ATG), size = 1.5) +
  facet_grid(
    facets = cell_line ~ CMV_PD,
    labeller = as_labeller(
      c(
        "-/-" = "CMV: -/-",
        "other P/D" = "CMV: other",
        "CD3" = "CD3",
        "CD4" = "CD4",
        "CD8" = "CD8"
      )
    )
  ) +
  labs(
    x = "Time since alloSCT (months)",
    y = expression(paste("cell count (x10"^"6","/l)")),
    linetype = "ITT",
    col = "Donor type"
  ) +
  log_axis_scales+
  coord_cartesian(xlim = c(0, 6), ylim = c(log(0.1), log(1500))) +
  scale_color_manual(
    labels = c("RD", "UD+ATG"),
    values = c("brown", "darkblue")
  ) +
  scale_linetype_manual(
    labels = c("No", "Yes (high risk)"),
    values = c("solid", "dotdash")
  ) +
  theme(legend.position = "bottom")
#ggsave("analysis/figures/preDLI_trajectories.png",dpi=300,width=6,height=8) # this is Figure 2


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
rbindlist(marg_preds_postDLI, idcol = "cell_line")[CMV_PD == "-/-"] |>
  ggplot(aes(x = intDLI1, y = pred, group = ATG)) + # was intDLI1 + 3
  geom_ribbon(aes(ymin = low, ymax = upp), fill = "gray", alpha = 0.5, col = NA) +
  geom_line(aes(col = ATG), size = 1.5, alpha = 0.75) +
  # Add repel labels?
  facet_grid(facets = cell_line ~ .
             # labeller = as_labeller(c("-/-"="CMV: -/-", "other P/D"="CMV: other",
             #                          "CD3"="CD3", "CD4"="CD4", "CD8"="CD8"
             #                          #"CD3"="total T-cells", "CD4"="CD4+ T-cells", "CD8"="CD8+ T-cells"
             #                          ))
             ) +
  labs(
    x = "Time since DLI (months)", # was alloSCT
    y = expression(paste("cell count (x10"^"6","/l)")),
    #linetype = "ITT",
    col = "Donor type"
  ) +
  scale_y_continuous(
    breaks = log(c(0.1, 1, 5, 25, 100, 500, 1500)),
    labels = c(0.1, 1, 5, 25, 100, 500, 1500)
  ) +
  coord_cartesian(xlim = c(0, 3), ylim = c(log(0.1), log(1500))) + # xlim was c(0,6)
  scale_color_manual(
    labels = c("RD", "UD+ATG"),
    values = c(colrs[[1]], colrs[[2]])
  ) +
  theme(legend.position = "bottom")
  #geom_vline(xintercept = 3, linetype = "dotted")
ggsave("analysis/figures/postDLI_trajectories.png",dpi=300,width=6,height=8) # this is FIgure 4


# NEW CMV FIGURE HERE: ----------------------------------------------------


rbindlist(marg_preds_postDLI, idcol = "cell_line") |>
  ggplot(aes(x = intDLI1, y = pred, group = CMV_PD, linetype = CMV_PD, col = CMV_PD)) +
  geom_line(linewidth = 1.5) +
  geom_ribbon(aes(ymin = low, ymax = upp), fill = "lightgray", alpha = 0.5, col = NA) +
  facet_grid(
    cell_line ~ ATG,
    labeller = as_labeller(c("UD" = "RD", "UD(+ATG)" = "UD(+ATG)", "CD3" = "CD3",
                             "CD4" = "CD4",
                             "CD8" = "CD8"))
  ) +
  labs(
    x = "Time since DLI (months)", # was alloSCT
    y = expression(paste("cell count (x10"^"6","/l)"))#,
    #linetype = "ITT",
    #col = "Donor type"
  ) +
  scale_y_continuous(
    breaks = log(c(0.1, 1, 5, 25, 100, 500, 1500)),
    labels = c(0.1, 1, 5, 25, 100, 500, 1500)
  ) +
  scale_color_manual(
    labels = c("-/-", "other P/D"),
    values = c(colrs[[1]], colrs[[2]])
  ) +
  guides(colour = guide_legend(reverse = TRUE), linetype = guide_legend(reverse = TRUE)) +
  #theme_cowplot() +
  #theme_light(base_size = 14) +
  coord_cartesian(xlim = c(0, 3), ylim = c(log(0.1), log(1500)))

# Pathetic for CMV
summary(preDLI_JM_value_corr_CD4)
summary(postDLI_JM_corr_CD4) # do it for post DLI


# Post DLI summaries ------------------------------------------------------


# -- Model summaries
ggplot(dat_wide_postDLI, aes(n_measurements)) +
  geom_histogram(bins = 13, col = "black", fill = "lightblue") +
  labs(x = "# Repeated measurements per patient")

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
ggplot(dat_long_postDLI, aes(intSCT2_7_reset, CD3_abs_log, group=IDAA)) + # CD3_abs_log instead of curr_val to get raw values
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
  mutate(
    endpoint=factor(
      ifelse(grepl("1",param),"GvHD",ifelse(grepl("2",param),"Relapse", "Other failure")),
      levels=c("GvHD","Relapse","Other failure")
    ),
    param2=factor(
      ifelse(grepl("ATG",param),"UD+ATG\n(vs. RD)",ifelse(grepl("hirisk",param),"High risk\n(vs. Non-high risk)", "Current value (t)")),
      levels=c("UD+ATG\n(vs. RD)","High risk\n(vs. Non-high risk)","Current value (t)")
    ),
    subset=factor(cell_line,levels=c("CD8","CD4","CD3")),
    label=paste0(ifelse(HR>=1,"  ", ""), round(HR,digits=2), " (", round(low,digits=2), "-", round(upp,digits=2), ")",
                 ifelse(HR<1, "   \n ", "\n "))
  )%>%
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
    size = 3,
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
        panel.grid.major.y=element_blank()) +
  scale_color_brewer(palette = "Dark2")
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
  mutate(
    endpoint=factor(ifelse(grepl("1",param),"GvHD","Relapse or other failure"),
                    levels=c("GvHD","Relapse or other failure")),
    param2=factor(ifelse(grepl("ATG",param),"UD+ATG\n(vs. RD)", "Current value (t)"),
                  levels=c("UD+ATG\n(vs. RD)", "Current value (t)")),
    subset=factor(cell_line,levels=c("CD8","CD4","CD3")),
    label=paste0(ifelse(HR>=1,"  ", ""), format(round(HR,digits=2)),
                 " (", round(low,digits=2), "-", round(upp,digits=2), ")",
                 ifelse(HR<1, "   \n \n ", "\n \n "))
  )%>%
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
    size = 3,
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
        panel.grid.major.y=element_blank(),plot.margin=margin(l=0)) +
  scale_color_brewer(palette = "Dark2")
ggsave("analysis/figures/postDLI_forest.png",dpi=300,height=6,width=11)



# Facet zoom customed -----------------------------------------------------

zoom_areas <- list(
  "CD4" = list("x" = c(2, 6), y = c(log(10), log(200))),
  "CD8" = list("x" = c(2, 6), y = c(log(25), log(500))),
  "CD3" = list("x" = c(2, 6), y = c(log(100), log(700)))
)

CD4_sub <- data.table(marg_preds_preDLI$CD4)[CMV_PD == "-/-"] |>
  ggplot(aes(x = intSCT2_7, y = pred, group = interaction(ATG, hirisk))) +
  geom_ribbon(
    aes(ymin = low, ymax = upp),
    fill = "lightgray",
    alpha = 0.5,
    col = NA
  ) +
  geom_line(aes(linetype = hirisk, col = ATG), size = 1.5) +
  labs(
    x = "Time since alloSCT (months)",
    y = expression(paste("CD4 cell count (x10"^"6","/l)")),
    linetype = "ITT",
    col = "Donor type"
  ) +
 # coord_cartesian(ylim = c(log()))
  scale_y_continuous(
    breaks = log(c(0.1, 1, 5, 25, 100, 500, 1500)),
    labels = c(0.1, 1, 5, 25, 100, 500, 1500)#,
    #limits = c(log(0.5), log(500))
  ) +
  #coord_fixed()
  #coord_cartesian(ylim = c(log(0.5), log(500))) +
  geom_rect(
    xmin = zoom_areas$CD4$x[1],
    xmax = zoom_areas$CD4$x[2],
    ymin = zoom_areas$CD4$y[1],
    ymax = zoom_areas$CD4$y[2],
    col = "black",
    linetype = "dashed",
    fill = NA
  ) +
  scale_color_manual(
    labels = c("-/-", "other P/D"),
    values = c(colrs[[2]], colrs[[1]])
  ) +
  facet_zoom2(
    xlim = zoom_areas$CD4$x,
    ylim = zoom_areas$CD4$y,
    show.area = FALSE,#shrink = T,
    zoom.size = 1L
  ) +
  theme(axis.title.x = element_blank())

CD8_sub <- data.table(marg_preds_preDLI$CD8)[CMV_PD == "-/-"] |>
  ggplot(aes(x = intSCT2_7, y = pred, group = interaction(ATG, hirisk))) +
  geom_ribbon(
    aes(ymin = low, ymax = upp),
    fill = "lightgray",
    alpha = 0.75,
    col = NA
  ) +
  geom_line(aes(linetype = hirisk, col = ATG), size = 1.5) +
  labs(
    x = "Time since alloSCT (months)",
    y = expression(paste("CD8 cell count (x10"^"6","/l)")),
    linetype = "ITT",
    col = "Donor type"
  ) +
  geom_rect(
    xmin = zoom_areas$CD8$x[1],
    xmax = zoom_areas$CD8$x[2],
    ymin = zoom_areas$CD8$y[1],
    ymax = zoom_areas$CD8$y[2],
    col = "black",
    linetype = "dashed",
    fill = NA
  ) +
  facet_zoom2(
    xlim = zoom_areas$CD8$x,
    ylim = zoom_areas$CD8$y,
    show.area = FALSE,
    zoom.size = 1L
  ) +
  scale_y_continuous(
    breaks = log(c(0.1, 1, 5, 25, 100, 500, 1500)),
    labels = c(0.1, 1, 5, 25, 100, 500, 1500)
  ) +
  theme(axis.title.x = element_blank())

#CD8_sub

# https://stackoverflow.com/questions/45221783/ggforce-facet-zoom-labels-only-on-zoomed-example/49927431#49927431

CD3_sub <- data.table(marg_preds_preDLI$CD3)[CMV_PD == "-/-"] |>
  ggplot(aes(x = intSCT2_7, y = pred, group = interaction(ATG, hirisk))) +
  geom_ribbon(
    aes(ymin = low, ymax = upp),
    fill = "lightgray",
    alpha = 0.75,
    col = NA
  ) +
  geom_line(aes(linetype = hirisk, col = ATG), size = 1.5) +
  labs(
    x = "Time since alloSCT (months)",
    y = expression(paste("CD3 cell count (x10"^"6","/l)")),
    linetype = "ITT",
    col = "Donor type"
  ) +
  geom_rect(
    xmin = zoom_areas$CD3$x[1],
    xmax = zoom_areas$CD3$x[2],
    ymin = zoom_areas$CD3$y[1],
    ymax = zoom_areas$CD3$y[2],
    col = "black",
    linetype = "dashed",
    fill = NA
  ) +
  facet_zoom2(
    xlim = zoom_areas$CD3$x,
    ylim = zoom_areas$CD3$y,
    show.area = FALSE,
    zoom.size = 1L
  ) +
  scale_y_continuous(
    breaks = log(c(0.1, 1, 5, 25, 100, 500, 1500)),
    labels = c(0.1, 1, 5, 25, 100, 500, 1500)
  ) +
  theme(axis.title.x = element_blank())

remove_zoom_rectangle <- function(p) {
  pb <- ggplot_build(p)
  pb$data[[3]][pb$data[[3]]$PANEL == 4, "colour"] <- NA
  pg <- ggplot_gtable(pb)
  ggplotify::as.ggplot(pg)
}

legend_b <- get_legend(CD4_sub + theme(legend.position = "top") +
                         scale_colour_discrete(labels = c("RD", "UD(+ATG)")) +
                         scale_linetype_discrete("Risk group",
                                                 labels = c("Non-high risk", "High risk"))
                         )

prow <- plot_grid(
  remove_zoom_rectangle(CD3_sub + theme(legend.position="none")),
  remove_zoom_rectangle(CD4_sub + theme(legend.position="none")),
  remove_zoom_rectangle(CD8_sub + theme(legend.position="none")),
  #align = 'vh',
  #labels = c("A", "B", "C"),
  #hjust = -1,
  nrow = 3
)
prow_leg <- plot_grid(
  legend_b, prow, ncol = 1, rel_heights = c(.1, 1)
)

final <- ggpubr::annotate_figure(
  prow_leg,
  bottom = ggpubr::text_grob("Time since alloSCT (months)",
                             hjust = 0.5)#, size = 10)
)

final
class(final)

#pdf(file = "figures/marginals-zoomed.pdf", width = 6, height = 8)
#final
#dev.off()
ggsave(
  "./analysis/figures/marginals-zoomed.png",
  plot = final,
  width = 6,
  height = 8,
  units = "in",
  bg = "transparent"
)


png(
  filename = "./analysis/figures/marginals-zoomed.png",
  width = 6, height = 8,
  units = "in", res = 300
)
final
dev.off()
#gridExtra::grid.arrange(
#  CD4_sub, CD3_sub, CD8_sub
#)


# New facet zoom, mirrored.. ----------------------------------------------


# Stack link here:
