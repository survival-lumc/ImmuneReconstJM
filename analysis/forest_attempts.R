tar_load(
  c(
    preDLI_JM_both_corr_CD3,
    preDLI_JM_both_corr_CD4,
    preDLI_JM_both_corr_CD8
  )
)


summ_CD3 <- summary(preDLI_JM_both_corr_CD3)
summ_CD3$`CoefTable-Event`
summ_CD3$D

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
