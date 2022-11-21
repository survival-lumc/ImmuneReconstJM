tar_target(reference_values, {
  cbind.data.frame(
    "cell_type" = c("CD3", "CD4", "CD8", "NK", "CD19"),
    "lower_limit" = c(860, 560, 260, 40, 60),
    "upper_limit" = c(2490, 1490, 990, 390, 1000)
  )
})
