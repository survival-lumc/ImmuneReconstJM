list(
  tar_target(
    NMA_preDLI_datasets,
    get_preDLI_datasets(
      dat_merged[TCD2 %in% c("NMA RD: ALT", "UD: ALT + ATG")], 
      admin_cens = 6
    )
  ),
  tar_target(
    NMA_postDLI_datasets,
    get_postDLI_datasets(
      dat_merged = dat_merged[TCD2 %in% c("NMA RD: ALT", "UD: ALT + ATG")],
      admin_cens_dli = 3
    )
  )
)
