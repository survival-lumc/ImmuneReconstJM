list(
  tar_target(
    lymphocytes_raw,
    data.table(readRDS(here("data-raw/2022-02-24_v9/lymphocytes.rds")))
  ),
  tar_target(
    variables_raw,
    data.table(readRDS(here("data-raw/2022-02-24_v9/variables.rds")))
  ),
  tar_target(dat_merged, prepare_raw_data(lymphocytes_raw, variables_raw))
)
