##**********************************##
## Script to process raw data files ##
##**********************************##

#' Cleans raw lymphocytes and variables data
prepare_raw_data <- function(lymphocytes, variables) {

  # Merge into one
  merged_dat <- data.table::merge.data.table(x = lymphocytes, y = variables)

  # Drop unused levels
  factors <- which(vapply(merged_dat, FUN = is.factor, FUN.VALUE = logical(1)))
  merged_dat[, (factors) := lapply(.SD, droplevels), .SDcols = factors]

  # Make ID and hirisk into factor
  merged_dat[, ':=' (
    IDAA = as.factor(IDAA),
    hirisk = factor(hirisk, levels = c(0, 1), labels = c("no", "yes"))
  )]

  # Replace 0 (non-detectable) values with 0.5
  cell_count_vars <- grep(x = colnames(merged_dat), pattern = "_abs$", value = TRUE)
  merged_dat[, (cell_count_vars) := lapply(
    .SD, function(col) ifelse(col == 0, 0.5, col)
  ), .SDcols = cell_count_vars]

  # Add log versions of cell counts
  merged_dat[, paste0(cell_count_vars, "_log") := lapply(.SD, log), .SDcols = cell_count_vars]

  # Add variable for transplant pre 2010 SCT
  merged_dat[, SCTyear_2010 := factor(
    as.numeric(SCTyear < 2010), levels = c(0, 1), labels = c("post_2010", "pre_2010")
  )]

  return(merged_dat)
}
