# (Edit this later..)

#' Summarises event submodel of joint model
#'
#' Maybe later develop into tidy(..., model = "event") or so?
#'
#' (Only works when no interaction in functional form)
summarise_event_submodel <- function(JMfit,
                                     digits = 3) {

  summary_raw <- data.table::data.table(
    round(summary(JMfit)$`CoefTable-Event`, digits),
    keep.rownames = TRUE
  )

  # Exclude baseline hazard spline coefficients
  data.table::setnames(summary_raw, "rn", "Coefficient")
  dt <- summary_raw[grep("^bs", x = Coefficient, invert = TRUE)]

  # Set up
  dt[, event_num := regmatches(
    x = Coefficient,
    m = regexpr(pattern = "[+1-9]$", text = Coefficient)
  )]

  dt[, event := factor(
    x = event_num,
    levels = seq_len(4),
    labels = c("Cell. intervention", "Relapse", "NRF: Other", "NRF: GVHD")
  )]

  dt[, Coefficient := gsub(x = Coefficient, pattern = "\\:.*", replacement = "")]
  dt[, Coefficient := gsub(x = Coefficient, pattern = "\\.[+1-9]$", replacement = "")]
  dt[, Coefficient := factor(
    x = Coefficient,
    levels = c("ATG", "hirisk", "SCTyear_2010", "Assoct"),
    labels = c("ATG", "Disease risk: high", "SCT pre-2010", "Assoct.")
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

  return(dt_edit)
  # dt_edit %>%
  #   kableExtra::kbl(format = "html") %>%
  #   kable_paper("striped", full_width = F) %>%
  #   pack_rows("Cell. intervention", 1, 3) %>%
  #   pack_rows("Relapse", 4, 6) %>%
  #   pack_rows("NRF: Other", 7, 7) %>%
  #   pack_rows("NRF: GVHD", 8, 8)
}
