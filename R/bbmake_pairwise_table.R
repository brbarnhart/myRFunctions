#' Create pairwise comparison table (preserves by= grouping + always cleans rowid)
#'
#' @param pw An `emmGrid` object (output from `pairs(emm)`)
#' @param model Original fitted model (optional)
#'
#' @return A tibble with grouping columns first, then contrast + effect sizes
#' @export
bbmake_pairwise_table <- function(pw, model = NULL) {

  # Auto-recover model if not supplied
  if (is.null(model)) {
    model <- tryCatch(pw@model.info$object, error = function(e) NULL)
    if (is.null(model)) {
      model <- tryCatch(get("model", envir = parent.frame(), inherits = TRUE), error = function(e) NULL)
    }
  }

  # Detect by-grouping variables (Sex, Diet, Satiety, etc.)
  by_vars <- if (!is.null(pw@misc$by.vars)) pw@misc$by.vars else character(0)

  # Full summary with CIs
  pw_summary <- summary(pw, infer = c(TRUE, TRUE)) |>
    as_tibble()

  # Robust CI column detection
  lcl_col <- ifelse("lower.CL" %in% names(pw_summary), "lower.CL", "asymp.LCL")
  ucl_col <- ifelse("upper.CL" %in% names(pw_summary), "upper.CL", "asymp.UCL")

  has_ratio <- "ratio"      %in% colnames(pw_summary)
  has_odds  <- "odds.ratio" %in% colnames(pw_summary)

  if (has_ratio) {
    pw_table <- pw_summary |>
      mutate(
        `Rate Ratio` = round(ratio, 2),
        `RR 95% CI`  = sprintf("[%.2f, %.2f]", .data[[lcl_col]], .data[[ucl_col]]),
        `% Change`   = sprintf("%+.1f%%", 100 * (ratio - 1))
      ) |>
      select(any_of(by_vars), contrast, `Rate Ratio`, `RR 95% CI`, `% Change`,
             SE, any_of(c("z.ratio", "t.ratio")), p.value)

  } else if (has_odds) {
    pw_table <- pw_summary |>
      mutate(
        `Odds Ratio` = round(odds.ratio, 2),
        `OR 95% CI`  = sprintf("[%.2f, %.2f]", .data[[lcl_col]], .data[[ucl_col]])
      ) |>
      select(any_of(by_vars), contrast, `Odds Ratio`, `OR 95% CI`,
             SE, any_of(c("z.ratio", "t.ratio")), p.value)

  } else {
    # ==================== GAUSSIAN / LINK-SCALE MODELS ====================
    pw_table <- pw_summary |>
      select(any_of(by_vars), everything()) |>
      rename(`Mean Difference` = estimate) |>
      select(any_of(by_vars), contrast, `Mean Difference`, SE, any_of("df"),
             any_of(c("t.ratio", "z.ratio")), p.value) |>
      rowid_to_column("rowid")

    # Cohen's d (only for true Gaussian models)
    if (!is.null(model)) {
      cohen_d <- tryCatch({
        emmeans::eff_size(
          pw,
          sigma = sigma(model),
          edf   = df.residual(model),
          method = "identity"
        ) |>
          summary(infer = TRUE) |>
          as_tibble() |>
          rowid_to_column("rowid") |>
          rename(`Cohen's d` = effect.size) |>
          mutate(`d 95% CI` = sprintf("[%.2f, %.2f]", lower.CL, upper.CL)) |>
          select(rowid, `Cohen's d`, `d 95% CI`)
      }, error = function(e) {
        warning("Cohen's d could not be calculated: ", e$message, call. = FALSE)
        NULL
      })

      if (!is.null(cohen_d)) {
        pw_table <- left_join(pw_table, cohen_d, by = "rowid") |>
          select(-rowid)
      } else {
        pw_table <- pw_table |> select(-rowid)
      }
    } else {
      pw_table <- pw_table |> select(-rowid)
    }
  }

  # Final polishing
  pw_table |>
    mutate(contrast = str_trim(contrast)) |>
    arrange(across(any_of(by_vars)), p.value)
}
