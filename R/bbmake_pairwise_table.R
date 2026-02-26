#' Create pairwise comparison table (now supports LM/LMM + GLM/GLMM + type = "response")
#'
#' Automatically chooses the right effect size:
#' • Gaussian (lm / lmer): Mean Difference + Cohen's d [95% CI]
#' • Count / rate models (glmmTMB negative binomial / Poisson with `type = "response"`): Rate Ratio + % Change + 95% CI
#' • Binary / logistic: Odds Ratio + 95% CI
#'
#' @param pw An `emmGrid` object (output from `pairs(emm)`)
#' @param model Original fitted model (optional – usually auto-detected)
#'
#' @return A tibble
#' @export
bbmake_pairwise_table <- function(pw, model = NULL) {

  # Auto-recover model if not supplied
  if (is.null(model)) {
    model <- tryCatch(pw@model.info$object, error = function(e) NULL)
    if (is.null(model)) {
      model <- tryCatch(get("model", envir = parent.frame(), inherits = TRUE), error = function(e) NULL)
    }
  }

  # Full summary with CIs (infer = TRUE guarantees confidence bounds)
  pw_summary <- summary(pw, infer = c(TRUE, TRUE)) |>
    as_tibble()

  # Robust CI column detection (handles type = "link" vs type = "response")
  lcl_col <- ifelse("lower.CL" %in% names(pw_summary), "lower.CL", "asymp.LCL")
  ucl_col <- ifelse("upper.CL" %in% names(pw_summary), "upper.CL", "asymp.UCL")

  has_ratio <- "ratio"      %in% colnames(pw_summary)
  has_odds  <- "odds.ratio" %in% colnames(pw_summary)

  if (has_ratio) {
    # ==================== COUNT MODELS (NB / Poisson GLMM) ====================
    pw_table <- pw_summary |>
      mutate(
        `Rate Ratio` = round(ratio, 2),
        `RR 95% CI`  = sprintf("[%.2f, %.2f]", .data[[lcl_col]], .data[[ucl_col]]),
        `% Change`   = sprintf("%+.1f%%", 100 * (ratio - 1))
      ) |>
      select(contrast, `Rate Ratio`, `RR 95% CI`, `% Change`,
             SE, any_of(c("z.ratio", "t.ratio")), p.value)

  } else if (has_odds) {
    # ==================== BINARY / LOGISTIC ====================
    pw_table <- pw_summary |>
      mutate(
        `Odds Ratio` = round(odds.ratio, 2),
        `OR 95% CI`  = sprintf("[%.2f, %.2f]", .data[[lcl_col]], .data[[ucl_col]])
      ) |>
      select(contrast, `Odds Ratio`, `OR 95% CI`,
             SE, any_of(c("z.ratio", "t.ratio")), p.value)

  } else {
    # ==================== GAUSSIAN MODELS (lm / lmer) ====================
    pw_table <- pw_summary |>
      rename(`Mean Difference` = estimate) |>
      select(contrast, `Mean Difference`, SE, any_of("df"),
             any_of(c("t.ratio", "z.ratio")), p.value) |>
      rowid_to_column("rowid")

    # Cohen's d (only for Gaussian)
    if (!is.null(model)) {
      tryCatch({
        cohen_d <- emmeans::eff_size(
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

        pw_table <- left_join(pw_table, cohen_d, by = "rowid") |>
          select(-rowid)
      }, error = function(e) {
        warning("Cohen's d could not be calculated: ", e$message, call. = FALSE)
      })
    }
  }

  # Final polishing
  pw_table |>
    mutate(contrast = str_trim(contrast)) |>
    arrange(p.value)
}
