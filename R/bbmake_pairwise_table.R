#' Create pairwise comparison table (preserves by= grouping + always cleans rowid)
#'
#' @param pw An `emmGrid` object (output from `pairs(emm)`)
#' @param model Original fitted model (optional). Recovered from `pw` when
#'   possible; used for Gaussian Cohen's d.
#' @param adjust Optional within-`by`-group multiplicity adjustment passed
#'   to [emmeans::summary.emmGrid()]. `NULL` (the default) keeps the
#'   adjustment already stored on `pw` (e.g. from `pairs(..., adjust = )`).
#'   Ignored when `cross.adjust` is not `"none"`.
#' @param cross.adjust If not `"none"`, treat **every** pairwise test in
#'   `pw` as one family: `summary(pw, by = NULL, adjust = cross.adjust)`.
#'   That is the total-number-of-comparisons correction (e.g. 4 Sex × Diet
#'   cells, or 6 tests from 3-level Stim × 2 Groups). Default `"none"`.
#'   This is *not* emmeans' own `cross.adjust` argument, which only
#'   adjusts matching contrasts across by-groups.
#'
#' @return A tibble with grouping columns first, then contrast + effect sizes
#' @export
bbmake_pairwise_table <- function(
  pw,
  model = NULL,
  adjust = NULL,
  cross.adjust = "none"
) {

  # Auto-recover model if not supplied
  if (is.null(model)) {
    model <- .bb_recover_model(pw)
    if (is.null(model)) {
      model <- tryCatch(
        get("model", envir = parent.frame(), inherits = TRUE),
        error = function(e) NULL
      )
    }
  }

  # Detect by-grouping variables (Sex, Diet, Satiety, etc.)
  by_vars <- if (!is.null(pw@misc$by.vars)) pw@misc$by.vars else character(0)

  pw_summary <- tibble::as_tibble(
    .bb_pairs_summary(pw, adjust = adjust, cross.adjust = cross.adjust)
  )

  # Robust CI column detection
  lcl_col <- ifelse("lower.CL" %in% names(pw_summary), "lower.CL", "asymp.LCL")
  ucl_col <- ifelse("upper.CL" %in% names(pw_summary), "upper.CL", "asymp.UCL")

  has_ratio <- "ratio"      %in% colnames(pw_summary)
  has_odds  <- "odds.ratio" %in% colnames(pw_summary)

  if (has_ratio) {
    pw_table <- pw_summary |>
      mutate(
        IRR        = round(ratio, 2),
        lower.CL   = .data[[lcl_col]],
        upper.CL   = .data[[ucl_col]],
        `% Change` = 100 * (ratio - 1)
      ) |>
      select(any_of(by_vars), contrast,
             any_of(c("z.ratio", "t.ratio")), p.value,
             IRR, lower.CL, upper.CL)

  } else if (has_odds) {
    pw_table <- pw_summary |>
      mutate(
        `Odds Ratio` = round(odds.ratio, 2),
        lower.CL   = .data[[lcl_col]],
        upper.CL   = .data[[ucl_col]]
      ) |>
      select(any_of(by_vars), contrast,
             any_of(c("z.ratio", "t.ratio")), p.value,
             `Odds Ratio`, lower.CL, upper.CL)

  } else {
    # ==================== GAUSSIAN / LINK-SCALE MODELS ====================
    pw_table <- pw_summary |>
      select(any_of(by_vars), everything()) |>
      rename(`Mean Difference` = estimate) |>
      select(
        any_of(by_vars),
        contrast,
        any_of(c("t.ratio", "z.ratio")),
        any_of("df"),
        p.value
      ) |>
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
          rename(`d` = effect.size) |>
          # mutate(`d 95% CI` = sprintf("[%.2f, %.2f]", lower.CL, upper.CL)) |>
          select(rowid, d, lower.CL, upper.CL)
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
    mutate(contrast = stringr::str_trim(contrast)) |>
    arrange(across(any_of(by_vars)), p.value)
}

#' @keywords internal
#' @noRd
.bb_pairs_summary <- function(
  pw,
  adjust = NULL,
  cross.adjust = "none",
  infer = c(TRUE, TRUE)
) {
  pool <- !is.null(cross.adjust) && !identical(cross.adjust, "none")
  if (isTRUE(pool)) {
    return(summary(pw, infer = infer, by = NULL, adjust = cross.adjust))
  }
  if (is.null(adjust)) {
    summary(pw, infer = infer)
  } else {
    summary(pw, infer = infer, adjust = adjust)
  }
}

#' @keywords internal
.bb_recover_model <- function(...) {
  classes <- c("lm", "glm", "glmmTMB", "merMod", "lmerModLmerTest", "glmerMod")
  for (obj in list(...)) {
    if (is.null(obj) || !inherits(obj, "emmGrid")) {
      next
    }
    m <- tryCatch(obj@model.info$object, error = function(e) NULL)
    if (inherits(m, classes)) {
      return(m)
    }
  }
  NULL
}
