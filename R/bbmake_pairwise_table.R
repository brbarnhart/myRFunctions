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
#' @return A tibble with grouping columns (when present), then `contrast`,
#'   the test statistic (`z.ratio` or `t.ratio`), `df` (when present), the
#'   effect (`IRR`, `Odds Ratio`, or `Mean Difference`), `SE`, `lower.CL`,
#'   `upper.CL`, and `p.value`. Gaussian tables may also include Cohen's `d`.
#'   Values are left at full precision for downstream formatting.
#' @export
bbmake_pairwise_table <- function(
  pw,
  model = NULL,
  adjust = NULL,
  cross.adjust = "none"
) {
  if (is.null(model)) {
    model <- .bb_recover_model(pw)
  }

  by_vars <- if (!is.null(pw@misc$by.vars)) pw@misc$by.vars else character(0)

  pw_summary <- tibble::as_tibble(
    .bb_pairs_summary(pw, adjust = adjust, cross.adjust = cross.adjust)
  )

  out <- .bb_tidy_pairwise_table(pw_summary, by_vars = by_vars)

  if (!is.null(model) && "Mean Difference" %in% names(out)) {
    cohen_d <- tryCatch(
      {
        d_tab <- emmeans::eff_size(
          pw,
          sigma = stats::sigma(model),
          edf = stats::df.residual(model),
          method = "identity"
        ) |>
          summary(infer = TRUE) |>
          tibble::as_tibble()
        tibble::tibble(d = d_tab$effect.size)
      },
      error = function(e) {
        warning("Cohen's d could not be calculated: ", e$message, call. = FALSE)
        NULL
      }
    )
    if (!is.null(cohen_d) && nrow(cohen_d) == nrow(out)) {
      out$d <- cohen_d$d
    }
  }

  dplyr::arrange(
    out,
    dplyr::across(dplyr::any_of(c(by_vars, "p.value")))
  )
}

#' @keywords internal
#' @noRd
.bb_tidy_pairwise_table <- function(pw_summary, by_vars) {
  lcl_col <- if ("lower.CL" %in% names(pw_summary)) {
    "lower.CL"
  } else if ("asymp.LCL" %in% names(pw_summary)) {
    "asymp.LCL"
  } else {
    NULL
  }
  ucl_col <- if ("upper.CL" %in% names(pw_summary)) {
    "upper.CL"
  } else if ("asymp.UCL" %in% names(pw_summary)) {
    "asymp.UCL"
  } else {
    NULL
  }

  if (!is.null(lcl_col) && lcl_col != "lower.CL") {
    names(pw_summary)[names(pw_summary) == lcl_col] <- "lower.CL"
  }
  if (!is.null(ucl_col) && ucl_col != "upper.CL") {
    names(pw_summary)[names(pw_summary) == ucl_col] <- "upper.CL"
  }

  if ("ratio" %in% names(pw_summary)) {
    pw_summary <- dplyr::rename(pw_summary, IRR = "ratio")
  } else if ("odds.ratio" %in% names(pw_summary)) {
    pw_summary <- dplyr::rename(pw_summary, `Odds Ratio` = "odds.ratio")
  } else if ("estimate" %in% names(pw_summary)) {
    pw_summary <- dplyr::rename(pw_summary, `Mean Difference` = "estimate")
  }

  if (!"contrast" %in% names(pw_summary) || !"p.value" %in% names(pw_summary)) {
    stop(
      "pairs() did not return contrast and p.value. Columns were: ",
      paste(names(pw_summary), collapse = ", "),
      call. = FALSE
    )
  }

  effect_cols <- intersect(
    c("IRR", "Odds Ratio", "Mean Difference"),
    names(pw_summary)
  )
  if (length(effect_cols) == 0L) {
    stop(
      "pairs() did not return ratio, odds.ratio, or estimate. Columns were: ",
      paste(names(pw_summary), collapse = ", "),
      call. = FALSE
    )
  }

  pw_summary[["contrast"]] <- stringr::str_trim(pw_summary[["contrast"]])

  keep <- c(
    by_vars,
    "contrast",
    "df",
    "z.ratio",
    "t.ratio",
    effect_cols,
    "lower.CL",
    "upper.CL",
    "SE",
    "p.value"
  )
  keep <- keep[keep %in% names(pw_summary)]
  dplyr::select(pw_summary, dplyr::all_of(keep))
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
