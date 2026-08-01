#' Normalize an emmeans grid into a plot-ready data frame
#'
#' Converts an `emmGrid` (or its data-frame summary) into a tibble with
#' consistent column names for plotting: `y` (estimate), `ymin` / `ymax`
#' (confidence limits), plus all factor columns from the grid.
#'
#' @param emm An `emmGrid` from [emmeans::emmeans()], or a data frame already
#'   produced by `as.data.frame(emm)` / `summary(emm)`.
#'
#' @return A tibble with columns `y`, `ymin`, `ymax`, and the original factor
#'   columns. Attributes:
#'   \describe{
#'     \item{`pri.vars`}{Primary (x-axis candidate) factor names}
#'     \item{`by.vars`}{By-variable names used for faceting}
#'     \item{`mean_col`}{Original mean column name (`emmean` or `response`)}
#'     \item{`lower_col`, `upper_col`}{Original CI column names}
#'   }
#' @export
#' @examples
#' \dontrun{
#' emm <- emmeans(mod, ~ Stim | Diet * Sex, type = "response")
#' df  <- bb_emm_df(emm)
#' ggplot(df, aes(x = Stim, y = y, ymin = ymin, ymax = ymax)) +
#'   geom_pointrange() +
#'   facet_grid(Diet ~ Sex)
#' }
bb_emm_df <- function(emm) {
  if (inherits(emm, "emmGrid")) {
    pri_vars <- emm@misc$pri.vars %||% character(0)
    by_vars  <- emm@misc$by.vars  %||% character(0)
    emm_df   <- as.data.frame(emm)
  } else if (is.data.frame(emm)) {
    pri_vars <- character(0)
    by_vars  <- character(0)
    emm_df   <- as.data.frame(emm)
  } else {
    stop("`emm` must be an emmGrid or a data frame.", call. = FALSE)
  }

  mean_col <- .bb_detect_mean_col(emm_df)
  ci_cols  <- .bb_detect_ci_cols(emm_df)

  out <- emm_df |>
    dplyr::mutate(
      y    = .data[[mean_col]],
      ymin = .data[[ci_cols$lower]],
      ymax = .data[[ci_cols$upper]]
    ) |>
    tibble::as_tibble()

  # If pri/by not available (data-frame input), try to recover factor columns
  if (length(pri_vars) == 0L && length(by_vars) == 0L) {
    meta_cols <- c(
      mean_col, ci_cols$lower, ci_cols$upper,
      "y", "ymin", "ymax", "SE", "df", "t.ratio", "z.ratio",
      "p.value", "asymp.LCL", "asymp.UCL", "lower.CL", "upper.CL",
      "emmean", "response", "rate", "prob", "null"
    )
    factor_cols <- setdiff(names(out), meta_cols)
    pri_vars <- factor_cols
  }

  attr(out, "pri.vars")   <- pri_vars
  attr(out, "by.vars")    <- by_vars
  attr(out, "mean_col")   <- mean_col
  attr(out, "lower_col")  <- ci_cols$lower
  attr(out, "upper_col")  <- ci_cols$upper
  out
}

# ── Internal helpers ──────────────────────────────────────────────────────────

#' @keywords internal
.bb_detect_mean_col <- function(df) {
  candidates <- c("response", "emmean", "rate", "prob")
  hit <- candidates[candidates %in% names(df)]
  if (length(hit) == 0L) {
    stop(
      "No mean column found. Expected one of: ",
      paste(candidates, collapse = ", "),
      call. = FALSE
    )
  }
  hit[[1]]
}

#' @keywords internal
.bb_detect_ci_cols <- function(df) {
  lower <- if ("lower.CL" %in% names(df)) {
    "lower.CL"
  } else if ("asymp.LCL" %in% names(df)) {
    "asymp.LCL"
  } else {
    stop("No lower CI column (lower.CL or asymp.LCL) found.", call. = FALSE)
  }
  upper <- if ("upper.CL" %in% names(df)) {
    "upper.CL"
  } else if ("asymp.UCL" %in% names(df)) {
    "asymp.UCL"
  } else {
    stop("No upper CI column (upper.CL or asymp.UCL) found.", call. = FALSE)
  }
  list(lower = lower, upper = upper)
}

#' @keywords internal
`%||%` <- function(x, y) {
  if (is.null(x)) y else x
}
