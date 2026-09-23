#' Format a pairwise comparison table
#'
#' Builds a [flextable::flextable()] from the tibble returned by
#' [bbmake_pairwise_table()] or [bbmake_pairwise_sensitivity_table()].
#' Value formatting follows [rempsyc::nice_table()]: estimates use two
#' decimals, `p.value` is drawn with [rempsyc::format_p()] (no leading
#' zero, `< .001` below that cutoff, significance stars), and the
#' confidence limits become one `[lower, upper]` column. That interval
#' stays where `lower.CL` and `upper.CL` were. Every other column stays
#' in the order those functions returned it.
#'
#' Headers are relabelled for the paper table: `contrast` becomes
#' `Contrast`, `z.ratio` and `t.ratio` become italic *z* and *t*, and
#' `p.value` becomes italic *p*. `SE`, `df`, and Cohen's `d` keep their
#' names and are italicised when those columns are still present.
#' Columns to the left of `Contrast` (the emmeans by-variables) are
#' merged vertically, then [flextable::theme_vanilla()] sets the borders,
#' bold header, and alignment. The result is a flextable, so
#' [flextable::save_as_docx()] and further flextable edits still apply.
#' The interval header is `95% CI`, which is emmeans' default level.
#'
#' @param data A data frame from [bbmake_pairwise_table()] or
#'   [bbmake_pairwise_sensitivity_table()].
#' @param drop Columns to remove before formatting, matched to the input
#'   names. Default `c("SE", "df")`. `NULL` or `character()` keeps every
#'   column. Names that are not in `data` are ignored.
#' @param merge Columns to merge vertically where the same value repeats
#'   in consecutive rows ([flextable::merge_v()]). `TRUE` (default) merges
#'   every column to the left of `contrast` (the emmeans by-variables:
#'   `Sex`, `Diet`, and so on). `FALSE` does not merge. A character vector
#'   names specific columns, using either the input names (`"contrast"`)
#'   or the printed names (`"Contrast"`).
#' @param stars If `TRUE` (default), append significance stars to *p*
#'   values (`*` `< .05`, `**` `< .01`, `***` `< .001`).
#'
#' @return A `flextable`. Grouping columns, `Contrast`, `Model`, the test
#'   statistic, the effect, `95% CI`, and `p` stay in the input order,
#'   with `SE` and `df` removed when `drop` says so.
#'
#' @seealso [bbmake_pairwise_table()], [bbmake_pairwise_sensitivity_table()],
#'   [rempsyc::nice_table()], [flextable::theme_vanilla()]
#' @export
#' @examples
#' set.seed(1)
#' dat <- expand.grid(
#'   Sex = factor(c("F", "M")),
#'   Stim = factor(c("A", "B")),
#'   id = 1:8
#' )
#' dat$y <- rpois(nrow(dat), lambda = 5)
#' mods <- list(
#'   Full = glm(y ~ Sex * Stim, data = dat, family = poisson),
#'   Reduced = glm(
#'     y ~ Sex * Stim,
#'     data = dat[dat$id != "1", ],
#'     family = poisson
#'   )
#' )
#' bbmake_pairwise_sensitivity_table(
#'   mods,
#'   ~ Stim | Sex,
#'   adjust = "none",
#'   cross.adjust = "holm"
#' ) |>
#'   bbnice_pairwise_table()
bbnice_pairwise_table <- function(
  data,
  drop = c("SE", "df"),
  merge = TRUE,
  stars = TRUE
) {
  data <- .bb_as_pairwise_table(data)
  if (is.null(drop)) {
    drop <- character()
  } else if (!is.character(drop)) {
    stop(
      "`drop` must be a character vector of column names, or NULL.",
      call. = FALSE
    )
  }
  if (!is.logical(stars) || length(stars) != 1L || is.na(stars)) {
    stop("`stars` must be TRUE or FALSE.", call. = FALSE)
  }

  data <- dplyr::select(data, -dplyr::any_of(drop))
  data <- .bb_collapse_ci(data)
  data <- .bb_rename_pairwise_headers(data)
  merge_cols <- .bb_merge_columns(merge, names(data))
  .bb_flextable_pairwise(data, merge_cols = merge_cols, stars = stars)
}

#' @keywords internal
#' @noRd
.bb_as_pairwise_table <- function(data) {
  if (!is.data.frame(data)) {
    stop(
      "`data` must be a data frame from bbmake_pairwise_table() or ",
      "bbmake_pairwise_sensitivity_table().",
      call. = FALSE
    )
  }
  data <- tibble::as_tibble(data)
  missing <- setdiff(c("contrast", "p.value"), names(data))
  if (length(missing)) {
    stop(
      "Expected a pairwise table with columns contrast and p.value. Missing: ",
      paste(missing, collapse = ", "),
      ". Columns were: ",
      paste(names(data), collapse = ", "),
      call. = FALSE
    )
  }
  data
}

#' @keywords internal
#' @noRd
.bb_collapse_ci <- function(data) {
  nm <- names(data)
  if (!all(c("lower.CL", "upper.CL") %in% nm)) {
    return(data)
  }
  # One interval column where the two limits started, so later columns
  # (usually p.value) stay to the right of the interval.
  pos <- min(match(c("lower.CL", "upper.CL"), nm))
  ci <- .bb_format_ci(data[["lower.CL"]], data[["upper.CL"]])
  data[["lower.CL"]] <- NULL
  data[["upper.CL"]] <- NULL
  tibble::add_column(data, `95% CI` = ci, .before = pos)
}

#' @keywords internal
#' @noRd
.bb_format_ci <- function(lower, upper, digits = 2L) {
  fmt <- function(x) {
    formatC(round(as.numeric(x), digits), digits, format = "f")
  }
  out <- paste0("[", fmt(lower), ", ", fmt(upper), "]")
  out[is.na(lower) | is.na(upper)] <- ""
  out
}

#' @keywords internal
#' @noRd
.bb_rename_pairwise_headers <- function(data) {
  map <- c(
    contrast = "Contrast",
    z.ratio = "z",
    t.ratio = "t",
    p.value = "p"
  )
  nm <- names(data)
  idx <- match(names(map), nm)
  ok <- !is.na(idx)
  nm[idx[ok]] <- unname(map[ok])
  if (anyDuplicated(nm)) {
    stop(
      "Header formatting produced duplicate column names: ",
      paste(unique(nm[duplicated(nm)]), collapse = ", "),
      call. = FALSE
    )
  }
  names(data) <- nm
  data
}

#' @keywords internal
#' @noRd
.bb_merge_columns <- function(merge, names_now) {
  if (isFALSE(merge)) {
    return(character())
  }
  aliases <- c(
    contrast = "Contrast",
    z.ratio = "z",
    t.ratio = "t",
    p.value = "p",
    lower.CL = "95% CI",
    upper.CL = "95% CI"
  )
  if (isTRUE(merge)) {
    pos <- match("Contrast", names_now)
    if (is.na(pos) || pos <= 1L) {
      return(character())
    }
    return(names_now[seq_len(pos - 1L)])
  }
  if (
    !is.character(merge) ||
      !length(merge) ||
      anyNA(merge) ||
      any(!nzchar(merge))
  ) {
    stop(
      "`merge` must be TRUE, FALSE, or a character vector of column names.",
      call. = FALSE
    )
  }
  resolved <- vapply(merge, function(col) {
    if (col %in% names(aliases)) unname(aliases[[col]]) else col
  }, character(1), USE.NAMES = FALSE)
  missing <- setdiff(resolved, names_now)
  if (length(missing)) {
    stop(
      "Column(s) not found for `merge`: ",
      paste(missing, collapse = ", "),
      ". Columns were: ",
      paste(names_now, collapse = ", "),
      call. = FALSE
    )
  }
  unique(resolved)
}

#' @keywords internal
#' @noRd
.bb_df_digits <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) > 0L && all(x %% 1 == 0)) 0L else 2L
}

#' @keywords internal
#' @noRd
.bb_flextable_pairwise <- function(data, merge_cols, stars) {
  ft <- flextable::flextable(data)
  ft <- flextable::font(ft, fontname = "Times New Roman", part = "all")
  ft <- flextable::fontsize(ft, size = 12, part = "all")
  ft <- flextable::line_spacing(ft, space = 2, part = "all")

  num_cols <- names(data)[vapply(data, is.numeric, logical(1))]
  num_cols <- setdiff(num_cols, "p")
  for (col in num_cols) {
    dgt <- if (identical(col, "df")) .bb_df_digits(data[[col]]) else 2L
    ft <- flextable::colformat_double(
      ft,
      j = col,
      digits = dgt,
      big.mark = ","
    )
  }

  if ("p" %in% names(data) && is.numeric(data[["p"]])) {
    fmt_p <- function(x) rempsyc::format_p(x, stars = stars)
    ft <- flextable::set_formatter(ft, p = fmt_p)
  }

  if (length(merge_cols) && nrow(data) > 0L) {
    ft <- flextable::merge_v(ft, j = merge_cols)
  }

  ft <- flextable::theme_vanilla(ft)
  ft <- flextable::valign(ft, valign = "center", part = "all")

  italic_cols <- intersect(c("z", "t", "p", "d", "SE", "df"), names(data))
  if (length(italic_cols)) {
    ft <- flextable::italic(ft, j = italic_cols, part = "header")
  }

  flextable::set_table_properties(ft, layout = "autofit")
}
