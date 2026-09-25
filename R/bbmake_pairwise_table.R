#' Create pairwise comparison table (preserves by= grouping + always cleans rowid)
#'
#' By default every comparison in `pw` is one family and p-values are
#' adjusted with Holm (`adjust = "none"`, `cross.adjust = "holm"`). A
#' message reports the method and how many comparisons that covers, for
#' example `Holm across 6 comparisons`. `pairs(emm)` returns one
#' `contrast` column. An interaction contrast such as
#' `contrast(emm, interaction = "pairwise")` keeps one column per factor
#' (`Stim_pairwise`, `Diet_pairwise`) plus the `by` variables.
#'
#' @param pw An `emmGrid` of contrasts: `pairs(emm)`, or an interaction
#'   contrast such as `contrast(emm, interaction = "pairwise", by = "Sex")`.
#' @param model Original fitted model (optional). Recovered from `pw` when
#'   possible; used for Gaussian Cohen's d.
#' @param adjust Within-`by`-group multiplicity adjustment passed to
#'   [emmeans::summary.emmGrid()]. Default `"none"`. `NULL` keeps the
#'   adjustment already stored on `pw` (e.g. from `pairs(..., adjust = )`).
#'   Ignored when `cross.adjust` is not `"none"`.
#' @param cross.adjust Adjustment applied to **every** pairwise test in
#'   `pw` as one family: `summary(pw, by = NULL, adjust = cross.adjust)`.
#'   That is the total-number-of-comparisons correction (e.g. 4 Sex × Diet
#'   cells, or 6 tests from 3-level Stim × 2 Groups). Default `"holm"`.
#'   Set `"none"` to skip it. This is *not* emmeans' own `cross.adjust`
#'   argument, which only adjusts matching contrasts across by-groups.
#'
#' @return A tibble with grouping columns (when present), then `contrast`
#'   or the interaction columns (`Stim_pairwise`, `Diet_pairwise`, and so
#'   on), `df` (when present), the test statistic (`z.ratio` or `t.ratio`),
#'   the effect (`IRR`, `Odds Ratio`, or `Mean Difference`), `lower.CL`,
#'   `upper.CL`, `SE`, and `p.value`. Gaussian tables that have a
#'   `contrast` column may also include Cohen's `d`. Values are left at
#'   full precision; format them with [bbnice_pairwise_table()].
#' @seealso [bbnice_pairwise_table()] to print this tibble as a flextable
#' @export
#' @examples
#' set.seed(1)
#' dat <- expand.grid(
#'   Sex = factor(c("F", "M")),
#'   Stim = factor(c("A", "B")),
#'   Diet = factor(c("C", "H")),
#'   id = 1:6
#' )
#' dat$y <- rpois(nrow(dat), lambda = 5)
#' mod <- glm(y ~ Sex * Stim * Diet, data = dat, family = poisson)
#' emmeans::emmeans(mod, ~ Stim | Diet * Sex, type = "response") |>
#'   emmeans::contrast(interaction = "pairwise", by = "Sex") |>
#'   bbmake_pairwise_table() |>
#'   bbnice_pairwise_table()
bbmake_pairwise_table <- function(
  pw,
  model = NULL,
  adjust = "none",
  cross.adjust = "holm"
) {
  if (is.null(model)) {
    model <- .bb_recover_model(pw)
  }

  by_vars <- if (!is.null(pw@misc$by.vars)) pw@misc$by.vars else character(0)

  pw_summary <- tibble::as_tibble(
    .bb_pairs_summary(pw, adjust = adjust, cross.adjust = cross.adjust)
  )

  out <- .bb_tidy_pairwise_table(pw_summary, by_vars = by_vars)

  # Interaction grids have no single contrast, so Cohen's d is left to
  # the simple pairwise case. eff_size() on a contrast-of-contrasts can
  # warn or return rows that do not line up with this table.
  if (
    !is.null(model) &&
      "Mean Difference" %in% names(out) &&
      "contrast" %in% names(out)
  ) {
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

  interaction_cols <- names(pw_summary)[
    .bb_is_interaction_header(names(pw_summary))
  ]
  interaction_cols <- setdiff(interaction_cols, by_vars)
  has_contrast <- "contrast" %in% names(pw_summary)
  if (
    !"p.value" %in% names(pw_summary) ||
      (!has_contrast && !length(interaction_cols))
  ) {
    stop(
      "pairs() did not return p.value and either contrast or a column ",
      "ending in _pairwise. Columns were: ",
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

  if (has_contrast) {
    pw_summary[["contrast"]] <- stringr::str_trim(pw_summary[["contrast"]])
  }
  for (col in interaction_cols) {
    pw_summary[[col]] <- stringr::str_trim(pw_summary[[col]])
  }

  keep <- c(
    by_vars,
    if (has_contrast) "contrast",
    interaction_cols,
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
    out <- summary(pw, infer = infer, by = NULL, adjust = cross.adjust)
    .bb_announce_adjustment(cross.adjust, nrow(out))
    return(out)
  }
  if (is.null(adjust)) {
    summary(pw, infer = infer)
  } else {
    summary(pw, infer = infer, adjust = adjust)
  }
}

#' @keywords internal
#' @noRd
.bb_announce_adjustment <- function(method, n, scope = NULL) {
  if (is.null(method) || identical(method, "none") || !is.numeric(n) || n < 1L) {
    return(invisible(NULL))
  }
  n <- as.integer(n)
  noun <- if (n == 1L) "comparison" else "comparisons"
  where <- if (is.null(scope) || !nzchar(scope)) {
    ""
  } else {
    paste0(" within each ", scope)
  }
  message(
    "P values adjusted with ", .bb_adjust_label(method),
    " across ", n, " ", noun, where, "."
  )
}

#' @keywords internal
#' @noRd
.bb_adjust_label <- function(method) {
  labels <- c(
    holm = "Holm",
    bonferroni = "Bonferroni",
    tukey = "Tukey",
    sidak = "Sidak",
    fdr = "FDR",
    hochberg = "Hochberg",
    hommel = "Hommel",
    BH = "BH",
    BY = "BY",
    scheffe = "Scheffe",
    mvt = "multivariate-t"
  )
  if (method %in% names(labels)) labels[[method]] else method
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
