#' Sensitivity table of planned contrasts across models
#'
#' Runs the same [emmeans::emmeans()] pairwise contrast on each fitted
#' model in a named list and stacks the results. The list names become
#' the `Model` column, so you can compare a full-data fit against outlier
#' removals (or any other refits you have already diagnosed).
#'
#' This is the **pairwise** counterpart to [bbmake_lrt_sensitivity_table()].
#' Pairwise tables ask whether planned cell comparisons change after a
#' refit; LRT tables ask whether Type II tests of model terms change.
#'
#' Fit and check diagnostics yourself, then pass the models in:
#'
#' ```
#' mods <- list(
#'   Full = mod_full,
#'   `Outliers removed` = mod_no_out
#' )
#' bbmake_pairwise_sensitivity_table(
#'   mods, ~ Stim | Diet * Sex, by = c("Sex", "Diet")
#' )
#' ```
#'
#' On the response scale, pairwise ratios are labelled `IRR` (counts,
#' rates) or `Odds Ratio` (binomial). Gaussian / link-scale pairs keep
#' `estimate`. Random effects are handled by emmeans in the usual way.
#'
#' @param models A **named** list of fitted models (e.g. `glmmTMB`, `glm`,
#'   `lm`, `lme4::merMod`). Names are used as-is in the `Model` column.
#' @param specs An emmeans specs formula, passed to [emmeans::emmeans()]
#'   (e.g. `~ Stim | Diet * Sex`).
#' @param by Optional grouping for pairwise contrasts. When `NULL` (the
#'   default), the by-variables already encoded in `specs` (after `|`) are
#'   used. Supply an explicit order such as `c("Sex", "Diet")` to control
#'   both pairing and column order.
#' @param reverse If `TRUE` (the default), reverse the contrast direction
#'   so later factor levels are in the numerator (the usual IRR direction).
#' @param type Passed to [emmeans::emmeans()]. Default `"response"` so count
#'   models return IRRs rather than log-scale differences.
#' @param adjust Multiplicity adjustment passed to [emmeans::contrast()]
#'   (the same argument as `pairs(..., adjust = )`). Default `"tukey"`,
#'   matching emmeans pairwise. Common values: `"tukey"`, `"bonferroni"`,
#'   `"holm"`, `"fdr"`, `"none"`. Adjustment is **within** each `by` group;
#'   with one contrast per cell (two-level `Stim`) Tukey and none coincide.
#' @param cross.adjust Additional p-value adjustment **across** `by` groups,
#'   passed to [emmeans::summary.emmGrid()]. Default `"none"`. Use this when
#'   each cell has one planned contrast (e.g. Stim within Sex × Diet) and
#'   you want those cells treated as one family:
#'   `adjust = "none", cross.adjust = "bonferroni"`. Valid methods are
#'   [stats::p.adjust.methods] plus `"sidak"`. Ignored unless there is more
#'   than one `by` group and the groups are the same size (emmeans rules).
#' @param digits Optional decimal places for the effect column and `p`.
#'   `NULL` (the default) leaves full precision so [rempsyc::nice_table()]
#'   can format the paper table.
#'
#' @return A tibble with grouping columns from `by` (when present), then
#'   `Contrast`, `Model`, the effect (`IRR`, `Odds Ratio`, or `estimate`),
#'   and `p`. Rows are ordered by grouping variables, contrast, then model
#'   list order.
#'
#' @seealso [bbmake_lrt_sensitivity_table()] for Type II LRTs across the
#'   same models, [bbmake_pairwise_table()] for a single-model pairwise
#'   table, [bbmake_lrt_table()] for a single-model LRT
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
#' mods <- list(
#'   Full = glm(y ~ Sex * Stim * Diet, data = dat, family = poisson),
#'   Reduced = glm(
#'     y ~ Sex * Stim * Diet,
#'     data = dat[dat$id != "1", ],
#'     family = poisson
#'   )
#' )
#' bbmake_pairwise_sensitivity_table(
#'   mods, ~ Stim | Diet * Sex, by = c("Sex", "Diet")
#' )
#' bbmake_pairwise_sensitivity_table(
#'   mods, ~ Stim | Diet * Sex,
#'   by = c("Sex", "Diet"),
#'   adjust = "none",
#'   cross.adjust = "bonferroni"
#' )
bbmake_pairwise_sensitivity_table <- function(
  models,
  specs,
  by = NULL,
  reverse = TRUE,
  type = "response",
  adjust = "tukey",
  cross.adjust = "none",
  digits = NULL
) {
  models <- .bb_as_named_model_list(models)
  if (missing(specs) || is.null(specs)) {
    stop(
      "`specs` is missing. Pass an emmeans formula such as `~ Stim | Diet * Sex`.",
      call. = FALSE
    )
  }
  if (!is.null(digits) &&
      (!is.numeric(digits) || length(digits) != 1L || is.na(digits) || digits < 0)) {
    stop("`digits` must be a single non-negative number or NULL.", call. = FALSE)
  }

  nms <- names(models)
  pieces <- vector("list", length(models))
  by_vars <- by

  for (i in seq_along(models)) {
    piece <- .bb_pairs_for_model(
      model = models[[i]],
      model_name = nms[[i]],
      specs = specs,
      by = by,
      reverse = reverse,
      type = type,
      adjust = adjust,
      cross.adjust = cross.adjust
    )
    if (i == 1L && (is.null(by_vars) || length(by_vars) == 0L)) {
      by_vars <- piece$by_vars
    }
    pieces[[i]] <- piece$data
  }

  out <- dplyr::bind_rows(pieces)
  out <- .bb_tidy_pairwise_sensitivity(out, by_vars = by_vars, model_levels = nms)

  if (!is.null(digits) && nrow(out) > 0L) {
    num_cols <- intersect(
      c("IRR", "Odds Ratio", "estimate", "p"),
      names(out)
    )
    if (length(num_cols) > 0L) {
      out <- dplyr::mutate(
        out,
        dplyr::across(dplyr::all_of(num_cols), function(x) round(x, digits))
      )
    }
  }

  tibble::as_tibble(out)
}

# ── Internal helpers ──────────────────────────────────────────────────────────

#' @keywords internal
.bb_as_named_model_list <- function(models) {
  if (missing(models) || is.null(models)) {
    stop(
      "`models` is missing. Pass a named list of fitted models, e.g. ",
      "list(Full = mod, `Outliers removed` = mod_no_out).",
      call. = FALSE
    )
  }
  if (inherits(models, c(
    "glmmTMB", "lm", "glm", "merMod", "lmerModLmerTest", "glmerMod"
  ))) {
    stop(
      "Pass a named list of models, not a single model. For example:\n",
      "list(Full = mod, `Outliers removed` = mod_no_out)",
      call. = FALSE
    )
  }
  if (!is.list(models) || length(models) < 1L) {
    stop("`models` must be a named list of fitted models.", call. = FALSE)
  }
  nms <- names(models)
  if (is.null(nms) || any(!nzchar(nms) | is.na(nms))) {
    stop(
      "`models` must be a named list; names become the Model column.",
      call. = FALSE
    )
  }
  models
}

#' @keywords internal
.bb_pairs_for_model <- function(
  model,
  model_name,
  specs,
  by,
  reverse,
  type,
  adjust,
  cross.adjust
) {
  emm <- tryCatch(
    emmeans::emmeans(model, specs = specs, type = type),
    error = function(e) {
      stop(
        "emmeans() failed for model '", model_name, "': ",
        conditionMessage(e),
        call. = FALSE
      )
    }
  )

  by_use <- by
  if (is.null(by_use) || length(by_use) == 0L) {
    by_use <- emm@misc$by.vars
  }

  method <- if (isTRUE(reverse)) "revpairwise" else "pairwise"
  pw <- tryCatch(
    if (is.null(by_use) || length(by_use) == 0L) {
      emmeans::contrast(emm, method = method, adjust = adjust)
    } else {
      emmeans::contrast(emm, method = method, by = by_use, adjust = adjust)
    },
    error = function(e) {
      stop(
        "contrast() failed for model '", model_name, "': ",
        conditionMessage(e),
        call. = FALSE
      )
    }
  )

  df <- as.data.frame(pw, adjust = adjust, cross.adjust = cross.adjust)
  df$Model <- model_name
  list(
    data = df,
    by_vars = if (is.null(by_use)) character() else as.character(by_use)
  )
}

#' @keywords internal
.bb_tidy_pairwise_sensitivity <- function(out, by_vars, model_levels) {
  if ("ratio" %in% names(out)) {
    out <- dplyr::rename(out, IRR = "ratio")
  } else if ("odds.ratio" %in% names(out)) {
    out <- dplyr::rename(out, `Odds Ratio` = "odds.ratio")
  }

  if (!"contrast" %in% names(out) || !"p.value" %in% names(out)) {
    stop(
      "contrast() did not return contrast and p.value. Columns were: ",
      paste(names(out), collapse = ", "),
      call. = FALSE
    )
  }

  effect_cols <- intersect(c("IRR", "Odds Ratio", "estimate"), names(out))
  if (length(effect_cols) == 0L) {
    stop(
      "contrast() did not return ratio, odds.ratio, or estimate. Columns were: ",
      paste(names(out), collapse = ", "),
      call. = FALSE
    )
  }

  out[["contrast"]] <- stringr::str_trim(out[["contrast"]])
  out[["Model"]] <- factor(out[["Model"]], levels = model_levels)
  out <- dplyr::rename(out, Contrast = "contrast", p = "p.value")

  keep <- c(by_vars, "Contrast", "Model", effect_cols, "p")
  keep <- keep[keep %in% names(out)]
  out <- dplyr::select(out, dplyr::all_of(keep))

  dplyr::arrange(
    out,
    dplyr::across(dplyr::any_of(c(by_vars, "Contrast", "Model")))
  )
}
