#' Sensitivity table of planned contrasts across models
#'
#' Runs the same [emmeans::emmeans()] pairwise contrast on each fitted
#' model in a named list, tidies each with [bbmake_pairwise_table()], and
#' stacks the results. The list names become the `Model` column, so you
#' can compare a full-data fit against outlier removals (or any other
#' refits you have already diagnosed).
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
#' rates) or `Odds Ratio` (binomial). Gaussian / link-scale pairs use
#' `Mean Difference`. Random effects are handled by emmeans in the usual
#' way. The default matches [bbmake_pairwise_table()]: Holm across every
#' comparison **within each model**, not across the stacked rows. A
#' message reports that count, for example
#' `Holm across 4 comparisons within each model`.
#'
#' `interaction = "pairwise"` (or `"revpairwise"`) contrasts the factors
#' with each other instead of running a simple pairwise test. Factors
#' named in `by` are held out of that interaction. `reverse` does not
#' change an `interaction` value you supply.
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
#'   Ignored when `interaction` is set; pass `interaction = "revpairwise"`
#'   for that direction.
#' @param type Passed to [emmeans::emmeans()]. Default `"response"` so count
#'   models return IRRs rather than log-scale differences.
#' @param adjust Within-`by`-group multiplicity adjustment passed to
#'   [emmeans::contrast()] (the same argument as `pairs(..., adjust = )`).
#'   Default `"none"`. Common values: `"tukey"`, `"bonferroni"`, `"holm"`,
#'   `"fdr"`, `"none"`. Adjustment is **within** each `by` group; with one
#'   contrast per cell (two-level `Stim`) Tukey and none coincide. Ignored
#'   when `cross.adjust` is not `"none"`.
#' @param cross.adjust Adjustment for **every** pairwise test in one model
#'   as a single family (`summary(pw, by = NULL, adjust = cross.adjust)`):
#'   all Stim comparisons in all Sex × Diet cells together. Default
#'   `"holm"`. Set `"none"` to skip it, or
#'   `adjust = "tukey", cross.adjust = "none"` for Tukey within each
#'   by-group. Valid methods are [stats::p.adjust.methods] plus `"sidak"`.
#'   Applied separately to each model, not to the stacked rows.
#' @param interaction `NULL` (the default) runs a simple pairwise contrast.
#'   `"pairwise"` or `"revpairwise"` is passed to [emmeans::contrast()]
#'   as `interaction`. Factors named in `by` are held out; the remaining
#'   factors are contrasted with each other. On `~ Stim | Diet * Sex`,
#'   `by = "Sex"` and `interaction = "pairwise"` contrasts Stim and Diet
#'   within Sex. `by = NULL` still uses every variable after `|`, so that
#'   formula would contrast only Stim and keep Diet as a grouping column.
#'
#' @return A tibble with grouping columns from `by` (when present), then
#'   `contrast` or the interaction columns (`Stim_pairwise`,
#'   `Diet_pairwise`, and so on), `Model`, the test statistic (`z.ratio`
#'   or `t.ratio`), `df` (when present), the effect (`IRR`, `Odds Ratio`,
#'   or `Mean Difference`), `SE`, `lower.CL`, `upper.CL`, and `p.value`.
#'   Rows are ordered by grouping variables, the contrast columns, then
#'   model list order. Values are left at full precision; format them with
#'   [bbnice_pairwise_table()].
#'
#' @seealso [bbnice_pairwise_table()] to print this tibble as a flextable,
#'   [bbmake_lrt_sensitivity_table()] for Type II LRTs across the
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
#'   cross.adjust = "bonferroni"
#' )
#' bbmake_pairwise_sensitivity_table(
#'   mods,
#'   ~ Stim | Diet * Sex,
#'   by = "Sex",
#'   interaction = "pairwise"
#' ) |>
#'   bbnice_pairwise_table()
bbmake_pairwise_sensitivity_table <- function(
  models,
  specs,
  by = NULL,
  reverse = TRUE,
  type = "response",
  adjust = "none",
  cross.adjust = "holm",
  interaction = NULL
) {
  models <- .bb_as_named_model_list(models)
  if (missing(specs) || is.null(specs)) {
    stop(
      "`specs` is missing. Pass an emmeans formula such as `~ Stim | Diet * Sex`.",
      call. = FALSE
    )
  }
  interaction <- .bb_as_interaction(interaction)

  nms <- names(models)
  pieces <- vector("list", length(models))
  by_vars <- by

  for (i in seq_along(models)) {
    # One announcement for the whole table, not one copy per model.
    piece <- suppressMessages(.bb_pairs_for_model(
      model = models[[i]],
      model_name = nms[[i]],
      specs = specs,
      by = by,
      reverse = reverse,
      type = type,
      adjust = adjust,
      cross.adjust = cross.adjust,
      interaction = interaction
    ))
    if (i == 1L && (is.null(by_vars) || length(by_vars) == 0L)) {
      by_vars <- piece$by_vars
    }
    pieces[[i]] <- piece$data
  }

  .bb_announce_sensitivity_adjustment(
    cross.adjust,
    vapply(pieces, nrow, integer(1)),
    nms
  )

  out <- dplyr::bind_rows(pieces)
  out[["Model"]] <- factor(out[["Model"]], levels = nms)

  interaction_cols <- names(out)[.bb_is_interaction_header(names(out))]
  interaction_cols <- setdiff(interaction_cols, c(by_vars, "Model"))

  keep <- c(
    by_vars,
    "contrast",
    interaction_cols,
    "Model",
    "z.ratio",
    "t.ratio",
    "df",
    "IRR",
    "Odds Ratio",
    "Mean Difference",
    "SE",
    "lower.CL",
    "upper.CL",
    "p.value"
  )
  keep <- keep[keep %in% names(out)]
  out <- dplyr::select(out, dplyr::all_of(keep))

  dplyr::arrange(
    tibble::as_tibble(out),
    dplyr::across(dplyr::any_of(c(
      by_vars, "contrast", interaction_cols, "Model"
    )))
  )
}

# ── Internal helpers ──────────────────────────────────────────────────────────

#' @keywords internal
#' @noRd
.bb_as_interaction <- function(interaction) {
  if (is.null(interaction)) {
    return(NULL)
  }
  ok <- is.character(interaction) &&
    length(interaction) == 1L &&
    !is.na(interaction) &&
    interaction %in% c("pairwise", "revpairwise")
  if (!ok) {
    stop(
      "`interaction` must be NULL, \"pairwise\", or \"revpairwise\".",
      call. = FALSE
    )
  }
  interaction
}

#' @keywords internal
#' @noRd
.bb_announce_sensitivity_adjustment <- function(method, counts, model_names) {
  if (is.null(method) || identical(method, "none") || length(counts) == 0L) {
    return(invisible(NULL))
  }
  if (length(unique(counts)) == 1L) {
    return(.bb_announce_adjustment(method, counts[[1]], scope = "model"))
  }
  detail <- paste(sprintf("%s: %d", model_names, counts), collapse = "; ")
  message(
    "P values adjusted with ", .bb_adjust_label(method),
    " within each model (", detail, " comparisons)."
  )
}

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
  cross.adjust,
  interaction = NULL
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

  pw <- tryCatch(
    {
      if (!is.null(interaction)) {
        if (is.null(by_use) || length(by_use) == 0L) {
          emmeans::contrast(emm, interaction = interaction, adjust = adjust)
        } else {
          emmeans::contrast(
            emm,
            interaction = interaction,
            by = by_use,
            adjust = adjust
          )
        }
      } else {
        method <- if (isTRUE(reverse)) "revpairwise" else "pairwise"
        if (is.null(by_use) || length(by_use) == 0L) {
          emmeans::contrast(emm, method = method, adjust = adjust)
        } else {
          emmeans::contrast(emm, method = method, by = by_use, adjust = adjust)
        }
      }
    },
    error = function(e) {
      stop(
        "contrast() failed for model '", model_name, "': ",
        conditionMessage(e),
        call. = FALSE
      )
    }
  )

  tab <- bbmake_pairwise_table(
    pw,
    model = NULL,
    adjust = adjust,
    cross.adjust = cross.adjust
  )
  tab$Model <- model_name
  list(
    data = tab,
    by_vars = if (is.null(by_use)) character() else as.character(by_use)
  )
}
