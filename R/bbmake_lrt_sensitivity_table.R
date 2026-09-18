#' Sensitivity table of Type II LRTs across models
#'
#' Runs [bbmake_lrt_table()] on each fitted model in a named list and
#' stacks the results. The list names become the `Model` column, so you
#' can compare a full-data fit against outlier removals (or any other
#' refits you have already diagnosed).
#'
#' This is the **omnibus** counterpart to
#' [bbmake_pairwise_sensitivity_table()]. Pairwise tables ask whether
#' planned cell comparisons change after a refit; LRT tables ask whether
#' Type II tests of model terms (main effects and interactions) change.
#' Both are appropriate and complementary. Likelihood-ratio tests are
#' used rather than Wald / *F* ANOVA because that is the package's
#' preferred term test for GLMMs; see [bbmake_lrt_table()].
#'
#' Fit and check diagnostics yourself, then pass the models in:
#'
#' ```
#' mods <- list(
#'   Full = mod_full,
#'   `Outliers removed` = mod_no_out
#' )
#' bbmake_lrt_sensitivity_table(mods)
#' ```
#'
#' @param models A **named** list of fitted models (e.g. `glmmTMB`, `glm`,
#'   `lm`, `lme4::merMod`). Names are used as-is in the `Model` column.
#'   Each model is passed to [bbmake_lrt_table()].
#' @param digits Optional number of decimal places for `LRT` and `p`.
#'   `NULL` (the default) leaves full precision so [rempsyc::nice_table()]
#'   can format the paper table.
#'
#' @return A tibble with columns `Term`, `Model`, `Df`, `LRT`, and `p`.
#'   Rows are ordered by the Type II term order of the first model, then
#'   by model list order. Attribute `note` records that the tests are
#'   Type II LRTs stacked across models.
#'
#' @seealso [bbmake_pairwise_sensitivity_table()] for planned contrasts
#'   across the same models, [bbmake_lrt_table()] for a single-model LRT,
#'   [bbmake_model_table()] for Wald / *F* tables
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
#' bbmake_lrt_sensitivity_table(mods)
bbmake_lrt_sensitivity_table <- function(models, digits = NULL) {
  models <- .bb_as_named_model_list(models)
  if (!is.null(digits) &&
      (!is.numeric(digits) || length(digits) != 1L || is.na(digits) || digits < 0)) {
    stop("`digits` must be a single non-negative number or NULL.", call. = FALSE)
  }

  nms <- names(models)
  pieces <- vector("list", length(models))

  for (i in seq_along(models)) {
    tab <- tryCatch(
      bbmake_lrt_table(models[[i]], digits = NULL),
      error = function(e) {
        stop(
          "bbmake_lrt_table() failed for model '", nms[[i]], "': ",
          conditionMessage(e),
          call. = FALSE
        )
      }
    )
    tab$Model <- nms[[i]]
    pieces[[i]] <- tab
  }

  out <- dplyr::bind_rows(pieces)
  term_order <- unique(unlist(lapply(pieces, function(x) x$Term), use.names = FALSE))
  out[["Term"]] <- factor(out[["Term"]], levels = term_order)
  out[["Model"]] <- factor(out[["Model"]], levels = nms)

  keep <- c("Term", "Model", "Df", "LRT", "p")
  keep <- keep[keep %in% names(out)]
  out <- dplyr::select(out, dplyr::all_of(keep))
  out <- dplyr::arrange(
    out,
    dplyr::across(dplyr::all_of(c("Term", "Model")))
  )

  if (!is.null(digits) && nrow(out) > 0L) {
    out <- dplyr::mutate(
      out,
      dplyr::across(c("LRT", "p"), function(x) round(x, digits))
    )
  }

  attr(out, "note") <- paste(
    "Type II likelihood-ratio tests (sequential drop1, Chisq);",
    "stacked across models"
  )
  tibble::as_tibble(out)
}
