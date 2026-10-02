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
#' `keep` is forwarded to [bbmake_lrt_table()] for each model, but only the
#' terms that model actually contains. A refit without the blocking factor
#' stays in the table and is fit without that term. A message names the
#' model on each per-model reminder, and names models that do not contain
#' a requested term.
#'
#' @param models A **named** list of fitted models (e.g. `glmmTMB`, `glm`,
#'   `lm`, `lme4::merMod`). Names are used as-is in the `Model` column.
#'   Each model is passed to [bbmake_lrt_table()].
#' @param digits Optional number of decimal places for `LRT` and `p`.
#'   `NULL` (the default) leaves full precision so [rempsyc::nice_table()]
#'   can format the paper table.
#' @param keep Optional character vector of fixed-effect term labels, as in
#'   [bbmake_lrt_table()], for example `keep = "Block"`. Each model retains
#'   the requested terms it contains. `NULL` (the default) removes every
#'   term `drop1()` returns before the next test.
#'
#' @return A tibble with columns `Term`, `Model`, `Df`, `LRT`, and `p`.
#'   Rows are ordered by the Type II term order of the first model, then
#'   by model list order. Models that lack a term have no row for it.
#'   Attribute `note` records that the tests are Type II LRTs stacked
#'   across models and, when `keep` is set, that those terms are retained
#'   in models that contain them.
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
#'
#' dat$Block <- factor(dat$id)
#' mods_block <- list(
#'   `With block` = glm(y ~ Block + Sex * Stim, data = dat, family = poisson),
#'   `No block` = glm(y ~ Sex * Stim, data = dat, family = poisson)
#' )
#' bbmake_lrt_sensitivity_table(mods_block, keep = "Block")
bbmake_lrt_sensitivity_table <- function(models, digits = NULL, keep = NULL) {
  models <- .bb_as_named_model_list(models)
  if (!is.null(digits) &&
      (!is.numeric(digits) || length(digits) != 1L || is.na(digits) || digits < 0)) {
    stop("`digits` must be a single non-negative number or NULL.", call. = FALSE)
  }
  keep <- .bb_lrt_check_keep(keep)
  if (length(keep) > 0L) {
    .bb_lrt_require_keep_somewhere(models, keep)
  }

  nms <- names(models)
  pieces <- vector("list", length(models))
  model_messages <- rep(list(character()), length(models))
  kept_by_model <- rep(list(character()), length(models))

  for (i in seq_along(models)) {
    keep_i <- .bb_lrt_keep_present(models[[i]], keep)
    kept_by_model[[i]] <- keep_i
    caught <- tryCatch(
      .bb_capture_lrt_table(models[[i]], keep_i),
      error = function(e) {
        stop(
          "bbmake_lrt_table() failed for model '", nms[[i]], "': ",
          conditionMessage(e),
          call. = FALSE
        )
      }
    )
    tab <- caught$tab
    tab$Model <- nms[[i]]
    pieces[[i]] <- tab
    model_messages[[i]] <- caught$messages
  }

  .bb_lrt_sensitivity_messages(nms, model_messages, keep, kept_by_model)

  out <- dplyr::bind_rows(pieces)
  term_order <- unique(unlist(lapply(pieces, function(x) x$Term), use.names = FALSE))
  out[["Term"]] <- factor(out[["Term"]], levels = term_order)
  out[["Model"]] <- factor(out[["Model"]], levels = nms)

  cols <- c("Term", "Model", "Df", "LRT", "p")
  cols <- cols[cols %in% names(out)]
  out <- dplyr::select(out, dplyr::all_of(cols))
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

  note <- paste(
    "Type II likelihood-ratio tests (sequential drop1, Chisq);",
    "stacked across models"
  )
  if (length(keep) > 0L) {
    note <- paste0(
      note, "; ", paste(keep, collapse = ", "),
      " retained in models that contain ",
      if (length(keep) == 1L) "it" else "them"
    )
  }
  attr(out, "note") <- note
  tibble::as_tibble(out)
}

#' @keywords internal
.bb_lrt_keep_present <- function(model, keep) {
  if (length(keep) == 0L) {
    return(character())
  }
  intersect(keep, .bb_fixed_term_labels(model))
}

#' @keywords internal
.bb_lrt_require_keep_somewhere <- function(models, keep) {
  present <- unique(unlist(
    lapply(models, function(model) .bb_lrt_keep_present(model, keep)),
    use.names = FALSE
  ))
  missing <- setdiff(keep, present)
  if (length(missing) == 0L) {
    return(invisible(NULL))
  }
  stop(
    "`keep` term(s) not in any model: ",
    paste(missing, collapse = ", "),
    ".",
    call. = FALSE
  )
}

#' @keywords internal
.bb_capture_lrt_table <- function(model, keep) {
  keep <- if (length(keep) == 0L) NULL else keep
  msgs <- character()
  tab <- withCallingHandlers(
    bbmake_lrt_table(model, digits = NULL, keep = keep),
    message = function(m) {
      msgs <<- c(msgs, conditionMessage(m))
      invokeRestart("muffleMessage")
    }
  )
  list(tab = tab, messages = msgs)
}

#' @keywords internal
.bb_lrt_sensitivity_messages <- function(model_names, model_messages, keep, kept_by_model) {
  lines <- character()
  for (i in seq_along(model_names)) {
    text <- trimws(paste(model_messages[[i]], collapse = " "))
    if (!nzchar(text)) {
      next
    }
    lines <- c(
      lines,
      paste0(.bb_lrt_quote_terms(model_names[[i]]), ": ", text)
    )
  }
  if (length(keep) > 0L) {
    for (term in keep) {
      lacks <- model_names[vapply(
        kept_by_model,
        function(x) !(term %in% x),
        logical(1)
      )]
      if (length(lacks) == 0L) {
        next
      }
      verb <- if (length(lacks) == 1L) "does not contain" else "do not contain"
      fitted <- if (length(lacks) == 1L) "that model is" else "those models are"
      lines <- c(lines, paste0(
        .bb_lrt_join_terms(lacks), " ", verb, " ",
        .bb_lrt_quote_terms(term), ", so ", fitted, " fit without it."
      ))
    }
  }
  if (length(lines) == 0L) {
    return(invisible(NULL))
  }
  message(paste(lines, collapse = "\n"))
}
