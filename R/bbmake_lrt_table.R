#' Type II likelihood-ratio table
#'
#' Sequential [stats::drop1()] chi-squared tests from the highest-order
#' interaction down to main effects. This is the LRT analogue of Type II
#' ANOVA: each term is tested given all other terms of equal or lower order,
#' respecting marginality.
#'
#' For `glmmTMB` (and other non-Gaussian GLMMs) this is usually preferable to
#' [car::Anova()] Wald tests when the sample is small. Wald *p*-values rely on
#' a quadratic approximation to the log-likelihood that is often poor with
#' few observations, beta/nbinom families, or few random-effect levels. The
#' LRT compares actual maximized likelihoods.
#'
#' Pass the **full** model (the one whose name you would have suffixed with
#' the highest interaction order, e.g. `mod_fc_pref2`). Reduced models are
#' fitted internally; you do not need to `update()` them yourself.
#'
#' Random-effect terms, offsets, and `glmmTMB` zero-inflation / dispersion
#' formulas are left untouched. Only fixed effects of the conditional model
#' are tested.
#'
#' @param model A fitted model with a `drop1()` method. Typical classes:
#'   `glmmTMB`, `glm`, `lm`, `lme4::merMod`. Must include the highest-order
#'   interaction you want tested.
#' @param digits Optional number of decimal places for `LRT` and `p`. `NULL`
#'   (the default) leaves full precision so [rempsyc::nice_table()] can
#'   format the paper table. Pass `3` to match a hand-rounded pipeline.
#'
#' @return A tibble with columns `Term`, `Df`, `LRT`, and `p` (highest-order
#'   terms first). Attribute `note` records that the tests are Type II LRTs.
#'
#' @seealso [bbmake_model_table()] for Wald / *F* tables,
#'   [stats::drop1()]
#' @export
#' @examples
#' set.seed(1)
#' dat <- expand.grid(
#'   Sex = factor(c("F", "M")),
#'   Stim = factor(c("A", "B")),
#'   id = 1:10
#' )
#' dat$y <- rpois(nrow(dat), lambda = 4)
#' mod <- glm(y ~ Sex * Stim, data = dat, family = poisson)
#' bbmake_lrt_table(mod)
#'
#' \dontrun{
#' mod <- glmmTMB::glmmTMB(
#'   pref ~ Sex * Stim,
#'   data = dat,
#'   family = glmmTMB::beta_family()
#' )
#' lrt <- bbmake_lrt_table(mod)
#' rempsyc::nice_table(lrt)
#' }
bbmake_lrt_table <- function(model, digits = NULL) {
  if (missing(model)) {
    stop(
      "`model` is missing. Pass the full fitted model ",
      "(highest-order interaction included).",
      call. = FALSE
    )
  }
  if (!is.null(digits) &&
      (!is.numeric(digits) || length(digits) != 1L || is.na(digits) || digits < 0)) {
    stop("`digits` must be a single non-negative number or NULL.", call. = FALSE)
  }

  model <- .bb_lrt_ensure_ml(model)

  labels <- .bb_fixed_term_labels(model)
  max_iter <- max(1L, length(labels))
  pieces <- vector("list", max_iter)
  n_out <- 0L
  current <- model

  for (i in seq_len(max_iter)) {
    piece <- .bb_lrt_drop1(current)
    if (nrow(piece) == 0L) {
      break
    }
    n_out <- n_out + 1L
    pieces[[n_out]] <- piece

    before <- .bb_fixed_term_labels(current)
    current <- .bb_drop_fixed_terms(current, piece$Term)
    after <- .bb_fixed_term_labels(current)
    if (identical(sort(before), sort(after))) {
      warning(
        "Could not drop term(s) ",
        paste(piece$Term, collapse = ", "),
        "; stopping the LRT sequence.",
        call. = FALSE
      )
      break
    }
  }

  if (n_out == 0L) {
    out <- .bb_empty_lrt()
  } else {
    out <- dplyr::bind_rows(pieces[seq_len(n_out)])
  }

  if (any(out$Df == 0, na.rm = TRUE)) {
    warning(
      "Some terms have 0 df (rank deficiency or unidentifiable). ",
      "Interpret those LRTs with caution.",
      call. = FALSE
    )
  }

  if (!is.null(digits) && nrow(out) > 0L) {
    out <- dplyr::mutate(
      out,
      dplyr::across(c("LRT", "p"), function(x) round(x, digits))
    )
  }

  attr(out, "note") <- "Type II likelihood-ratio tests (sequential drop1, Chisq)"
  tibble::as_tibble(out)
}

# ── Internal helpers ──────────────────────────────────────────────────────────

#' @keywords internal
.bb_empty_lrt <- function() {
  tibble::tibble(
    Term = character(),
    Df   = integer(),
    LRT  = numeric(),
    p    = numeric()
  )
}

#' @keywords internal
.bb_fixed_term_labels <- function(model) {
  tl <- attr(stats::terms(model), "term.labels")
  if (is.null(tl)) character() else tl
}

#' @keywords internal
.bb_lrt_ensure_ml <- function(model) {
  reml <- tryCatch(
    inherits(model, "glmmTMB") && isTRUE(model$modelInfo$REML),
    error = function(e) FALSE
  )
  if (!isTRUE(reml)) {
    return(model)
  }
  warning(
    "Model was fitted with REML. Refitting with ML for likelihood-ratio ",
    "tests of fixed effects.",
    call. = FALSE
  )
  .bb_update_model(model, REML = FALSE)
}

#' @keywords internal
.bb_lrt_drop1 <- function(model) {
  d1 <- tryCatch(
    stats::drop1(model, test = "Chisq"),
    error = function(e) {
      stop(
        "drop1() failed for ",
        paste(deparse(stats::formula(model)), collapse = " "),
        ": ",
        conditionMessage(e),
        call. = FALSE
      )
    }
  )

  out <- tibble::rownames_to_column(as.data.frame(d1), "Term")
  out <- out[!out$Term %in% c("<none>", "1", "(Intercept)"), , drop = FALSE]
  if (nrow(out) == 0L) {
    return(.bb_empty_lrt())
  }

  df_col  <- .bb_first_col(out, c("Df", "df"))
  lrt_col <- .bb_first_col(out, c("LRT", "Chisq", "scaled dev."))
  p_col   <- .bb_first_col(out, c("Pr(>Chi)", "Pr(>Chisq)"))

  if (is.na(df_col) || is.na(p_col)) {
    stop(
      "drop1() did not return Df and a chi-squared p-value. Columns were: ",
      paste(names(out), collapse = ", "),
      call. = FALSE
    )
  }

  lrt <- if (is.na(lrt_col)) {
    vapply(out$Term, function(term) {
      reduced <- .bb_drop_fixed_terms(model, term)
      as.numeric(2 * (stats::logLik(model) - stats::logLik(reduced)))
    }, numeric(1))
  } else {
    out[[lrt_col]]
  }

  tibble::tibble(
    Term = out$Term,
    Df   = as.integer(out[[df_col]]),
    LRT  = as.numeric(lrt),
    p    = as.numeric(out[[p_col]])
  )
}

#' @keywords internal
.bb_first_col <- function(df, candidates) {
  hit <- candidates[candidates %in% names(df)]
  if (length(hit) == 0L) NA_character_ else hit[[1L]]
}

#' @keywords internal
.bb_drop_fixed_terms <- function(model, terms) {
  terms <- terms[!is.na(terms) & nzchar(terms)]
  if (length(terms) == 0L) {
    return(model)
  }
  rhs <- paste(paste0("-", terms), collapse = " ")
  f <- stats::as.formula(
    paste("~ .", rhs),
    env = environment(stats::formula(model))
  )
  tryCatch(
    .bb_update_model(model, f),
    error = function(e) {
      stop(
        "Failed to drop term(s) ",
        paste(terms, collapse = ", "),
        ": ",
        conditionMessage(e),
        call. = FALSE
      )
    }
  )
}

#' @keywords internal
.bb_update_model <- function(model, ...) {
  cl <- stats::update(model, ..., evaluate = FALSE)
  env <- environment(stats::formula(model))
  if (is.null(env) || identical(env, emptyenv())) {
    env <- attr(stats::terms(model), ".Environment")
  }
  tryCatch(
    eval(cl, envir = env),
    error = function(e) {
      mf <- tryCatch(stats::model.frame(model), error = function(e2) NULL)
      if (is.null(mf)) {
        stop(e)
      }
      stats::update(model, ..., data = mf, subset = NULL)
    }
  )
}
