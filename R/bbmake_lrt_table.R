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
#' An additive fixed effect, such as a blocking factor, is not a margin of
#' the interactions. The first [stats::drop1()] therefore tests it in the
#' same step as the highest-order interaction, and every later test is fit
#' without it. Name that term in `keep` to leave it in every reduced model.
#' Its own likelihood-ratio test is still reported once, when marginality
#' allows. A message names terms that stay in and additive terms that are
#' removed before the lower-order tests. Silence it with `suppressMessages()`.
#'
#' @param model A fitted model with a `drop1()` method. Typical classes:
#'   `glmmTMB`, `glm`, `lm`, `lme4::merMod`. Must include the highest-order
#'   interaction you want tested.
#' @param digits Optional number of decimal places for `LRT` and `p`. `NULL`
#'   (the default) leaves full precision so [rempsyc::nice_table()] can
#'   format the paper table. Pass `3` to match a hand-rounded pipeline.
#' @param keep Optional character vector of fixed-effect term labels to
#'   retain in every reduced model, for example `keep = "Block"` or
#'   `keep = c("Block", "Cohort")`. Labels must match the model terms
#'   (`Sex:Stim`, not `Sex * Stim`). `NULL` (the default) removes every term
#'   `drop1()` returns before the next test.
#'
#' @return A tibble with columns `Term`, `Df`, `LRT`, and `p` (highest-order
#'   terms first). Attribute `note` records that the tests are Type II LRTs
#'   and, when `keep` is set, which terms stay in every reduced model.
#'
#' @seealso [bbmake_lrt_sensitivity_table()] to stack Type II LRTs across
#'   models, [bbmake_model_table()] for Wald / *F* tables,
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
#' dat$Block <- factor(dat$id)
#' mod_block <- glm(y ~ Block + Sex * Stim, data = dat, family = poisson)
#' bbmake_lrt_table(mod_block, keep = "Block")
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
bbmake_lrt_table <- function(model, digits = NULL, keep = NULL) {
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

  keep <- .bb_lrt_resolve_keep(model, keep)
  model <- .bb_lrt_ensure_ml(model)

  labels <- .bb_fixed_term_labels(model)
  max_iter <- max(1L, length(labels))
  pieces <- vector("list", max_iter)
  n_out <- 0L
  current <- model
  # A kept term is recorded the first time drop1() can test it, then left in.
  seen_keep <- character()

  for (i in seq_len(max_iter)) {
    piece <- .bb_lrt_drop1(current)
    if (nrow(piece) == 0L) {
      break
    }
    if (i == 1L) {
      .bb_lrt_announce_keep(labels, keep, piece$Term)
    }

    piece <- piece[!piece$Term %in% seen_keep, , drop = FALSE]
    if (nrow(piece) == 0L) {
      break
    }
    n_out <- n_out + 1L
    pieces[[n_out]] <- piece
    seen_keep <- union(seen_keep, intersect(piece$Term, keep))

    drop_now <- setdiff(piece$Term, keep)
    if (length(drop_now) == 0L) {
      break
    }

    before <- .bb_fixed_term_labels(current)
    current <- .bb_drop_fixed_terms(current, drop_now)
    after <- .bb_fixed_term_labels(current)
    if (identical(sort(before), sort(after))) {
      warning(
        "Could not drop term(s) ",
        paste(drop_now, collapse = ", "),
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

  note <- "Type II likelihood-ratio tests (sequential drop1, Chisq)"
  if (length(keep) > 0L) {
    note <- paste0(
      note, "; ", paste(keep, collapse = ", "),
      " retained in every reduced model"
    )
  }
  attr(out, "note") <- note
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
.bb_lrt_check_keep <- function(keep) {
  if (is.null(keep)) {
    return(character())
  }
  if (!is.character(keep) || anyNA(keep) || any(!nzchar(keep))) {
    stop(
      "`keep` must be a character vector of fixed-effect term labels, or NULL.",
      call. = FALSE
    )
  }
  if (length(keep) == 0L) {
    return(character())
  }
  unique(keep)
}

#' @keywords internal
.bb_lrt_resolve_keep <- function(model, keep) {
  keep <- .bb_lrt_check_keep(keep)
  if (length(keep) == 0L) {
    return(character())
  }
  labels <- .bb_fixed_term_labels(model)
  missing <- setdiff(keep, labels)
  if (length(missing) > 0L) {
    available <- if (length(labels) == 0L) {
      "(none)"
    } else {
      paste(labels, collapse = ", ")
    }
    stop(
      "`keep` term(s) not in the model: ",
      paste(missing, collapse = ", "),
      ". Fixed-effect terms: ",
      available,
      ".",
      call. = FALSE
    )
  }
  keep
}

#' @keywords internal
.bb_lrt_quote_terms <- function(x) {
  paste0("`", x, "`")
}

#' @keywords internal
.bb_lrt_join_terms <- function(x) {
  x <- .bb_lrt_quote_terms(x)
  n <- length(x)
  if (n <= 1L) {
    return(x)
  }
  if (n == 2L) {
    return(paste(x, collapse = " and "))
  }
  paste0(paste(x[-n], collapse = ", "), ", and ", x[[n]])
}

#' @keywords internal
.bb_lrt_keep_phrase <- function(keep) {
  joined <- .bb_lrt_join_terms(keep)
  if (length(keep) == 1L) {
    paste0(
      joined, " stays in every reduced model. ",
      "Its likelihood-ratio test is still reported when drop1() can test it. ",
      "Later LRTs are fit with ", joined, "."
    )
  } else {
    paste0(
      joined, " stay in every reduced model. ",
      "Their likelihood-ratio tests are still reported when drop1() can test them. ",
      "Later LRTs are fit with ", joined, "."
    )
  }
}

#' @keywords internal
.bb_lrt_early_phrase <- function(early, companions, keep_already) {
  early_txt <- .bb_lrt_join_terms(early)
  verb <- if (length(early) == 1L) "is" else "are"
  pronoun <- if (length(early) == 1L) "it" else "them"
  with_txt <- if (length(companions) == 0L) {
    ""
  } else {
    paste0(" together with ", .bb_lrt_join_terms(companions))
  }
  how <- if (isTRUE(keep_already)) {
    paste0("Add ", early_txt, " to `keep`")
  } else if (length(early) == 1L) {
    paste0("Pass keep = \"", early, "\"")
  } else {
    paste0(
      "Pass keep = c(",
      paste(sprintf("\"%s\"", early), collapse = ", "),
      ")"
    )
  }
  paste0(
    early_txt, " ", verb, " removed", with_txt,
    " before lower-order terms are tested. ",
    "Later LRTs are fit without ", early_txt, ". ",
    how, " to retain ", pronoun, " in every reduced model."
  )
}

#' @keywords internal
.bb_lrt_announce_keep <- function(labels, keep, offered) {
  orders <- lengths(strsplit(offered, ":", fixed = TRUE))
  early <- if (length(offered) == 0L) {
    character()
  } else {
    offered[orders < max(orders)]
  }
  # Every current term is tested in this single drop1(), so nothing is
  # removed before a later test.
  if (length(setdiff(labels, offered)) == 0L) {
    early <- character()
  }
  early_open <- setdiff(early, keep)
  companions <- setdiff(setdiff(offered, keep), early_open)

  parts <- character()
  if (length(keep) > 0L) {
    parts <- c(parts, .bb_lrt_keep_phrase(keep))
  }
  if (length(early_open) > 0L) {
    parts <- c(
      parts,
      .bb_lrt_early_phrase(
        early_open,
        companions,
        keep_already = length(keep) > 0L
      )
    )
  }
  if (length(parts) == 0L) {
    return(invisible(NULL))
  }
  message(paste(parts, collapse = " "))
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
