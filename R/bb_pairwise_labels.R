#' Build pairwise significance labels for plotting
#'
#' Creates a data frame suitable for [ggpubr::stat_pvalue_manual()] from
#' pairwise contrasts: `group1`, `group2`, `p.signif`, `y.position`, plus any
#' by-variables needed for faceting. Also attaches an `effect_annotation`
#' column for optional corner labels.
#'
#' @param pw An `emmGrid` of pairwise contrasts (e.g. `pairs(emm)`), or `NULL`
#'   if `pw_table` is supplied.
#' @param emm An `emmGrid` (or data frame / result of [bb_emm_df()]) used to
#'   place significance brackets above the tallest CI in each panel.
#' @param pw_table Optional output of [bbmake_pairwise_table()]. Built from
#'   `pw` when omitted.
#' @param model Optional fitted model forwarded to [bbmake_pairwise_table()]
#'   for effect sizes (e.g. Cohen's d).
#' @param y.adjust Vertical nudge added to each bracket's base y-position.
#'   A scalar is recycled; a vector is applied after stacking order within
#'   each panel (length must match the number of contrasts or the number of
#'   panels when there is one contrast per panel).
#' @param step Fraction of the overall y-range used to stack multiple
#'   brackets within the same panel. Default `0.08`. Set to `0` to keep all
#'   brackets at the same height (plus `y.adjust`).
#' @param hide.ns If `TRUE`, drop rows with `p.signif == "ns"`.
#'
#' @return A tibble with columns `group1`, `group2`, `p.value`, `p.signif`,
#'   `y.position`, `effect_annotation`, plus by-variables when present.
#'   Attribute `by.vars` records faceting columns.
#' @export
#' @examples
#' \dontrun{
#' emm <- emmeans(mod, ~ Stim | Group)
#' pw  <- pairs(emm)
#' sig <- bb_pairwise_labels(pw, emm = emm, model = mod)
#'
#' df <- bb_emm_df(emm)
#' ggplot(df, aes(x = Stim, y = y, ymin = ymin, ymax = ymax)) +
#'   geom_pointrange() +
#'   facet_wrap(~ Group) +
#'   ggpubr::stat_pvalue_manual(
#'     sig, label = "p.signif",
#'     xmin = "group1", xmax = "group2",
#'     y.position = "y.position"
#'   )
#' }
bb_pairwise_labels <- function(
    pw = NULL,
    emm = NULL,
    pw_table = NULL,
    model = NULL,
    y.adjust = 0,
    step = 0.08,
    hide.ns = FALSE
) {
  if (is.null(pw) && is.null(pw_table)) {
    stop("Provide at least one of `pw` or `pw_table`.", call. = FALSE)
  }

  # ── Resolve pw_table ────────────────────────────────────────────────────
  if (is.null(pw_table)) {
    if (!inherits(pw, "emmGrid")) {
      stop("`pw` must be an emmGrid from pairs() when `pw_table` is NULL.",
           call. = FALSE)
    }
    pw_table <- bbmake_pairwise_table(pw, model = model)
  }
  pw_table <- tibble::as_tibble(pw_table)

  if (!"contrast" %in% names(pw_table) || !"p.value" %in% names(pw_table)) {
    stop("`pw_table` must contain `contrast` and `p.value` columns.",
         call. = FALSE)
  }

  # ── By-variables ────────────────────────────────────────────────────────
  by_vars <- character(0)
  if (!is.null(pw) && inherits(pw, "emmGrid")) {
    by_vars <- pw@misc$by.vars %||% character(0)
  }
  if (length(by_vars) == 0L && !is.null(emm) && inherits(emm, "emmGrid")) {
    by_vars <- emm@misc$by.vars %||% character(0)
  }
  if (length(by_vars) == 0L) {
    # Fall back to columns present in both table and emm_df (non-contrast meta)
    meta <- c(
      "contrast", "p.value", "t.ratio", "z.ratio", "df", "SE",
      "estimate", "ratio", "odds.ratio", "Mean Difference",
      "IRR", "Odds Ratio", "d", "Cohen's d", "lower.CL", "upper.CL",
      "rowid", "% Change"
    )
    candidates <- setdiff(names(pw_table), meta)
    if (!is.null(emm)) {
      emm_df_tmp <- if (is.data.frame(emm) && !inherits(emm, "emmGrid")) {
        emm
      } else {
        bb_emm_df(emm)
      }
      candidates <- intersect(candidates, names(emm_df_tmp))
    }
    by_vars <- candidates
  }
  by_vars <- by_vars[by_vars %in% names(pw_table)]

  # ── Significance stars + group split ────────────────────────────────────
  labels <- pw_table |>
    dplyr::mutate(
      p.signif = dplyr::case_when(
        .data$p.value < 0.001 ~ "***",
        .data$p.value < 0.01  ~ "**",
        .data$p.value < 0.05  ~ "*",
        TRUE                  ~ "ns"
      )
    )

  labels <- .bb_split_contrast_groups(labels)
  labels <- .bb_add_effect_annotation(labels)

  if (isTRUE(hide.ns)) {
    labels <- dplyr::filter(labels, .data$p.signif != "ns")
  }

  # ── y.position from emm ─────────────────────────────────────────────────
  if (is.null(emm)) {
    labels <- dplyr::mutate(labels, y.position = NA_real_)
  } else {
    emm_df <- if (is.data.frame(emm) && "ymax" %in% names(emm)) {
      tibble::as_tibble(emm)
    } else {
      bb_emm_df(emm)
    }

    y_base <- .bb_ypos_lookup(emm_df, by_vars)

    if (length(by_vars) > 0L) {
      labels <- labels |>
        dplyr::left_join(y_base, by = by_vars)
    } else {
      labels <- labels |>
        dplyr::mutate(y.position = y_base$y.position[[1]])
    }
  }

  # ── Stack multiple brackets within a panel ──────────────────────────────
  y_range <- if (!is.null(emm) && is.data.frame(emm_df) && nrow(emm_df) > 0L) {
    diff(range(emm_df$ymax, emm_df$ymin, na.rm = TRUE))
  } else {
    1
  }
  if (!is.finite(y_range) || y_range == 0) y_range <- 1
  step_abs <- step * y_range

  if (length(by_vars) > 0L) {
    labels <- labels |>
      dplyr::group_by(dplyr::across(dplyr::all_of(by_vars))) |>
      dplyr::mutate(
        .stack = dplyr::row_number() - 1L,
        y.position = .data$y.position + .stack * step_abs
      ) |>
      dplyr::ungroup() |>
      dplyr::select(-".stack")
  } else {
    labels <- labels |>
      dplyr::mutate(
        y.position = .data$y.position + (dplyr::row_number() - 1L) * step_abs
      )
  }

  # ── Apply y.adjust ──────────────────────────────────────────────────────
  n <- nrow(labels)
  if (length(y.adjust) == 1L) {
    labels <- dplyr::mutate(labels, y.position = .data$y.position + y.adjust)
  } else if (length(y.adjust) == n) {
    labels <- dplyr::mutate(labels, y.position = .data$y.position + y.adjust)
  } else {
    stop(
      "`y.adjust` must be length 1 or match the number of contrast rows (",
      n, ").",
      call. = FALSE
    )
  }

  attr(labels, "by.vars") <- by_vars
  labels
}

# ── Internal helpers ──────────────────────────────────────────────────────────

#' @keywords internal
.bb_split_contrast_groups <- function(df) {
  contrasts <- as.character(df$contrast)
  # Prefer ratio-style separator when present in any row
  if (any(grepl(" / ", contrasts, fixed = TRUE))) {
    sep <- " / "
  } else {
    sep <- " - "
  }

  df |>
    dplyr::mutate(
      group1 = stringr::str_trim(stringr::str_split_i(.data$contrast, sep, 1)),
      group2 = stringr::str_trim(stringr::str_split_i(.data$contrast, sep, 2))
    )
}

#' @keywords internal
.bb_add_effect_annotation <- function(df) {
  if ("d" %in% names(df)) {
    df |>
      dplyr::mutate(
        effect_annotation = sprintf("Cohen's d = %.2f", abs(.data$d))
      )
  } else if ("Cohen's d" %in% names(df)) {
    df |>
      dplyr::mutate(
        effect_annotation = sprintf("Cohen's d = %.2f", abs(.data[["Cohen's d"]]))
      )
  } else if ("IRR" %in% names(df)) {
    df |>
      dplyr::mutate(
        effect_annotation = sprintf("IRR = %.2f", .data$IRR)
      )
  } else if ("Odds Ratio" %in% names(df)) {
    df |>
      dplyr::mutate(
        effect_annotation = sprintf("Odds Ratio = %.2f", .data[["Odds Ratio"]])
      )
  } else {
    df |>
      dplyr::mutate(effect_annotation = "")
  }
}

#' @keywords internal
.bb_ypos_lookup <- function(emm_df, by_vars) {
  if (length(by_vars) == 0L) {
    tibble::tibble(
      y.position = max(emm_df$ymax, na.rm = TRUE)
    )
  } else {
    emm_df |>
      dplyr::group_by(dplyr::across(dplyr::all_of(by_vars))) |>
      dplyr::summarise(
        y.position = max(.data$ymax, na.rm = TRUE),
        .groups = "drop"
      )
  }
}
