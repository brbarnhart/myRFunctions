#' Create a pairwise comparison plot from emmeans
#'
#' High-level convenience wrapper around [bb_emm_df()] and
#' [bb_add_pairwise()]. Builds a ggplot with points, T-capped 95% CI
#' error bars, optional connecting lines, faceting by by-variables,
#' significance brackets, and a single effect-size annotation per panel.
#'
#' P-values default to Holm across every pairwise test
#' (`adjust = "none"`, `cross.adjust = "holm"`). The same message as
#' [bbmake_pairwise_table()] reports the count, for example
#' `Holm across 6 comparisons`. A prebuilt `pw_table` is used as-is.
#'
#' For custom layouts, call the helpers yourself and add intervals /
#' brackets with [bb_add_errorbar()] and [bb_add_pairwise()].
#'
#' @param emm An `emmGrid` from [emmeans::emmeans()].
#' @param x Bare name or string of the x-axis factor. Default: the sole
#'   primary variable on `emm` (error if ambiguous).
#' @param facets A facet formula (e.g. `Diet ~ Sex` or `~ Group`), or `NULL`
#'   to auto-facet by by-variables (`facet_wrap` for one by-var, `facet_grid`
#'   for two, `facet_wrap` with multiple vars otherwise). Set to `FALSE` to
#'   suppress faceting.
#' @param pw Optional result of `pairs(emm)`. Computed from `emm` when both
#'   `pw` and `pw_table` are `NULL`, using `reverse`, `adjust`, and
#'   `cross.adjust`. Pass this only for a custom contrast family.
#' @param pw_table Optional output of [bbmake_pairwise_table()]. Built from
#'   `pw` when omitted. If you already have a table, p-values are used as-is
#'   (`adjust` / `cross.adjust` are not re-applied).
#' @param model Optional fitted model used for Gaussian effect sizes
#'   (Cohen's d). Recovered from `emm` / `pw` when possible; not needed for
#'   IRR / odds-ratio plots.
#' @param reverse If `TRUE` (the default), reverse pairwise direction so
#'   later factor levels are in the numerator. Ignored when `pw` or
#'   `pw_table` is supplied.
#' @param adjust Within-`by`-group multiplicity adjustment used when `pw`
#'   is computed from `emm`. Default `"none"`. Ignored when `pw` or
#'   `pw_table` is supplied, and ignored when `cross.adjust` is not `"none"`.
#' @param cross.adjust Adjustment for **all** pairwise tests on the plot
#'   as one family (every Stim pair in every Sex × Diet cell). Default
#'   `"holm"`. Set `"none"` for no family-wide correction, or
#'   `adjust = "tukey", cross.adjust = "none"` for Tukey within each
#'   by-group. Ignored when `pw_table` is supplied.
#' @param connect If `TRUE`, draw lines connecting points within each panel
#'   along `x` (grouped by by-variables).
#' @param hide.ns Passed to [bb_add_pairwise()] and [bb_pairwise_labels()].
#' @param y.adjust Vertical nudge for significance brackets.
#' @param step Stacking step for multiple brackets within a panel (fraction
#'   of y-range). See [bb_add_pairwise()].
#' @param annotate_effect If `TRUE`, draw one effect-size label per panel
#'   (first contrast in that panel). When `y_expand` is `NULL`, this also
#'   uses extra top padding so the annotation is not cramped against the
#'   panel edge.
#' @param y_expand Y-axis range expansion passed to
#'   [ggplot2::scale_y_continuous()]. `NULL` (default) uses
#'   `expansion(mult = c(0.02, 0.28))` when `annotate_effect` is `TRUE`
#'   and `expansion(mult = c(0.02, 0.18))` otherwise. A numeric of length
#'   1 or 2 is treated as multiplicative padding
#'   (`expansion(mult = y_expand)`). Otherwise pass the result of
#'   [ggplot2::expansion()].
#' @param hjust,vjust Position adjustments for the effect-size annotation
#'   (placed at `x = Inf`, `y = Inf`). Text size is not hardcoded: it
#'   follows the plot theme (ggplot2 4.0+: `theme_*(base_size)` or
#'   `theme(geom = element_geom(fontsize = ...))`).
#' @param point_size Size passed to [ggplot2::geom_point()] (`size`).
#' @param linewidth Line width for error bars and connecting lines.
#' @param errorbar_width Horizontal width of the T-caps on
#'   [ggplot2::geom_errorbar()].
#' @param interval `"ci"` (default) uses emmeans confidence limits (the
#'   backtransformed limits when the grid was created with
#'   `type = "response"`). `"sem"` uses the estimate ± `sem_mult` times
#'   the emmeans SE. See [bb_add_errorbar()].
#' @param sem_mult Multiplier for `interval = "sem"`. Default `1`.
#' @param show_zero_line If `TRUE`, draw a horizontal line at y = 0.
#' @param color Color for points / lines / error bars.
#' @param ... Additional arguments reserved for future use (currently unused).
#'
#' @return A ggplot object. Add labels/themes with `+` as usual.
#' @export
#' @examples
#' \dontrun{
#' emm <- emmeans(mod, ~ Stim | Diet * Sex, type = "response")
#'
#' bbmake_pairwise_plot(emm) +
#'   labs(x = "Stimulation", y = "Breakpoint (active pokes)")
#'
#' # Unadjusted p values
#' bbmake_pairwise_plot(emm, cross.adjust = "none")
#'
#' # Tukey within each by-group
#' bbmake_pairwise_plot(emm, adjust = "tukey", cross.adjust = "none")
#'
#' # Custom pairs still work. The default Holm adjustment is applied
#' # unless cross.adjust = "none" or pw_table is supplied.
#' bbmake_pairwise_plot(emm, pw = pairs(emm, reverse = TRUE, adjust = "tukey"))
#'
#' # Extra top padding (multiplicative c(bottom, top))
#' bbmake_pairwise_plot(emm, y_expand = c(0.02, 0.4))
#'
#' # Deeper control with helpers
#' df <- bb_emm_df(emm)
#' ggplot(df, aes(x = Stim, y = y, fill = Diet)) +
#'   geom_point() +
#'   facet_wrap(~ Sex) +
#'   bb_add_errorbar(emm) +
#'   bb_add_pairwise(pairs(emm))
#' }
bbmake_pairwise_plot <- function(
    emm,
    x = NULL,
    facets = NULL,
    pw = NULL,
    pw_table = NULL,
    model = NULL,
    reverse = TRUE,
    adjust = "none",
    cross.adjust = "holm",
    connect = TRUE,
    hide.ns = FALSE,
    y.adjust = 0,
    step = 0.08,
    annotate_effect = TRUE,
    y_expand = NULL,
    hjust = 1.05,
    vjust = 1.3,
    point_size = 2,
    linewidth = 0.6,
    errorbar_width = 0.15,
    interval = c("ci", "sem"),
    sem_mult = 1,
    show_zero_line = TRUE,
    color = "black",
    ...
) {
  if (!inherits(emm, "emmGrid")) {
    stop("`emm` must be an emmGrid from emmeans().", call. = FALSE)
  }

  interval <- match.arg(interval)
  emm_df   <- .bb_apply_interval(bb_emm_df(emm), interval, sem_mult)
  pri_vars <- attr(emm_df, "pri.vars") %||% character(0)
  by_vars  <- attr(emm_df, "by.vars")  %||% character(0)

  # ── Resolve x (NSE: bare name or string) ────────────────────────────────
  x_name <- .bb_resolve_x(rlang::enquo(x), pri_vars, emm_df)
  if (!x_name %in% names(emm_df)) {
    stop("x variable '", x_name, "' not found in emmeans grid.", call. = FALSE)
  }

  # Precompute line groups so aes() stays simple
  if (length(by_vars) > 0L) {
    emm_df$.line_group <- interaction(emm_df[by_vars], drop = TRUE)
  } else {
    emm_df$.line_group <- 1L
  }

  # ── Resolve pairwise contrasts ──────────────────────────────────────────
  # NULL keeps whatever adjustment is already stored on a user-supplied pw
  # when cross.adjust is "none". The table default is adjust = "none".
  adjust_for_table <- NULL
  if (is.null(pw) && is.null(pw_table)) {
    method <- if (isTRUE(reverse)) "revpairwise" else "pairwise"
    pw <- emmeans::contrast(emm, method = method, adjust = adjust)
    adjust_for_table <- adjust
  }
  if (is.null(model)) {
    model <- .bb_recover_model(emm, pw)
  }
  if (is.null(pw_table) && inherits(pw, "emmGrid")) {
    pw_table <- bbmake_pairwise_table(
      pw,
      model = model,
      adjust = adjust_for_table,
      cross.adjust = cross.adjust
    )
  }
  pw_for_add <- if (!is.null(pw_table)) pw_table else pw

  # ── Base plot ───────────────────────────────────────────────────────────
  p <- ggplot2::ggplot(emm_df) +
    ggplot2::geom_errorbar(
      ggplot2::aes(
        x    = .data[[x_name]],
        ymin = .data$ymin,
        ymax = .data$ymax
      ),
      width     = errorbar_width,
      linewidth = linewidth,
      color     = color
    ) +
    ggplot2::geom_point(
      ggplot2::aes(
        x = .data[[x_name]],
        y = .data$y
      ),
      size  = point_size,
      color = color
    )

  if (isTRUE(connect)) {
    p <- p +
      ggplot2::geom_line(
        ggplot2::aes(
          x     = .data[[x_name]],
          y     = .data$y,
          group = .data$.line_group
        ),
        linewidth = linewidth,
        color     = color
      )
  }

  if (isTRUE(show_zero_line)) {
    p <- p +
      ggplot2::geom_hline(
        yintercept = 0,
        linetype   = "solid",
        color      = "gray50",
        linewidth  = 0.6
      )
  }

  # ── Facets ──────────────────────────────────────────────────────────────
  p <- .bb_add_facets(p, facets = facets, by_vars = by_vars)

  # ── Significance brackets ───────────────────────────────────────────────
  p <- p +
    bb_add_pairwise(
      pw_for_add,
      model    = model,
      hide.ns  = hide.ns,
      y.adjust = y.adjust,
      step     = step
    )

  # ── Effect annotation (one per panel) ───────────────────────────────────
  if (isTRUE(annotate_effect)) {
    sig <- bb_pairwise_labels(
      pw       = pw,
      emm      = emm_df,
      pw_table = pw_table,
      model    = model,
      y.adjust = y.adjust,
      step     = step,
      hide.ns  = hide.ns
    )
    sig_by <- attr(sig, "by.vars") %||% character(0)
    if (nrow(sig) > 0L) {
      ann <- .bb_effect_annotation_df(sig, sig_by)
      if (nrow(ann) > 0L && any(nzchar(ann$effect_annotation))) {
        p <- p +
          ggplot2::geom_text(
            data = ann,
            ggplot2::aes(label = .data$effect_annotation),
            x = Inf, y = Inf,
            hjust = hjust, vjust = vjust,
            inherit.aes = FALSE
          )
      }
    }
  }

  p +
    ggplot2::scale_y_continuous(
      expand = .bb_resolve_y_expand(y_expand, annotate_effect)
    )
}

# ── Internal helpers ──────────────────────────────────────────────────────────

#' @keywords internal
.bb_resolve_y_expand <- function(y_expand, annotate_effect) {
  if (is.null(y_expand)) {
    return(ggplot2::expansion(
      mult = c(0.02, if (isTRUE(annotate_effect)) 0.28 else 0.18)
    ))
  }
  if (is.numeric(y_expand) && length(y_expand) %in% c(1L, 2L)) {
    return(ggplot2::expansion(mult = y_expand))
  }
  y_expand
}

#' @keywords internal
.bb_resolve_x <- function(x_quo, pri_vars, emm_df) {
  if (!rlang::quo_is_null(x_quo)) {
    if (rlang::quo_is_symbol(x_quo)) {
      return(rlang::as_name(x_quo))
    }
    val <- rlang::eval_tidy(x_quo)
    if (is.character(val) && length(val) == 1L) {
      return(val)
    }
    stop("Could not interpret `x`. Pass a bare name or a string.", call. = FALSE)
  }

  if (length(pri_vars) == 1L) return(pri_vars[[1]])
  if (length(pri_vars) == 0L) {
    candidates <- setdiff(
      names(emm_df),
      c("y", "ymin", "ymax", "SE", "df", "emmean", "response",
        "lower.CL", "upper.CL", "asymp.LCL", "asymp.UCL",
        "t.ratio", "z.ratio", "p.value", ".line_group")
    )
    if (length(candidates) >= 1L) return(candidates[[1]])
    stop("Could not infer x-axis variable; please supply `x`.", call. = FALSE)
  }
  stop(
    "Multiple primary variables on emm (",
    paste(pri_vars, collapse = ", "),
    "); please supply `x`.",
    call. = FALSE
  )
}

#' @keywords internal
.bb_add_facets <- function(p, facets, by_vars) {
  if (isFALSE(facets)) return(p)

  if (!is.null(facets)) {
    if (inherits(facets, "formula")) {
      rhs_only <- length(facets) == 2L
      if (rhs_only) {
        return(p + ggplot2::facet_wrap(facets))
      }
      return(p + ggplot2::facet_grid(facets))
    }
    stop("`facets` must be a formula, NULL, or FALSE.", call. = FALSE)
  }

  # Auto from by_vars
  if (length(by_vars) == 0L) return(p)
  if (length(by_vars) == 1L) {
    return(p + ggplot2::facet_wrap(stats::as.formula(paste("~", by_vars[[1]]))))
  }
  if (length(by_vars) == 2L) {
    fml <- stats::as.formula(paste(by_vars[[1]], "~", by_vars[[2]]))
    return(p + ggplot2::facet_grid(fml))
  }
  # 3+ by-vars: wrap on all
  fml <- stats::as.formula(paste("~", paste(by_vars, collapse = " + ")))
  p + ggplot2::facet_wrap(fml)
}

#' @keywords internal
.bb_effect_annotation_df <- function(sig, by_vars) {
  if (length(by_vars) == 0L) {
    sig |>
      dplyr::slice(1) |>
      dplyr::select("effect_annotation")
  } else {
    sig |>
      dplyr::group_by(dplyr::across(dplyr::all_of(by_vars))) |>
      dplyr::slice(1) |>
      dplyr::ungroup() |>
      dplyr::select(dplyr::all_of(c(by_vars, "effect_annotation")))
  }
}
