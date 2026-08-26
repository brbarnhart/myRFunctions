#' Create a pairwise comparison plot from emmeans
#'
#' High-level convenience wrapper around [bb_emm_df()] and
#' [bb_pairwise_labels()]. Builds a ggplot with points, T-capped 95% CI
#' error bars, optional connecting lines, faceting by by-variables,
#' significance brackets, and a single effect-size annotation per panel.
#'
#' For custom layouts, call the helpers yourself and layer geoms manually.
#'
#' @param emm An `emmGrid` from [emmeans::emmeans()].
#' @param x Bare name or string of the x-axis factor. Default: the sole
#'   primary variable on `emm` (error if ambiguous).
#' @param facets A facet formula (e.g. `Diet ~ Sex` or `~ Group`), or `NULL`
#'   to auto-facet by by-variables (`facet_wrap` for one by-var, `facet_grid`
#'   for two, `facet_wrap` with multiple vars otherwise). Set to `FALSE` to
#'   suppress faceting.
#' @param pw Optional result of `pairs(emm)`. Computed automatically when
#'   both `pw` and `pw_table` are `NULL`.
#' @param pw_table Optional output of [bbmake_pairwise_table()]. Built from
#'   `pw` when omitted.
#' @param model Optional fitted model forwarded for effect-size calculation.
#' @param connect If `TRUE`, draw lines connecting points within each panel
#'   along `x` (grouped by by-variables).
#' @param hide.ns Passed to [bb_pairwise_labels()].
#' @param y.adjust Vertical nudge for significance brackets.
#' @param step Stacking step for multiple brackets within a panel (fraction
#'   of y-range). See [bb_pairwise_labels()].
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
#'   (placed at `x = Inf`, `y = Inf`).
#' @param point_size Size passed to [ggplot2::geom_point()] (`size`).
#' @param linewidth Line width for error bars and connecting lines.
#' @param errorbar_width Horizontal width of the T-caps on
#'   [ggplot2::geom_errorbar()].
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
#' # Custom pairs
#' bbmake_pairwise_plot(emm, pw = pairs(emm, reverse = TRUE, adjust = "tukey"))
#'
#' # Extra top padding (multiplicative c(bottom, top))
#' bbmake_pairwise_plot(emm, y_expand = c(0.02, 0.4))
#'
#' # Deeper control with helpers
#' df  <- bb_emm_df(emm)
#' sig <- bb_pairwise_labels(pairs(emm), emm = emm)
#' ggplot(df, aes(x = Stim, y = y, ymin = ymin, ymax = ymax)) +
#'   geom_errorbar(width = 0.2) +
#'   geom_point() +
#'   facet_grid(Diet ~ Sex) +
#'   ggpubr::stat_pvalue_manual(
#'     sig, label = "p.signif",
#'     xmin = "group1", xmax = "group2",
#'     y.position = "y.position"
#'   )
#' }
bbmake_pairwise_plot <- function(
    emm,
    x = NULL,
    facets = NULL,
    pw = NULL,
    pw_table = NULL,
    model = NULL,
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
    show_zero_line = TRUE,
    color = "black",
    ...
) {
  if (!inherits(emm, "emmGrid")) {
    stop("`emm` must be an emmGrid from emmeans().", call. = FALSE)
  }

  emm_df   <- bb_emm_df(emm)
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

  # ── Resolve pairwise labels ─────────────────────────────────────────────
  if (is.null(pw) && is.null(pw_table)) {
    # pairs() is an S3 method registered by emmeans, not an exported object
    pw <- emmeans::contrast(emm, method = "pairwise")
  }

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
  if (nrow(sig) > 0L && !all(is.na(sig$y.position))) {
    p <- p +
      ggpubr::stat_pvalue_manual(
        data       = sig,
        label      = "p.signif",
        xmin       = "group1",
        xmax       = "group2",
        y.position = "y.position",
        hide.ns    = FALSE,
        tip.length = 0.01
      )
  }

  # ── Effect annotation (one per panel) ───────────────────────────────────
  if (isTRUE(annotate_effect) && nrow(sig) > 0L) {
    ann <- .bb_effect_annotation_df(sig, sig_by)
    if (nrow(ann) > 0L && any(nzchar(ann$effect_annotation))) {
      p <- p +
        ggplot2::geom_text(
          data = ann,
          ggplot2::aes(label = .data$effect_annotation),
          x = Inf, y = Inf,
          hjust = hjust, vjust = vjust,
          size = 4.2,
          inherit.aes = FALSE
        )
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
