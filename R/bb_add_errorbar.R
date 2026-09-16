#' Add emmeans confidence-interval or SEM error bars to a ggplot
#'
#' Overlay model-based intervals from [emmeans::emmeans()] on an existing
#' ggplot with `+`. This is the interval counterpart to [bb_add_pairwise()]:
#' dodged `fill` / `colour` groups get dodged bars, and by-variables (for
#' example `Sex`) are matched to facets automatically.
#'
#' **Backtransformation.** If you created the grid with
#' `emmeans(..., type = "response")` on a log, logit, or other link-scale
#' model, `interval = "ci"` (the default) uses the confidence limits emmeans
#' computed on the link scale and then backtransformed. Those intervals are
#' generally asymmetric on the response scale. That is the same source
#' [bbmake_pairwise_plot()] uses. `interval = "sem"` is the estimate ±
#' `sem_mult * SE` on the plotted scale (delta-method SE after
#' `type = "response"`) and is symmetric; it is not a backtransformed CI.
#'
#' @param emm An `emmGrid` from [emmeans::emmeans()], or a data frame /
#'   result of [bb_emm_df()] with `y`, `ymin`, `ymax` (and `SE` when
#'   `interval = "sem"`).
#' @param interval `"ci"` (default) uses emmeans confidence limits.
#'   `"sem"` uses the estimate ± `sem_mult` times the emmeans SE.
#' @param sem_mult Multiplier for `interval = "sem"`. Default `1` (one
#'   standard error). Ignored for `"ci"`.
#' @param data Optional data used only to copy factor levels onto the
#'   emmeans grid so `interaction()` and dodge order match the plot.
#'   Defaults to the plot data.
#' @param mapping Optional aesthetic mapping used to resolve `x` and the
#'   dodge grouping (`fill`, then `colour`, then `group`). Defaults to the
#'   plot mapping.
#' @param group Optional name of the dodged grouping variable (e.g.
#'   `"Stim"`). Inferred from `fill` / `colour` / `group` when omitted.
#' @param dodge_width Width of the dodge, matching
#'   `position_dodge(width = ...)`. `NULL` (default) infers it from existing
#'   layers, else `0.8`.
#' @param width Horizontal width of the T-caps on
#'   [ggplot2::geom_errorbar()].
#' @param linewidth Line width for the error bars.
#' @param color Color for the error bars.
#' @param ... Additional arguments passed to [ggplot2::geom_errorbar()].
#'
#' @return An object that can be added to a ggplot with `+`.
#' @export
#' @seealso [bb_add_pairwise()], [bb_emm_df()], [bbmake_pairwise_plot()]
#' @examples
#' \dontrun{
#' emm <- emmeans(mod, ~ Stim | Sex * Diet * Satiety, type = "response")
#'
#' ggplot(df, aes(x = interaction(Diet, Satiety), y = Breakpoint, fill = Stim)) +
#'   geom_bar(
#'     stat = "summary", fun = mean,
#'     width = 0.8, position = position_dodge(width = 0.8)
#'   ) +
#'   facet_wrap(~ Sex) +
#'   bb_add_errorbar(emm) +
#'   bb_add_pairwise(pairs(emm, reverse = TRUE), hide.ns = TRUE)
#' }
bb_add_errorbar <- function(
    emm,
    interval = c("ci", "sem"),
    sem_mult = 1,
    data = NULL,
    mapping = NULL,
    group = NULL,
    dodge_width = NULL,
    width = 0.2,
    linewidth = 0.6,
    color = "black",
    ...
) {
  interval <- match.arg(interval)
  structure(
    list(
      emm = emm,
      interval = interval,
      sem_mult = sem_mult,
      data = data,
      mapping = mapping,
      group = rlang::enquo(group),
      dodge_width = dodge_width,
      width = width,
      linewidth = linewidth,
      color = color,
      extra = list(...)
    ),
    class = "bb_errorbar_layer"
  )
}

#' @export
print.bb_errorbar_layer <- function(x, ...) {
  cat("<bb_errorbar_layer>\n")
  cat("Add to a ggplot with `+ bb_add_errorbar(...)`.\n")
  invisible(x)
}

#' @importFrom ggplot2 ggplot_add
#' @export
ggplot_add.bb_errorbar_layer <- function(object, plot, ...) {
  df <- .bb_errorbar_layer_data(object, plot)
  mapping <- .bb_plot_mapping(plot, object$mapping)
  if (is.null(mapping$x)) {
    stop(
      "`bb_add_errorbar()` needs an `x` aesthetic on the plot or in `mapping`.",
      call. = FALSE
    )
  }

  group_map <- .bb_group_map(object$group, mapping)
  dodge <- !is.null(group_map) && !.bb_same_mapping(mapping$x, group_map)
  dodge_width <- object$dodge_width %||% .bb_infer_dodge_width(plot)

  if (dodge) {
    layer_aes <- ggplot2::aes(
      x = !!mapping$x,
      ymin = .data$ymin,
      ymax = .data$ymax,
      group = !!group_map
    )
    pos <- ggplot2::position_dodge(width = dodge_width)
  } else {
    layer_aes <- ggplot2::aes(
      x = !!mapping$x,
      ymin = .data$ymin,
      ymax = .data$ymax
    )
    pos <- ggplot2::position_identity()
  }

  layer_args <- c(
    list(
      data = df,
      mapping = layer_aes,
      position = pos,
      inherit.aes = FALSE,
      width = object$width,
      linewidth = object$linewidth,
      color = object$color,
      na.rm = TRUE
    ),
    object$extra
  )

  plot + do.call(ggplot2::geom_errorbar, layer_args)
}

# ── Internal helpers ──────────────────────────────────────────────────────────

#' @keywords internal
.bb_errorbar_layer_data <- function(object, plot) {
  df <- bb_emm_df(object$emm)
  df <- .bb_apply_interval(df, object$interval, object$sem_mult)

  plot_data <- object$data %||% plot$data
  .bb_align_to_plot_data(df, plot_data)
}

#' @keywords internal
#' @noRd
.bb_apply_interval <- function(df, interval = "ci", sem_mult = 1) {
  interval <- match.arg(interval, c("ci", "sem"))
  if (identical(interval, "ci")) {
    return(df)
  }
  if (!"SE" %in% names(df)) {
    stop(
      "`interval = \"sem\"` needs an `SE` column from emmeans.",
      call. = FALSE
    )
  }
  if (!is.numeric(sem_mult) || length(sem_mult) != 1L ||
      is.na(sem_mult) || sem_mult < 0) {
    stop("`sem_mult` must be a single non-negative number.", call. = FALSE)
  }
  keep <- attributes(df)
  df$ymin <- df$y - sem_mult * df$SE
  df$ymax <- df$y + sem_mult * df$SE
  for (nm in c("pri.vars", "by.vars", "mean_col", "lower_col", "upper_col")) {
    if (!is.null(keep[[nm]])) {
      attr(df, nm) <- keep[[nm]]
    }
  }
  df
}

#' @keywords internal
.bb_group_map <- function(group_quo, mapping) {
  if (!rlang::quo_is_null(group_quo)) {
    if (rlang::quo_is_symbol(group_quo) || rlang::quo_is_call(group_quo)) {
      return(group_quo)
    }
    val <- rlang::eval_tidy(group_quo)
    if (is.character(val) && length(val) == 1L) {
      return(rlang::new_quosure(rlang::expr(.data[[!!val]])))
    }
    stop(
      "Could not interpret `group`. Pass a column name or a bare name.",
      call. = FALSE
    )
  }
  mapping$fill %||% mapping$colour %||% mapping$group
}

#' @keywords internal
.bb_same_mapping <- function(a, b) {
  if (is.null(a) || is.null(b)) {
    return(FALSE)
  }
  identical(rlang::quo_get_expr(a), rlang::quo_get_expr(b))
}

#' @keywords internal
.bb_align_to_plot_data <- function(emm_df, plot_data) {
  if (!is.data.frame(plot_data) || nrow(plot_data) == 0L) {
    return(emm_df)
  }
  for (nm in intersect(names(emm_df), names(plot_data))) {
    src <- plot_data[[nm]]
    if (is.factor(src)) {
      emm_df[[nm]] <- factor(as.character(emm_df[[nm]]), levels = levels(src))
    }
  }
  emm_df
}
