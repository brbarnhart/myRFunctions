#' Add emmeans pairwise comparisons to a ggplot
#'
#' Add significance brackets from `pairs(emmeans(...))` to an existing ggplot
#' with `+`. Comparisons on the fill/colour grouping are drawn between dodged
#' bars at each x position; comparisons on the x-axis span the matching
#' x categories. By-variables in `pw` (for example `Sex`) are matched to
#' facets automatically.
#'
#' `dodge_width` must match the plot's `position_dodge(width = ...)`. When
#' omitted, it is inferred from dodge / jitter-dodge layers.
#'
#' Bracket height is `y.fun` of the plotted `y` values in each cell (default
#' [max], so brackets sit above jittered points), plus padding. Use
#' `y.fun = mean` if the plot shows bars of means without points.
#'
#' @param pw An `emmGrid` of pairwise contrasts (e.g. `pairs(emm)`), or a
#'   data frame with `contrast`, `p.value`, and any by-variable columns.
#' @param data Optional data used to compute y-positions and x / group
#'   levels. Defaults to the plot data.
#' @param mapping Optional aesthetic mapping used to resolve `x`, `y`, and
#'   the dodge grouping (`fill`, then `colour`, then `group`). Defaults to
#'   the plot mapping.
#' @param group Optional name of the dodged grouping variable (the factor
#'   being compared, e.g. `"Stim"`). Inferred from `fill` / `colour` /
#'   `group` aesthetics when omitted.
#' @param dodge_width Width of the dodge, matching
#'   `position_dodge(width = ...)`. `NULL` (default) infers it from existing
#'   layers, else `0.8`.
#' @param y.fun Function used to compute the bracket baseline from the
#'   plotted `y` values in each cell. Default [max].
#' @param y.adjust Extra vertical padding added after `y.fun`. `NULL`
#'   (default) uses 5% of the data y-range.
#' @param step Fraction of the y-range used to stack multiple brackets at
#'   the same x position. Default `0.08`.
#' @param hide.ns If `TRUE`, drop rows with `p.signif == "ns"`.
#' @param tip.length Passed to [ggpubr::stat_pvalue_manual()].
#' @param label Column used as the bracket label. Default `"p.signif"`.
#' @param model Optional fitted model forwarded to [bb_pairwise_labels()].
#' @param ... Additional arguments passed to [ggpubr::stat_pvalue_manual()].
#'
#' @return An object that can be added to a ggplot with `+`.
#' @export
#' @seealso [bb_pairwise_labels()], [bbmake_pairwise_plot()]
#' @examples
#' \dontrun{
#' pw <- emmeans(mod, ~ Stim | Sex * Diet * Satiety, type = "response") |>
#'   pairs(reverse = TRUE)
#'
#' ggplot(df, aes(x = interaction(Diet, Satiety), y = Breakpoint, fill = Stim)) +
#'   geom_bar(
#'     stat = "summary", fun = mean,
#'     width = 0.8, position = position_dodge(width = 0.8)
#'   ) +
#'   geom_point(position = position_jitterdodge(dodge.width = 0.8)) +
#'   facet_wrap(~ Sex) +
#'   bb_add_pairwise(pw, hide.ns = TRUE)
#' }
bb_add_pairwise <- function(
    pw,
    data = NULL,
    mapping = NULL,
    group = NULL,
    dodge_width = NULL,
    y.fun = max,
    y.adjust = NULL,
    step = 0.08,
    hide.ns = FALSE,
    tip.length = 0.01,
    label = "p.signif",
    model = NULL,
    ...
) {
  structure(
    list(
      pw = pw,
      data = data,
      mapping = mapping,
      group = rlang::enquo(group),
      dodge_width = dodge_width,
      y.fun = y.fun,
      y.adjust = y.adjust,
      step = step,
      hide.ns = hide.ns,
      tip.length = tip.length,
      label = label,
      model = model,
      extra = list(...)
    ),
    class = "bb_pairwise_layer"
  )
}

#' @export
print.bb_pairwise_layer <- function(x, ...) {
  cat("<bb_pairwise_layer>\n")
  cat("Add to a ggplot with `+ bb_add_pairwise(...)`.\n")
  invisible(x)
}

#' @importFrom ggplot2 ggplot_add
#' @export
ggplot_add.bb_pairwise_layer <- function(object, plot, ...) {
  sig <- .bb_pairwise_layer_data(object, plot)

  if (nrow(sig) == 0L) {
    return(plot)
  }

  layer_args <- c(
    list(
      data = sig,
      label = object$label,
      xmin = "xmin",
      xmax = "xmax",
      y.position = "y.position",
      tip.length = object$tip.length,
      inherit.aes = FALSE,
      hide.ns = FALSE
    ),
    object$extra
  )

  plot + do.call(ggpubr::stat_pvalue_manual, layer_args)
}

# ── Internal helpers ──────────────────────────────────────────────────────────

#' @keywords internal
.bb_pairwise_layer_data <- function(object, plot) {
  data <- object$data %||% plot$data
  if (!is.data.frame(data) || nrow(data) == 0L) {
    stop(
      "`bb_add_pairwise()` needs a data frame. Pass `data` or map data in ggplot().",
      call. = FALSE
    )
  }

  mapping <- .bb_plot_mapping(plot, object$mapping)
  if (is.null(mapping$x) || is.null(mapping$y)) {
    stop(
      "`bb_add_pairwise()` needs `x` and `y` aesthetics on the plot or in `mapping`.",
      call. = FALSE
    )
  }

  pw <- object$pw
  pw_table <- NULL
  if (is.data.frame(pw) && !inherits(pw, "emmGrid")) {
    pw_table <- pw
    pw <- NULL
  }

  labels <- bb_pairwise_labels(
    pw       = pw,
    pw_table = pw_table,
    model    = object$model,
    y.adjust = 0,
    step     = 0,
    hide.ns  = object$hide.ns
  )
  by_vars <- attr(labels, "by.vars") %||% character(0)
  by_vars <- by_vars[by_vars %in% names(labels)]

  if (nrow(labels) == 0L) {
    return(labels)
  }

  x_vals <- rlang::eval_tidy(mapping$x, data)
  y_vals <- rlang::eval_tidy(mapping$y, data)
  group_vals <- .bb_eval_group(object$group, mapping, data)

  x_levels <- .bb_discrete_levels(x_vals)
  group_levels <- .bb_discrete_levels(group_vals)

  contrast_levels <- unique(c(
    as.character(labels$group1),
    as.character(labels$group2)
  ))
  in_x <- length(x_levels) > 0L && all(contrast_levels %in% x_levels)
  in_group <- length(group_levels) > 0L && all(contrast_levels %in% group_levels)

  if (in_x) {
    mode <- "x"
  } else if (in_group) {
    mode <- "dodge"
  } else {
    stop(
      "Contrast groups (", paste(contrast_levels, collapse = ", "),
      ") were not found in the x-axis levels or the fill/colour/group aesthetic.",
      call. = FALSE
    )
  }

  x_pw <- tryCatch(
    rlang::eval_tidy(mapping$x, labels),
    error = function(e) NULL
  )

  if (mode == "dodge") {
    if (is.null(x_pw)) {
      stop(
        "Could not evaluate the plot's x mapping on the pairwise table. ",
        "The x aesthetic must use columns present in `pw` (for example ",
        "interaction(Diet, Satiety) when those are by-variables).",
        call. = FALSE
      )
    }
    x_num <- match(as.character(x_pw), x_levels)
    if (anyNA(x_num)) {
      missing <- unique(as.character(x_pw)[is.na(x_num)])
      stop(
        "Could not match pairwise rows to x-axis values: ",
        paste(missing, collapse = ", "),
        call. = FALSE
      )
    }

    dodge_width <- object$dodge_width %||% .bb_infer_dodge_width(plot)
    n <- length(group_levels)
    idx1 <- match(as.character(labels$group1), group_levels)
    idx2 <- match(as.character(labels$group2), group_levels)
    xmin <- x_num + .bb_dodge_offset(idx1, n, dodge_width)
    xmax <- x_num + .bb_dodge_offset(idx2, n, dodge_width)
    labels$xmin <- pmin(xmin, xmax)
    labels$xmax <- pmax(xmin, xmax)
    labels$.bb_x <- as.character(x_pw)
  } else {
    labels$xmin <- as.character(labels$group1)
    labels$xmax <- as.character(labels$group2)
    labels$.bb_x <- NA_character_
  }

  labels <- .bb_align_by_vars(labels, data, by_vars)
  labels <- .bb_add_y_position(
    labels,
    data = data,
    y_vals = y_vals,
    x_vals = x_vals,
    by_vars = by_vars,
    mode = mode,
    y.fun = object$y.fun,
    y.adjust = object$y.adjust,
    step = object$step
  )

  labels
}

#' @keywords internal
.bb_plot_mapping <- function(plot, override = NULL) {
  mapping <- plot$mapping
  if (is.null(mapping)) {
    mapping <- ggplot2::aes()
  }
  for (layer in plot$layers) {
    layer_map <- layer$mapping
    if (is.null(layer_map) || length(names(layer_map)) == 0L) {
      next
    }
    for (nm in names(layer_map)) {
      if (is.null(mapping[[nm]])) {
        mapping[[nm]] <- layer_map[[nm]]
      }
    }
  }
  if (!is.null(override) && length(override) > 0L) {
    for (nm in names(override)) {
      mapping[[nm]] <- override[[nm]]
    }
  }
  mapping
}

#' @keywords internal
.bb_eval_group <- function(group_quo, mapping, data) {
  if (!rlang::quo_is_null(group_quo)) {
    if (rlang::quo_is_symbol(group_quo) || rlang::quo_is_call(group_quo)) {
      return(rlang::eval_tidy(group_quo, data))
    }
    val <- rlang::eval_tidy(group_quo)
    if (is.character(val) && length(val) == 1L && val %in% names(data)) {
      return(data[[val]])
    }
    stop(
      "Could not interpret `group`. Pass a column name or a bare name.",
      call. = FALSE
    )
  }

  group_map <- mapping$fill %||% mapping$colour %||% mapping$group
  if (is.null(group_map)) {
    return(NULL)
  }
  rlang::eval_tidy(group_map, data)
}

#' @keywords internal
.bb_discrete_levels <- function(x) {
  if (is.null(x)) {
    return(character(0))
  }
  if (is.factor(x)) {
    lv <- levels(x)
    return(lv[lv %in% unique(as.character(x))])
  }
  unique(as.character(x))
}

#' @keywords internal
.bb_dodge_offset <- function(index, n, dodge_width) {
  (index - (n + 1) / 2) * (dodge_width / n)
}

#' @keywords internal
.bb_infer_dodge_width <- function(plot, default = 0.8) {
  for (layer in rev(plot$layers)) {
    pos <- layer$position
    if (inherits(pos, "PositionDodge") || inherits(pos, "PositionDodge2")) {
      w <- pos$width
      if (!is.null(w) && length(w) == 1L && is.finite(w)) {
        return(w)
      }
    }
    if (inherits(pos, "PositionJitterdodge")) {
      w <- pos$dodge.width
      if (!is.null(w) && length(w) == 1L && is.finite(w)) {
        return(w)
      }
    }
  }
  default
}

#' @keywords internal
.bb_align_by_vars <- function(labels, data, by_vars) {
  for (v in by_vars) {
    if (!v %in% names(data) || !v %in% names(labels)) {
      next
    }
    if (is.factor(data[[v]])) {
      labels[[v]] <- factor(as.character(labels[[v]]), levels = levels(data[[v]]))
    } else if (is.character(data[[v]])) {
      labels[[v]] <- as.character(labels[[v]])
    }
  }
  labels
}

#' @keywords internal
.bb_add_y_position <- function(
    labels,
    data,
    y_vals,
    x_vals,
    by_vars,
    mode,
    y.fun,
    y.adjust,
    step
) {
  y_vals <- as.numeric(y_vals)
  y_range <- diff(range(y_vals, na.rm = TRUE))
  if (!is.finite(y_range) || y_range == 0) {
    y_range <- 1
  }
  pad <- if (is.null(y.adjust)) 0.05 * y_range else y.adjust
  step_abs <- step * y_range

  if ("y.position" %in% names(labels)) {
    labels$y.position <- NULL
  }

  plot_df <- tibble::as_tibble(data)
  plot_df$.bb_y <- y_vals
  plot_df$.bb_x <- as.character(x_vals)

  keys <- if (mode == "dodge") {
    unique(c(".bb_x", by_vars[by_vars %in% names(plot_df)]))
  } else {
    by_vars[by_vars %in% names(plot_df)]
  }

  y_fun_na <- function(v) y.fun(v[!is.na(v)])

  y_base <- if (length(keys) == 0L) {
    tibble::tibble(y.position = y_fun_na(plot_df$.bb_y))
  } else {
    plot_df |>
      dplyr::group_by(dplyr::across(dplyr::all_of(keys))) |>
      dplyr::summarise(
        y.position = y_fun_na(.data$.bb_y),
        .groups = "drop"
      )
  }

  if (length(keys) == 0L) {
    labels$y.position <- y_base$y.position[[1]]
  } else {
    join_keys <- keys[keys %in% names(labels)]
    if (length(join_keys) == 0L) {
      labels$y.position <- max(y_base$y.position, na.rm = TRUE)
    } else {
      labels <- dplyr::left_join(
        labels,
        y_base,
        by = join_keys
      )
    }
  }

  if (anyNA(labels$y.position)) {
    fallback <- max(y_vals, na.rm = TRUE)
    labels$y.position[is.na(labels$y.position)] <- fallback
  }

  labels$y.position <- labels$y.position + pad

  stack_vars <- if (mode == "dodge") {
    unique(c(".bb_x", by_vars[by_vars %in% names(labels)]))
  } else {
    by_vars[by_vars %in% names(labels)]
  }
  stack_vars <- stack_vars[stack_vars %in% names(labels)]

  if (length(stack_vars) > 0L) {
    labels <- labels |>
      dplyr::group_by(dplyr::across(dplyr::all_of(stack_vars))) |>
      dplyr::mutate(
        y.position = .data$y.position + (dplyr::row_number() - 1L) * step_abs
      ) |>
      dplyr::ungroup()
  } else {
    labels <- labels |>
      dplyr::mutate(
        y.position = .data$y.position + (dplyr::row_number() - 1L) * step_abs
      )
  }

  dplyr::select(labels, -dplyr::any_of(".bb_x"))
}
