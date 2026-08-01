#' Bar plot with mean CI, optional jittered points, and zero line
#'
#' Thin wrapper around [ggpubr::ggbarplot()] that applies a consistent theme
#' and fill palette, optionally overlays individual points (jitter-dodged to
#' match the bars), and optionally draws a horizontal line at y = 0.
#'
#' @param data A data frame.
#' @param x,y,fill Character names of columns mapped to x, y, and fill
#'   (same convention as [ggpubr::ggbarplot()]).
#' @param dodge_width Width passed to [ggplot2::position_dodge()] and to the
#'   bar width.
#' @param jitter.width Horizontal jitter for individual points
#'   ([ggplot2::position_jitterdodge()]).
#' @param point.size,point.alpha,point.color Appearance of individual points.
#' @param add_hline If `TRUE`, draw a horizontal line at y = 0.
#' @param add_points If `TRUE`, overlay jittered individual points aligned to
#'   the dodged bars.
#' @param ... Additional arguments passed to [ggpubr::ggbarplot()].
#'
#' @return A ggplot object.
#' @export
#' @examples
#' \dontrun{
#' bb_ggbarplot(df, x = "group", y = "score", fill = "condition") +
#'   labs(y = "Score")
#' }
bb_ggbarplot <- function(
    data,
    x,
    y,
    fill,
    dodge_width = 0.9,
    jitter.width = 0.2,
    point.size = 1.5,
    point.alpha = 0.6,
    point.color = "gray10",
    add_hline = TRUE,
    add_points = TRUE,
    ...
) {
  p <- ggpubr::ggbarplot(
    data = data,
    x = x,
    y = y,
    fill = fill,
    add = "mean_ci",
    error.plot = "errorbar",
    position = ggplot2::position_dodge(dodge_width),
    width = dodge_width,
    add.params = list(width = 0.2),
    ggtheme = NULL,
    palette = NULL,
    ...
  ) +
    ggplot2::theme_bw() +
    ggplot2::theme(text = ggplot2::element_text(size = 18)) +
    ggplot2::scale_fill_brewer(palette = "Set1")

  # Jittered points (aligned to the dodged bars)
  if (isTRUE(add_points)) {
    p <- p +
      ggplot2::geom_point(
        ggplot2::aes(group = .data[[fill]]),
        position = ggplot2::position_jitterdodge(
          dodge.width = dodge_width,
          jitter.width = jitter.width,
          jitter.height = 0
        ),
        color = point.color,
        size = point.size,
        alpha = point.alpha
      )
  }

  # Horizontal line at zero
  if (isTRUE(add_hline)) {
    p <- p +
      ggplot2::geom_hline(
        yintercept = 0,
        linetype = "solid",
        color = "black"
      )
  }

  p
}
