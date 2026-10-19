# Effect-size labels and significance stars are drawn as geoms, so a
# theme added after the plot is built (`+ theme_bw(base_size = 18)`)
# has to be read while the layer is set up. Built-in complete themes
# size axis tick labels at 0.8 * base_size. These geoms copy that
# resolved tick size. When tick labels are blank, they use the same
# 0.8 ratio of the theme geom fontsize, which complete themes set
# from base_size.

.bb_theme_sets_geom_text <- function() {
  utils::packageVersion("ggplot2") >= "4.0.0"
}

# Millimetres, the unit geom_text() uses for size.
.bb_ticklabel_size_mm <- function(theme) {
  if (is.null(theme)) {
    return(NULL)
  }

  axis <- tryCatch(
    ggplot2::calc_element("axis.text", theme),
    error = function(e) NULL
  )
  axis_pt <- if (inherits(axis, "element_text")) axis$size else NULL
  if (is.numeric(axis_pt) && !inherits(axis_pt, "rel") &&
      length(axis_pt) == 1L && is.finite(axis_pt)) {
    return(axis_pt / ggplot2::.pt)
  }

  geom_el <- tryCatch(
    ggplot2::calc_element("geom", theme),
    error = function(e) NULL
  )
  fontsize <- if (is.null(geom_el)) NULL else geom_el@fontsize
  if (is.numeric(fontsize) && length(fontsize) == 1L && is.finite(fontsize)) {
    return(0.8 * fontsize)
  }
  NULL
}

.bb_apply_ticklabel_size <- function(data, theme, aesthetic) {
  if (is.null(data) || nrow(data) == 0L) {
    return(data)
  }
  size_mm <- .bb_ticklabel_size_mm(theme)
  if (is.null(size_mm)) {
    return(data)
  }
  data[[aesthetic]] <- size_mm
  data
}

.bb_GeomTicklabelText <- ggplot2::ggproto(
  "GeomTicklabelText",
  ggplot2::GeomText,
  use_defaults = function(self, data, params = list(), modifiers = ggplot2::aes(), ...) {
    data <- ggplot2::ggproto_parent(ggplot2::GeomText, self)$use_defaults(
      data, params, modifiers, ...
    )
    if (!is.null(params$size) || "size" %in% names(modifiers)) {
      return(data)
    }
    .bb_apply_ticklabel_size(data, list(...)$theme, "size")
  }
)

.bb_GeomBracket <- getFromNamespace("GeomBracket", "ggpubr")

.bb_GeomTicklabelBracket <- ggplot2::ggproto(
  "GeomTicklabelBracket",
  .bb_GeomBracket,
  use_defaults = function(self, data, params = list(), modifiers = ggplot2::aes(), ...) {
    data <- ggplot2::ggproto_parent(.bb_GeomBracket, self)$use_defaults(
      data, params, modifiers, ...
    )
    if (!is.null(params$label.size) || "label.size" %in% names(modifiers)) {
      return(data)
    }
    .bb_apply_ticklabel_size(data, list(...)$theme, "label.size")
  }
)

.bb_effect_annotation_layer <- function(data, hjust, vjust) {
  geom <- if (.bb_theme_sets_geom_text()) {
    .bb_GeomTicklabelText
  } else {
    ggplot2::GeomText
  }
  ggplot2::layer(
    geom = geom,
    stat = "identity",
    data = data,
    mapping = ggplot2::aes(label = .data$effect_annotation),
    position = "identity",
    show.legend = FALSE,
    inherit.aes = FALSE,
    params = list(
      x = Inf,
      y = Inf,
      hjust = hjust,
      vjust = vjust,
      parse = FALSE,
      check_overlap = FALSE,
      na.rm = FALSE
    )
  )
}

# ggpubr always passes label.size = 3.88 mm. Drop that fixed value so
# the bracket geom can read the theme while the plot is built.
.bb_use_ticklabel_brackets <- function(plot) {
  if (!.bb_theme_sets_geom_text()) {
    return(plot)
  }
  n <- length(plot$layers)
  if (n < 1L) {
    return(plot)
  }
  layer <- plot$layers[[n]]
  if (!inherits(layer$geom, "GeomBracket")) {
    return(plot)
  }
  layer$aes_params$label.size <- NULL
  layer$geom <- .bb_GeomTicklabelBracket
  plot$layers[[n]] <- layer
  plot
}
