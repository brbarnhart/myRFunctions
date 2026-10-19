# tests/testthat/test-bbmake_pairwise_plot.R

library(testthat)
library(emmeans)
library(lme4)
library(lmerTest)
library(dplyr)
library(tibble)
library(ggplot2)

setup_plot_data <- function() {
  set.seed(42)
  dat <- expand.grid(
    subj  = factor(1:12),
    Stim  = factor(c("Low", "Med", "High"), levels = c("Low", "Med", "High")),
    Group = factor(c("Control", "Treatment"))
  ) |>
    as_tibble() |>
    mutate(
      y = 10 +
        3 * (Stim == "Med") +
        6 * (Stim == "High") +
        4 * (Group == "Treatment") +
        2.5 * (Stim == "High" & Group == "Treatment") +
        rnorm(n(), sd = 2.5)
    )

  mod <- suppressWarnings(lmer(y ~ Stim * Group + (1 | subj), data = dat))
  emm <- emmeans(mod, ~ Stim | Group)
  pw  <- pairs(emm, adjust = "tukey")
  list(mod = mod, emm = emm, pw = pw)
}

test_that("bbmake_pairwise_plot returns a ggplot with auto pairs + facets", {
  skip_if_not_installed("lmerTest")
  skip_if_not_installed("emmeans")
  skip_if_not_installed("ggplot2")
  skip_if_not_installed("ggpubr")

  s <- setup_plot_data()
  expect_message(
    p <- bbmake_pairwise_plot(s$emm),
    "Holm across 6 comparisons"
  )

  expect_s3_class(p, "ggplot")
  expect_true(inherits(p, "gg"))

  layer_classes <- vapply(p$layers, function(l) class(l$geom)[[1]], character(1))
  expect_true(any(grepl("Errorbar", layer_classes)))
  expect_true(any(grepl("Point", layer_classes)))
  # Facet by Group should be present
  expect_true(!is.null(p$facet))
  expect_false(inherits(p$facet, "FacetNull"))
})

test_that("bbmake_pairwise_plot cross.adjust is applied without a prebuilt pw", {
  skip_if_not_installed("lmerTest")
  skip_if_not_installed("emmeans")
  skip_if_not_installed("ggpubr")

  s <- setup_plot_data()
  expect_no_error({
    p <- bbmake_pairwise_plot(
      s$emm,
      adjust = "none",
      cross.adjust = "holm"
    )
  })
  expect_s3_class(p, "ggplot")

  none <- bbmake_pairwise_table(
    pairs(s$emm, reverse = TRUE, adjust = "none"),
    cross.adjust = "none"
  )
  holm <- bbmake_pairwise_table(
    pairs(s$emm, reverse = TRUE, adjust = "none"),
    cross.adjust = "holm"
  )
  expect_true(all(holm$p.value >= none$p.value - 1e-12))
  expect_true(any(holm$p.value > none$p.value + 1e-12))
})

test_that("bbmake_pairwise_plot accepts explicit pw and x as string", {
  skip_if_not_installed("lmerTest")
  skip_if_not_installed("emmeans")
  skip_if_not_installed("ggpubr")

  s <- setup_plot_data()
  expect_no_error({
    p <- bbmake_pairwise_plot(
      s$emm,
      x = "Stim",
      pw = s$pw,
      model = s$mod,
      y.adjust = 0.5
    )
  })
  expect_s3_class(p, "ggplot")
})

test_that("bbmake_pairwise_plot accepts bare-name x and custom facets", {
  skip_if_not_installed("lmerTest")
  skip_if_not_installed("emmeans")
  skip_if_not_installed("ggpubr")

  s <- setup_plot_data()
  p <- bbmake_pairwise_plot(
    s$emm,
    x = Stim,
    facets = ~ Group,
    pw = s$pw,
    connect = TRUE,
    show_zero_line = FALSE,
    annotate_effect = FALSE
  )
  expect_s3_class(p, "ggplot")

  layer_classes <- vapply(p$layers, function(l) class(l$geom)[[1]], character(1))
  expect_true(any(grepl("Line", layer_classes)))
})

test_that("bbmake_pairwise_plot works with prebuilt pw_table", {
  skip_if_not_installed("lmerTest")
  skip_if_not_installed("emmeans")
  skip_if_not_installed("ggpubr")

  s <- setup_plot_data()
  tab <- bbmake_pairwise_table(s$pw, model = s$mod)

  p <- bbmake_pairwise_plot(s$emm, pw_table = tab, model = s$mod)
  expect_s3_class(p, "ggplot")
})

annotation_text_size <- function(p) {
  built <- ggplot2::ggplot_build(p)
  sizes <- vapply(built$data, function(d) {
    if ("label" %in% names(d) && "size" %in% names(d) && any(is.infinite(d$x))) {
      as.numeric(d$size[[1]])
    } else {
      NA_real_
    }
  }, numeric(1))
  sizes <- sizes[!is.na(sizes)]
  if (length(sizes) == 0L) {
    return(NA_real_)
  }
  sizes[[1]]
}

axis_text_size_mm <- function(p) {
  built <- ggplot2::ggplot_build(p)
  ggplot2::calc_element("axis.text", built$plot@theme)$size / ggplot2::.pt
}

bracket_label_size <- function(p) {
  built <- ggplot2::ggplot_build(p)
  sizes <- vapply(built$data, function(d) {
    if ("label.size" %in% names(d) && "label" %in% names(d)) {
      as.numeric(d$label.size[[1]])
    } else {
      NA_real_
    }
  }, numeric(1))
  sizes <- sizes[!is.na(sizes)]
  if (length(sizes) == 0L) {
    return(NA_real_)
  }
  sizes[[1]]
}

test_that("effect annotation and stars match axis tick labels", {
  skip_if_not_installed("lmerTest")
  skip_if_not_installed("emmeans")
  skip_if_not_installed("ggpubr")
  skip_if(
    utils::packageVersion("ggplot2") < "4.0.0",
    "tick-label geom text size requires ggplot2 >= 4.0.0"
  )

  s <- setup_plot_data()
  p <- bbmake_pairwise_plot(s$emm, pw = s$pw, model = s$mod)

  small <- p + ggplot2::theme_grey(base_size = 11)
  large <- p + ggplot2::theme_grey(base_size = 22)
  ticks <- p +
    ggplot2::theme_grey(base_size = 16) +
    ggplot2::theme(axis.text = ggplot2::element_text(size = 9))

  size_small <- annotation_text_size(small)
  size_large <- annotation_text_size(large)

  expect_false(is.na(size_small))
  expect_equal(size_small, axis_text_size_mm(small), tolerance = 1e-6)
  expect_equal(size_large, axis_text_size_mm(large), tolerance = 1e-6)
  expect_equal(size_large, size_small * 2, tolerance = 1e-6)
  expect_equal(bracket_label_size(small), size_small, tolerance = 1e-6)
  expect_equal(bracket_label_size(large), size_large, tolerance = 1e-6)
  expect_equal(annotation_text_size(ticks) * ggplot2::.pt, 9, tolerance = 1e-6)
  expect_equal(bracket_label_size(ticks) * ggplot2::.pt, 9, tolerance = 1e-6)
})

test_that("bbmake_pairwise_plot y_expand defaults and overrides", {
  skip_if_not_installed("lmerTest")
  skip_if_not_installed("emmeans")
  skip_if_not_installed("ggpubr")

  s <- setup_plot_data()
  p_ann <- bbmake_pairwise_plot(s$emm, pw = s$pw, annotate_effect = TRUE)
  p_off <- bbmake_pairwise_plot(s$emm, pw = s$pw, annotate_effect = FALSE)
  p_vec <- bbmake_pairwise_plot(
    s$emm, pw = s$pw, annotate_effect = TRUE, y_expand = c(0.05, 0.4)
  )
  p_exp <- bbmake_pairwise_plot(
    s$emm,
    pw = s$pw,
    y_expand = ggplot2::expansion(mult = c(0.01, 0.5), add = c(0, 1))
  )

  expect_equal(
    p_ann$scales$get_scales("y")$expand,
    ggplot2::expansion(mult = c(0.02, 0.28))
  )
  expect_equal(
    p_off$scales$get_scales("y")$expand,
    ggplot2::expansion(mult = c(0.02, 0.18))
  )
  expect_equal(
    p_vec$scales$get_scales("y")$expand,
    ggplot2::expansion(mult = c(0.05, 0.4))
  )
  expect_equal(
    p_exp$scales$get_scales("y")$expand,
    ggplot2::expansion(mult = c(0.01, 0.5), add = c(0, 1))
  )
})

test_that("bbmake_pairwise_plot can suppress facets", {
  skip_if_not_installed("lmerTest")
  skip_if_not_installed("emmeans")
  skip_if_not_installed("ggpubr")

  s <- setup_plot_data()
  p <- bbmake_pairwise_plot(s$emm, facets = FALSE, pw = s$pw)
  expect_s3_class(p, "ggplot")
  expect_true(inherits(p$facet, "FacetNull"))
})

test_that("helpers compose into a custom ggplot like the user snippet", {
  skip_if_not_installed("lmerTest")
  skip_if_not_installed("emmeans")
  skip_if_not_installed("ggpubr")

  s <- setup_plot_data()
  df <- bb_emm_df(s$emm)

  p <- ggplot(df, aes(x = Stim, y = y, ymin = ymin, ymax = ymax)) +
    geom_errorbar(width = 0.2, linewidth = 0.9, color = "black") +
    geom_point(size = 2, color = "black") +
    facet_wrap(~ Group) +
    bb_add_pairwise(s$pw) +
    labs(x = "Stimulation", y = "Response")

  expect_s3_class(p, "ggplot")
  built <- ggplot_build(p)
  br <- built$data[[length(built$data)]]
  expect_true("annotation" %in% names(br) || "label" %in% names(br))
})
