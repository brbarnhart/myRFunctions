# tests/testthat/test-bb_ggbarplot.R

library(testthat)
library(ggplot2)

setup_bar_data <- function() {
  set.seed(1)
  expand.grid(
    group = factor(c("A", "B")),
    condition = factor(c("Cont", "Stim")),
    rep = 1:8
  ) |>
    transform(
      score = 2 +
        as.numeric(condition) +
        as.numeric(group) +
        rnorm(length(group), sd = 0.5)
    )
}

test_that("bb_ggbarplot returns a ggplot with bars, points, and zero line", {
  skip_if_not_installed("ggpubr")
  skip_if_not_installed("ggplot2")

  dat <- setup_bar_data()
  p <- bb_ggbarplot(dat, x = "group", y = "score", fill = "condition")

  expect_s3_class(p, "ggplot")
  layer_classes <- vapply(p$layers, function(l) class(l$geom)[[1]], character(1))
  expect_true(any(grepl("Point", layer_classes)))
  expect_true(any(grepl("Hline", layer_classes)))
})

test_that("bb_ggbarplot can suppress points and hline", {
  skip_if_not_installed("ggpubr")
  skip_if_not_installed("ggplot2")

  dat <- setup_bar_data()
  p <- bb_ggbarplot(
    dat,
    x = "group",
    y = "score",
    fill = "condition",
    add_points = FALSE,
    add_hline = FALSE
  )

  layer_classes <- vapply(p$layers, function(l) class(l$geom)[[1]], character(1))
  expect_false(any(grepl("Point", layer_classes)))
  expect_false(any(grepl("Hline", layer_classes)))
})
