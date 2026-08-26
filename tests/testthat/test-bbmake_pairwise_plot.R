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
  p <- bbmake_pairwise_plot(s$emm, model = s$mod)

  expect_s3_class(p, "ggplot")
  expect_true(inherits(p, "gg"))

  layer_classes <- vapply(p$layers, function(l) class(l$geom)[[1]], character(1))
  expect_true(any(grepl("Errorbar", layer_classes)))
  expect_true(any(grepl("Point", layer_classes)))
  # Facet by Group should be present
  expect_true(!is.null(p$facet))
  expect_false(inherits(p$facet, "FacetNull"))
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
  df  <- bb_emm_df(s$emm)
  sig <- bb_pairwise_labels(s$pw, emm = s$emm, model = s$mod)

  p <- ggplot(df, aes(x = Stim, y = y, ymin = ymin, ymax = ymax)) +
    geom_errorbar(width = 0.2, linewidth = 0.9, color = "black") +
    geom_point(size = 2, color = "black") +
    facet_wrap(~ Group) +
    ggpubr::stat_pvalue_manual(
      sig,
      label = "p.signif",
      xmin = "group1",
      xmax = "group2",
      y.position = "y.position",
      tip.length = 0.01
    ) +
    labs(x = "Stimulation", y = "Response")

  expect_s3_class(p, "ggplot")
})
