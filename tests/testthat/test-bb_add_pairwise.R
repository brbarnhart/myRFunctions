# tests/testthat/test-bb_add_pairwise.R

library(testthat)
library(emmeans)
library(dplyr)
library(tibble)
library(ggplot2)

setup_dodge_plot <- function() {
  set.seed(42)
  dat <- expand.grid(
    subj = factor(1:8),
    Stim = factor(c("Cont", "Opto")),
    Sex = factor(c("Male", "Female")),
    Diet = factor(c("Chow", "HFHS")),
    Satiety = factor(c("Fed", "Fasted"))
  ) |>
    as_tibble() |>
    mutate(
      Breakpoint = 8 +
        6 * (Stim == "Opto" & Sex == "Female" & Diet == "Chow" & Satiety == "Fed") +
        rnorm(n(), sd = 1.2)
    )

  mod <- lm(Breakpoint ~ Stim * Sex * Diet * Satiety, data = dat)
  emm <- emmeans(mod, ~ Stim | Sex * Diet * Satiety)
  pw <- pairs(emm, reverse = TRUE)
  list(dat = dat, mod = mod, emm = emm, pw = pw)
}

base_dodge_plot <- function(dat, col_width = 0.8) {
  ggplot(dat, aes(x = interaction(Diet, Satiety), y = Breakpoint, fill = Stim)) +
    geom_bar(
      stat = "summary",
      fun = mean,
      width = col_width,
      position = position_dodge(width = col_width)
    ) +
    geom_point(
      position = position_jitterdodge(
        dodge.width = col_width,
        jitter.width = col_width / 3
      )
    ) +
    facet_wrap(~Sex)
}

bracket_data <- function(p) {
  built <- ggplot2::ggplot_build(p)
  built$data[[length(built$data)]]
}

test_that("bb_add_pairwise prints a ggplot_add placeholder", {
  obj <- bb_add_pairwise(data.frame(contrast = "A - B", p.value = 0.01))
  expect_s3_class(obj, "bb_pairwise_layer")
  expect_output(print(obj), "bb_add_pairwise")
})

test_that("bb_add_pairwise adds dodged brackets on interaction x + facet", {
  skip_if_not_installed("emmeans")
  skip_if_not_installed("ggplot2")
  skip_if_not_installed("ggpubr")

  s <- setup_dodge_plot()
  p <- base_dodge_plot(s$dat) + bb_add_pairwise(s$pw, dodge_width = 0.8)

  expect_s3_class(p, "ggplot")

  built <- ggplot_build(p)
  expect_gt(length(built$data), 2L)

  bars <- built$data[[1]]
  br <- built$data[[length(built$data)]]

  expect_true("annotation" %in% names(br) || "label" %in% names(br))
  expect_true(all(c("xmin", "xmax") %in% names(br)))

  # First interaction level is dodged to 0.8 and 1.2 with width 0.8, n = 2
  bar_x <- sort(unique(round(bars$x[bars$x < 1.5], 6)))
  expect_equal(bar_x, c(0.8, 1.2))

  # At least one bracket sits on that first x category
  br_pairs <- unique(round(cbind(br$xmin, br$xmax), 6))
  expect_true(any(br_pairs[, 1] == 0.8 & br_pairs[, 2] == 1.2))
})

test_that("bb_add_pairwise infers dodge_width from position_dodge", {
  skip_if_not_installed("emmeans")
  skip_if_not_installed("ggplot2")
  skip_if_not_installed("ggpubr")

  s <- setup_dodge_plot()
  p <- base_dodge_plot(s$dat, col_width = 0.8) + bb_add_pairwise(s$pw)

  dummy <- ggplot(s$dat, aes(x = interaction(Diet, Satiety), y = Breakpoint, fill = Stim)) +
    geom_bar(stat = "summary", fun = mean, position = position_dodge(width = 0.8))
  obj <- bb_add_pairwise(s$pw)
  sig <- myRFunctions:::.bb_pairwise_layer_data(obj, dummy)

  n <- n_distinct(s$dat$Stim)
  expect_equal(n, 2L)
  expect_equal(sort(unique(round(sig$xmax - sig$xmin, 6))), 0.4)
})

test_that("bb_add_pairwise hide.ns drops ns rows", {
  skip_if_not_installed("emmeans")
  skip_if_not_installed("ggplot2")
  skip_if_not_installed("ggpubr")

  s <- setup_dodge_plot()
  dummy <- base_dodge_plot(s$dat)
  sig_all <- myRFunctions:::.bb_pairwise_layer_data(
    bb_add_pairwise(s$pw, hide.ns = FALSE, dodge_width = 0.8),
    dummy
  )
  sig_sig <- myRFunctions:::.bb_pairwise_layer_data(
    bb_add_pairwise(s$pw, hide.ns = TRUE, dodge_width = 0.8),
    dummy
  )

  expect_equal(nrow(sig_all), 8L)
  expect_true(nrow(sig_sig) < nrow(sig_all))
  expect_true(all(sig_sig$p.signif != "ns"))
  expect_true("Sex" %in% names(sig_all))
  expect_true(all(c("xmin", "xmax", "y.position", "p.signif") %in% names(sig_all)))
})

test_that("bb_add_pairwise y.position sits above the cell max", {
  skip_if_not_installed("emmeans")
  skip_if_not_installed("ggplot2")
  skip_if_not_installed("ggpubr")

  s <- setup_dodge_plot()
  dummy <- base_dodge_plot(s$dat)
  sig <- myRFunctions:::.bb_pairwise_layer_data(
    bb_add_pairwise(s$pw, dodge_width = 0.8, y.adjust = 0.5),
    dummy
  )

  cell <- s$dat |>
    filter(Sex == "Female", Diet == "Chow", Satiety == "Fed")
  row <- sig |>
    filter(Sex == "Female", Diet == "Chow", Satiety == "Fed")

  expect_equal(nrow(row), 1L)
  expect_gt(row$y.position[[1]], max(cell$Breakpoint))
})

test_that("bb_add_pairwise group override matches fill variable", {
  skip_if_not_installed("emmeans")
  skip_if_not_installed("ggplot2")
  skip_if_not_installed("ggpubr")

  s <- setup_dodge_plot()
  dummy <- base_dodge_plot(s$dat)
  sig <- myRFunctions:::.bb_pairwise_layer_data(
    bb_add_pairwise(s$pw, group = "Stim", dodge_width = 0.8),
    dummy
  )
  expect_equal(nrow(sig), 8L)
  expect_true(all(sig$group1 %in% c("Cont", "Opto")))
  expect_true(all(sig$group2 %in% c("Cont", "Opto")))
})

test_that("bb_add_pairwise x-axis mode spans factor levels", {
  skip_if_not_installed("emmeans")
  skip_if_not_installed("ggplot2")
  skip_if_not_installed("ggpubr")

  set.seed(1)
  dat <- expand.grid(
    Stim = factor(c("Low", "Med", "High"), levels = c("Low", "Med", "High")),
    Group = factor(c("Control", "Treatment")),
    rep = 1:6
  ) |>
    as_tibble() |>
    mutate(y = 10 + as.numeric(Stim) + as.numeric(Group) + rnorm(n()))

  mod <- lm(y ~ Stim * Group, data = dat)
  emm <- emmeans(mod, ~ Stim | Group)
  pw <- pairs(emm)

  p <- ggplot(dat, aes(x = Stim, y = y)) +
    stat_summary(fun = mean, geom = "point") +
    facet_wrap(~Group) +
    bb_add_pairwise(pw)

  expect_s3_class(p, "ggplot")

  dummy <- ggplot(dat, aes(x = Stim, y = y)) + facet_wrap(~Group)
  sig <- myRFunctions:::.bb_pairwise_layer_data(bb_add_pairwise(pw), dummy)

  expect_true(all(sig$xmin %in% c("Low", "Med", "High")))
  expect_true(all(sig$xmax %in% c("Low", "Med", "High")))
  expect_true("Group" %in% names(sig))
  expect_gt(nrow(sig), 1L)
})

test_that("bb_add_pairwise errors when contrast groups cannot be mapped", {
  skip_if_not_installed("emmeans")
  skip_if_not_installed("ggplot2")

  s <- setup_dodge_plot()
  p <- ggplot(s$dat, aes(x = Diet, y = Breakpoint)) +
    geom_point() +
    facet_wrap(~Sex)

  expect_error(p + bb_add_pairwise(s$pw), "fill/colour/group")
})

test_that("bb_add_pairwise uses ymax when CI columns are present", {
  skip_if_not_installed("emmeans")
  skip_if_not_installed("ggplot2")
  skip_if_not_installed("ggpubr")

  set.seed(1)
  dat <- expand.grid(
    Stim = factor(c("Low", "Med", "High"), levels = c("Low", "Med", "High")),
    Group = factor(c("Control", "Treatment")),
    rep = 1:6
  ) |>
    as_tibble() |>
    mutate(y = 10 + as.numeric(Stim) + as.numeric(Group) + rnorm(n()))

  mod <- lm(y ~ Stim * Group, data = dat)
  emm <- emmeans(mod, ~ Stim | Group)
  pw <- pairs(emm)
  df <- bb_emm_df(emm)

  dummy <- ggplot(df, aes(x = Stim, y = y)) +
    geom_point() +
    facet_wrap(~Group)
  sig <- myRFunctions:::.bb_pairwise_layer_data(
    bb_add_pairwise(pw, y.adjust = 0, step = 0),
    dummy
  )

  for (g in unique(as.character(df$Group))) {
    expected <- max(df$ymax[df$Group == g], na.rm = TRUE)
    got <- unique(sig$y.position[as.character(sig$Group) == g])
    expect_equal(got, expected)
  }
})
