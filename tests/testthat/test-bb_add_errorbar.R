# tests/testthat/test-bb_add_errorbar.R

library(testthat)
library(emmeans)
library(dplyr)
library(tibble)
library(ggplot2)

setup_errorbar <- function() {
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
  list(dat = dat, mod = mod, emm = emm)
}

base_dodge_plot <- function(dat, col_width = 0.8) {
  ggplot(dat, aes(x = interaction(Diet, Satiety), y = Breakpoint, fill = Stim)) +
    geom_bar(
      stat = "summary",
      fun = mean,
      width = col_width,
      position = position_dodge(width = col_width)
    ) +
    facet_wrap(~Sex)
}

errorbar_layer <- function(p) {
  built <- ggplot2::ggplot_build(p)
  hits <- vapply(built$data, function(d) {
    all(c("ymin", "ymax", "x") %in% names(d)) &&
      !all(d$ymin == 0 | is.na(d$ymin)) &&
      !"annotation" %in% names(d)
  }, logical(1))
  idx <- which(hits)
  if (length(idx) == 0L) {
    stop("No errorbar layer found")
  }
  built$data[[idx[[length(idx)]]]]
}

test_that("bb_add_errorbar prints a ggplot_add placeholder", {
  obj <- bb_add_errorbar(
    tibble(Stim = "A", emmean = 1, SE = 0.1, lower.CL = 0.8, upper.CL = 1.2)
  )
  expect_s3_class(obj, "bb_errorbar_layer")
  expect_output(print(obj), "bb_add_errorbar")
})

test_that("bb_add_errorbar dodges to match bar x positions", {
  skip_if_not_installed("emmeans")
  skip_if_not_installed("ggplot2")

  s <- setup_errorbar()
  p <- base_dodge_plot(s$dat) + bb_add_errorbar(s$emm, dodge_width = 0.8)

  expect_s3_class(p, "ggplot")
  built <- ggplot_build(p)
  bars <- built$data[[1]]
  eb <- errorbar_layer(p)

  bar_x <- sort(unique(round(bars$x, 6)))
  eb_x <- sort(unique(round(eb$x, 6)))
  expect_equal(eb_x, bar_x)
})

test_that("bb_add_errorbar CI limits match emmeans", {
  skip_if_not_installed("emmeans")
  skip_if_not_installed("ggplot2")

  s <- setup_errorbar()
  dummy <- base_dodge_plot(s$dat)
  df <- myRFunctions:::.bb_errorbar_layer_data(
    bb_add_errorbar(s$emm, interval = "ci"),
    dummy
  )
  emm_df <- bb_emm_df(s$emm)

  expect_equal(sort(df$ymin), sort(emm_df$ymin))
  expect_equal(sort(df$ymax), sort(emm_df$ymax))
  expect_equal(sort(df$y), sort(emm_df$y))
})

test_that("bb_add_errorbar SEM uses estimate ± SE", {
  skip_if_not_installed("emmeans")
  skip_if_not_installed("ggplot2")

  s <- setup_errorbar()
  dummy <- base_dodge_plot(s$dat)
  df <- myRFunctions:::.bb_errorbar_layer_data(
    bb_add_errorbar(s$emm, interval = "sem"),
    dummy
  )
  emm_df <- bb_emm_df(s$emm)

  expect_equal(df$ymin, df$y - df$SE)
  expect_equal(df$ymax, df$y + df$SE)
  expect_equal(sort(df$y), sort(emm_df$y))
  expect_false(isTRUE(all.equal(sort(df$ymin), sort(emm_df$ymin))))
})

test_that("backtransformed CIs are not ± z * SE on the response scale", {
  skip_if_not_installed("emmeans")
  skip_if_not_installed("ggplot2")

  set.seed(1)
  dat <- expand.grid(
    Stim = factor(c("Cont", "Opto")),
    Group = factor(c("A", "B")),
    rep = 1:12
  )
  dat$y <- rpois(nrow(dat), lambda = ifelse(dat$Stim == "Opto", 8, 3))
  mod <- glm(y ~ Stim * Group, data = dat, family = poisson)
  emm <- emmeans(mod, ~ Stim | Group, type = "response")

  dummy <- ggplot(dat, aes(x = Group, y = y, fill = Stim)) +
    geom_bar(stat = "summary", fun = mean, position = position_dodge(0.8))
  df_ci <- myRFunctions:::.bb_errorbar_layer_data(
    bb_add_errorbar(emm, interval = "ci"),
    dummy
  )
  df_sem <- myRFunctions:::.bb_errorbar_layer_data(
    bb_add_errorbar(emm, interval = "sem"),
    dummy
  )

  expect_true(attr(bb_emm_df(emm), "mean_col") %in% c("response", "rate"))
  expect_true(all(df_ci$ymin > 0))
  # Backtransformed CI is not the symmetric delta-method SEM band
  expect_false(isTRUE(all.equal(df_ci$ymin, df_sem$ymin)))
  expect_false(isTRUE(all.equal(
    df_ci$ymin,
    df_ci$y - 1.96 * df_ci$SE
  )))
})

test_that("bb_add_errorbar x-axis mode does not dodge", {
  skip_if_not_installed("emmeans")
  skip_if_not_installed("ggplot2")

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

  p <- ggplot(dat, aes(x = Stim, y = y)) +
    stat_summary(fun = mean, geom = "point") +
    facet_wrap(~Group) +
    bb_add_errorbar(emm)

  expect_s3_class(p, "ggplot")
  eb <- errorbar_layer(p)
  expect_equal(sort(unique(round(eb$x, 6))), c(1, 2, 3))
})

test_that("bb_add_errorbar errors without an x aesthetic", {
  skip_if_not_installed("emmeans")
  skip_if_not_installed("ggplot2")

  s <- setup_errorbar()
  p <- ggplot(s$dat, aes(y = Breakpoint)) + geom_point()
  expect_error(p + bb_add_errorbar(s$emm), "`x` aesthetic")
})

test_that("bbmake_pairwise_plot interval = 'sem' uses ± SE", {
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
  emm_df <- bb_emm_df(emm)

  p <- bbmake_pairwise_plot(emm, pw = pairs(emm), interval = "sem")
  built <- ggplot_build(p)
  eb <- built$data[[1]]

  expect_equal(sort(eb$ymin), sort(emm_df$y - emm_df$SE), tolerance = 1e-8)
  expect_equal(sort(eb$ymax), sort(emm_df$y + emm_df$SE), tolerance = 1e-8)
})
