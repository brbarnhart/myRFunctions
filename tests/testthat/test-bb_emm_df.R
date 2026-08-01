# tests/testthat/test-bb_emm_df.R

library(testthat)
library(emmeans)
library(lme4)
library(lmerTest)
library(dplyr)
library(tibble)

setup_emm <- function() {
  set.seed(42)
  dat <- expand.grid(
    subj  = factor(1:10),
    Stim  = factor(c("Low", "Med", "High"), levels = c("Low", "Med", "High")),
    Group = factor(c("Control", "Treatment"))
  ) |>
    as_tibble() |>
    mutate(
      y = 10 +
        3 * (Stim == "Med") +
        6 * (Stim == "High") +
        4 * (Group == "Treatment") +
        rnorm(n(), sd = 2)
    )

  mod <- suppressWarnings(lmer(y ~ Stim * Group + (1 | subj), data = dat))
  emm <- emmeans(mod, ~ Stim | Group)
  list(mod = mod, emm = emm, dat = dat)
}

test_that("bb_emm_df normalizes emmean + lower.CL/upper.CL columns", {
  skip_if_not_installed("lmerTest")
  skip_if_not_installed("emmeans")

  s <- setup_emm()
  df <- bb_emm_df(s$emm)

  expect_s3_class(df, "tbl_df")
  expect_true(all(c("y", "ymin", "ymax") %in% names(df)))
  expect_true(all(c("Stim", "Group") %in% names(df)))
  expect_true(all(is.finite(df$y)))
  expect_true(all(df$ymin <= df$y, na.rm = TRUE))
  expect_true(all(df$ymax >= df$y, na.rm = TRUE))

  expect_equal(attr(df, "mean_col"), "emmean")
  expect_true(attr(df, "lower_col") %in% c("lower.CL", "asymp.LCL"))
  expect_equal(attr(df, "pri.vars"), "Stim")
  expect_equal(attr(df, "by.vars"), "Group")
})

test_that("bb_emm_df works on a plain data frame with response + asymp CIs", {
  raw <- tibble(
    Stim = factor(c("Low", "High")),
    response = c(1.2, 2.4),
    asymp.LCL = c(0.8, 1.9),
    asymp.UCL = c(1.6, 2.9)
  )

  df <- bb_emm_df(raw)
  expect_equal(df$y, raw$response)
  expect_equal(df$ymin, raw$asymp.LCL)
  expect_equal(df$ymax, raw$asymp.UCL)
  expect_equal(attr(df, "mean_col"), "response")
})

test_that("bb_emm_df errors when mean column is missing", {
  expect_error(bb_emm_df(tibble(x = 1, lower.CL = 0, upper.CL = 2)), "mean column")
})
