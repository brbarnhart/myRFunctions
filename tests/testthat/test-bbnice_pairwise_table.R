# tests/testthat/test-bbnice_pairwise_table.R

library(testthat)
library(dplyr)
library(tibble)

needs_nice <- function() {
  skip_if_not_installed("flextable")
  skip_if_not_installed("rempsyc")
}

# Sensitivity-table shape: by, contrast, Model, statistic, effect, limits, p
sensitivity_like <- function() {
  tibble(
    Sex = factor(c("F", "F", "M", "M"), levels = c("F", "M")),
    contrast = c("B - A", "B - A", "B - A", "B - A"),
    Model = factor(
      c("Full", "Reduced", "Full", "Reduced"),
      levels = c("Full", "Reduced")
    ),
    z.ratio = c(1.23456, 2.2, 0.1, -3.456),
    df = c(12, 12, 11.5, 11.5),
    IRR = c(1.234, 1.2, 0.996, 2000.4),
    SE = c(0.111, 0.222, 0.333, 0.444),
    lower.CL = c(1.111, 0.996, 0.5, -1.234),
    upper.CL = c(2.222, 1.5, 1.204, 3),
    p.value = c(0.0004, 0.03, 0.2, 0.049)
  )
}

# Single-model Gaussian shape, including Cohen's d after p
gaussian_like <- function() {
  tibble(
    Diet = factor(c("Ctrl", "HF")),
    contrast = c("HF - Ctrl", "High - Low"),
    df = c(18.2, 18.2),
    t.ratio = c(2.345, -0.456),
    `Mean Difference` = c(4.567, -0.004),
    lower.CL = c(1.2, -1.5),
    upper.CL = c(7.8, 1.1),
    SE = c(1.5, 0.6),
    p.value = c(0.012, 0.4),
    d = c(0.8123, -0.0456)
  )
}

# merge_v stores vertical spans in body$spans$columns (rows is colspan)
col_spans <- function(ft, col) {
  j <- match(col, ft$col_keys)
  as.numeric(ft$body$spans$columns[, j])
}

ft_html <- function(ft) {
  obj <- flextable::htmltools_value(ft)
  html <- obj[[which(vapply(obj, is.character, logical(1)))[1]]]
  as.character(html)
}

# ==================================================================
# Column order and header names
# ==================================================================

test_that("default formatting keeps column order and renames headers", {
  needs_nice()
  tab <- sensitivity_like()
  ft <- bbnice_pairwise_table(tab)

  expect_s3_class(ft, "flextable")
  expect_equal(
    ft$col_keys,
    c("Sex", "Contrast", "Model", "z", "IRR", "95% CI", "p")
  )
  expect_equal(as.character(ft$body$dataset$Model), as.character(tab$Model))
  expect_equal(ft$body$dataset$IRR, tab$IRR)
  expect_equal(
    ft$body$dataset[["95% CI"]],
    c("[1.11, 2.22]", "[1.00, 1.50]", "[0.50, 1.20]", "[-1.23, 3.00]")
  )
})

test_that("drop = NULL keeps SE and df where the input put them", {
  needs_nice()
  ft <- bbnice_pairwise_table(sensitivity_like(), drop = NULL)
  expect_equal(
    ft$col_keys,
    c("Sex", "Contrast", "Model", "z", "df", "IRR", "SE", "95% CI", "p")
  )
})

test_that("gaussian tables keep Mean Difference, t, and Cohen's d in place", {
  needs_nice()
  ft <- bbnice_pairwise_table(gaussian_like())
  expect_equal(
    ft$col_keys,
    c("Diet", "Contrast", "t", "Mean Difference", "95% CI", "p", "d")
  )
  expect_equal(
    ft$body$dataset[["95% CI"]],
    c("[1.20, 7.80]", "[-1.50, 1.10]")
  )
})

test_that("a shuffled input is not forced back into a canonical order", {
  needs_nice()
  tab <- sensitivity_like() |>
    select(p.value, IRR, contrast, upper.CL, Sex, lower.CL, z.ratio)
  ft <- bbnice_pairwise_table(tab, drop = NULL, merge = FALSE)
  expect_equal(
    ft$col_keys,
    c("p", "IRR", "Contrast", "95% CI", "Sex", "z")
  )
})

test_that("rows stay in the input order", {
  needs_nice()
  tab <- sensitivity_like()[c(3, 1, 4, 2), ]
  ft <- bbnice_pairwise_table(tab, merge = FALSE)
  expect_equal(ft$body$dataset$IRR, tab$IRR)
  expect_equal(ft$body$dataset$p, tab$p.value)
})

# ==================================================================
# p values, intervals, and vertical merging
# ==================================================================

test_that("printed p values match rempsyc::format_p", {
  needs_nice()
  tab <- sensitivity_like()
  ft <- bbnice_pairwise_table(tab, stars = TRUE)
  html <- ft_html(ft)
  for (p in tab$p.value) {
    label <- rempsyc::format_p(p, stars = TRUE)
    escaped <- gsub("<", "&lt;", label, fixed = TRUE)
    expect_true(
      grepl(label, html, fixed = TRUE) || grepl(escaped, html, fixed = TRUE),
      info = label
    )
  }

  plain <- bbnice_pairwise_table(tab, stars = FALSE)
  plain_html <- ft_html(plain)
  bare <- rempsyc::format_p(0.03, stars = FALSE)
  starred <- rempsyc::format_p(0.03, stars = TRUE)
  expect_true(grepl(bare, plain_html, fixed = TRUE))
  expect_false(grepl(starred, plain_html, fixed = TRUE))
})

test_that("missing confidence limits print as a blank interval", {
  needs_nice()
  tab <- sensitivity_like()
  tab$lower.CL[2] <- NA_real_
  ft <- bbnice_pairwise_table(tab)
  expect_equal(ft$body$dataset[["95% CI"]][2], "")
})

test_that("by-variables are merged and contrast is left repeated", {
  needs_nice()
  tab <- sensitivity_like()
  ft <- bbnice_pairwise_table(tab)
  expect_equal(col_spans(ft, "Sex"), c(2, 0, 2, 0))
  expect_equal(col_spans(ft, "Contrast"), rep(1, 4))
  expect_equal(col_spans(ft, "Model"), rep(1, 4))

  both <- tibble(
    Sex = rep(c("F", "M"), each = 4),
    Diet = rep(rep(c("C", "H"), each = 2), 2),
    contrast = "B - A",
    z.ratio = 1,
    IRR = 1,
    lower.CL = 0.5,
    upper.CL = 1.5,
    p.value = 0.2
  )
  ft2 <- bbnice_pairwise_table(both)
  expect_equal(col_spans(ft2, "Sex"), c(4, 0, 0, 0, 4, 0, 0, 0))
  expect_equal(col_spans(ft2, "Diet"), rep(c(2, 0), 4))
})

test_that("merge accepts the input contrast name and can be turned off", {
  needs_nice()
  tab <- sensitivity_like()
  merged <- bbnice_pairwise_table(tab, merge = "contrast")
  expect_equal(col_spans(merged, "Contrast"), c(4, 0, 0, 0))
  expect_equal(col_spans(merged, "Sex"), rep(1, 4))

  plain <- bbnice_pairwise_table(tab, merge = FALSE)
  expect_equal(col_spans(plain, "Sex"), rep(1, 4))
})

test_that("stat headers are italic and the body uses theme_vanilla", {
  needs_nice()
  ft <- bbnice_pairwise_table(gaussian_like(), drop = NULL)
  html <- ft_html(ft)
  expect_match(html, "font-style:italic")
  expect_match(html, "font-weight:bold")
  expect_match(html, "Times New Roman")
})

# ==================================================================
# Pipes from the table builders, and input checks
# ==================================================================

test_that("sensitivity and pairwise builders pipe into the flextable", {
  needs_nice()
  skip_if_not_installed("emmeans")

  set.seed(1)
  dat <- expand.grid(
    Sex = factor(c("F", "M")),
    Stim = factor(c("A", "B")),
    Diet = factor(c("C", "H")),
    id = factor(1:6)
  )
  dat$y <- rpois(nrow(dat), lambda = 5)
  mods <- list(
    Full = glm(y ~ Sex * Stim * Diet, data = dat, family = poisson),
    Reduced = glm(
      y ~ Sex * Stim * Diet,
      data = dat[dat$id != "1", ],
      family = poisson
    )
  )

  sens <- bbmake_pairwise_sensitivity_table(
    mods,
    ~ Stim | Diet * Sex,
    by = c("Sex", "Diet"),
    adjust = "none",
    cross.adjust = "holm"
  )
  ft <- bbnice_pairwise_table(sens)
  expect_equal(
    ft$col_keys,
    c("Sex", "Diet", "Contrast", "Model", "z", "IRR", "95% CI", "p")
  )
  expect_equal(col_spans(ft, "Sex")[1], 4)

  emm <- emmeans::emmeans(mods$Full, ~ Stim | Sex, type = "response")
  pw <- pairs(emm, reverse = TRUE, adjust = "none")
  one <- bbmake_pairwise_table(pw, cross.adjust = "holm")
  ft1 <- bbnice_pairwise_table(one)
  expect_equal(ft1$col_keys[1:3], c("Sex", "Contrast", "z"))
  expect_equal(ft1$col_keys[length(ft1$col_keys)], "p")
})

test_that("invalid inputs error clearly", {
  needs_nice()
  tab <- sensitivity_like()
  expect_error(bbnice_pairwise_table(1), "data frame")
  expect_error(bbnice_pairwise_table(mtcars), "contrast")
  expect_error(bbnice_pairwise_table(tab, drop = 1), "`drop`")
  expect_error(bbnice_pairwise_table(tab, stars = "yes"), "`stars`")
  expect_error(bbnice_pairwise_table(tab, merge = 1), "`merge`")
  expect_error(bbnice_pairwise_table(tab, merge = "Group"), "not found")
  expect_error(
    bbnice_pairwise_table(tab[, names(tab) != "p.value"]),
    "p.value"
  )
})
