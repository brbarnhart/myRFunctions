# tests/testthat/test-bbmake_pairwise_plot.R

library(testthat)
library(emmeans)
library(lme4)
library(lmerTest)     # for anova() with p-values
library(tidyverse)
library(ggpubr)       # for stat_pvalue_manual

# We'll assume the function is already available in the environment
# (e.g. loaded via devtools::load_all() or source())

setup_data <- function() {
  # ── Minimal reproducible dataset ───────────────────────────────────────
  set.seed(42)
  dat <- expand.grid(
    subj   = factor(1:12),
    Stim   = factor(c("Low", "Med", "High")),
    Group  = factor(c("Control", "Treatment"))
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

  mod <- lmer(y ~ Stim * Group + (1 | subj), data = dat)
  emm <- emmeans(mod, ~ Stim | Group)
  pw  <- pairs(emm, adjust = "tukey")

  pw_table <- summary(pw) |>
    as_tibble() |>
    select(contrast, estimate, SE, df, t.ratio, p.value) |>
    mutate(
      `Cohen's d`    = abs(estimate / SE) * 0.8,
      lower.CL       = NA_real_,
      upper.CL       = NA_real_
    )

  list(emm = emm, pw_table = pw_table)
}


test_that("bbmake_pairwise_plot returns a ggplot object", {
  dat <- setup_data()
  emm      <- dat$emm
  pw_table <- dat$pw_table

  pw  <- pairs(emm, adjust = "tukey")

  # Fake pw_table (mimics output of bbmake_pairwise_table)
  pw_table <- summary(pw) |>
    as_tibble() |>
    select(contrast, estimate, SE, df, t.ratio, p.value) |>
    mutate(
      `Cohen's d`    = abs(estimate / SE) * 0.8,   # rough fake
      lower.CL       = NA_real_,
      upper.CL       = NA_real_
    )

  # ── Test basic execution ───────────────────────────────────────────────
  p <- bbmake_pairwise_plot(
    emm      = emm,
    pw_table = pw_table,
    formula  = ~ Stim | Group,
    # group1   = "Low",
    # group2   = "Med",
    y.adjust = 3
  )

  expect_s3_class(p, "ggplot")
  expect_true(inherits(p, "gg"))
  expect_true(length(p$layers) >= 4)          # at least hline + line + point + errorbar

  # Optional: check that key geoms are present (loose check)
  layer_classes <- sapply(p$layers, function(l) class(l$geom)[1])
  expect_true(any(grepl("GeomErrorbar", layer_classes)))
  expect_true(any(grepl("GeomLine",     layer_classes)))

  # ── Test it doesn't error with different inputs ────────────────────────
  expect_no_error(
    bbmake_pairwise_plot(
      emm      = emm,
      pw_table = pw_table,
      formula  = ~ Stim | Group,
      # group1   = "Med",
      # group2   = "High",
      y.adjust = 5
    )
  )

  # You can add more specific expectations later, e.g.
  # expect_equal(p$labels$y, "emmean")   # or whatever your yvar is called
})


# Optional second test block for multiple comparisons / vector y.adjust
test_that("bbmake_pairwise_plot handles multiple contrasts", {
  dat <- setup_data()
  emm      <- dat$emm
  pw_table <- dat$pw_table

  pw_table_multi <- pw_table |>
    slice(1:2)   # pretend we have two rows

  expect_no_error(
    bbmake_pairwise_plot(
      emm      = emm,
      pw_table = pw_table_multi,
      formula  = ~ Stim | Group,
      # group1   = c("Low", "Med"),
      # group2   = c("Med", "High"),
      y.adjust = c(2, 5)   # test vector input (note: your current code uses only first value)
    )
  )
})
