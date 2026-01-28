# tests/testthat/test-bbmake_anova_table.R

library(testthat)
library(lme4)
library(lmerTest)     # for p-values and DenDF in anova(lmer)
library(effectsize)
library(dplyr)
library(tibble)

# ── Helper: create test data and models ─────────────────────────────────────
setup_models <- function() {
  set.seed(42)

  # Balanced small dataset
  dat <- expand.grid(
    ID   = factor(1:16),
    Sex  = factor(c("M", "F")),
    Diet = factor(c("Ctrl", "HF")),
    Stim = factor(c("Low", "High"))
  ) |>
    as_tibble() |>
    mutate(
      Total_Correct_Pokes = 20 +
        4 * (Sex == "M") +
        7 * (Diet == "HF") +
        5.5 * (Stim == "High") +
        3 * (Sex == "M" & Diet == "HF") +
        rnorm(n(), sd = 3.2)
    )

  # Fixed-effects only model
  mod_lm <- lm(Total_Correct_Pokes ~ Sex * Diet * Stim, data = dat)

  # Mixed model (with random intercept)
  mod_lmer <- lmer(Total_Correct_Pokes ~ Sex * Diet * Stim + (1 | ID), data = dat)

  list(
    lm    = mod_lm,
    lmer  = mod_lmer,
    data  = dat
  )
}


test_that("bbmake_anova_table works on lm() model and returns expected structure", {
  mods <- setup_models()
  model <- mods$lm

  tab <- bbmake_anova_table(model)

  # Basic class & structure
  expect_s3_class(tab, "data.frame")
  expect_s3_class(tab, "tbl_df")           # since you use tibble
  expect_true(nrow(tab) >= 1)

  # Required columns (exact names as in your function)
  expected_cols <- c(
    "Effect", "NumDF", "DenDF", "F value", "Pr(>F)",
    "Omega2 (partial)", "Omega2 (95% CI)"
  )
  expect_true(all(expected_cols %in% names(tab)))

  # For lm(): DenDF should be NA
  expect_true(all(is.na(tab$DenDF)))

  # Omega columns should exist and be character/numeric as expected
  expect_true(is.numeric(tab$`Omega2 (partial)`))
  expect_true(is.character(tab$`Omega2 (95% CI)`))

  # At least some effects should have omega values (not all NA)
  expect_true(any(!is.na(tab$`Omega2 (partial)`)))

  # p-values should be formatted (your function uses format.pval)
  expect_true(any(grepl("^<", tab$`Pr(>F)`) | grepl("^[0-9]", tab$`Pr(>F)`)))
})


test_that("bbmake_anova_table works on lmer() model and handles DenDF correctly", {
  skip_if_not_installed("lmerTest")   # safety net

  mods <- setup_models()
  model <- mods$lmer

  tab <- bbmake_anova_table(model)

  # Basic structure same as above
  expect_s3_class(tab, "data.frame")
  expect_s3_class(tab, "tbl_df")
  expect_true(nrow(tab) >= 1)

  expected_cols <- c(
    "Effect", "NumDF", "DenDF", "F value", "Pr(>F)",
    "Omega2 (partial)", "Omega2 (95% CI)"
  )
  expect_true(all(expected_cols %in% names(tab)))

  # For lmer (with lmerTest): DenDF should be present & numeric
  expect_true(is.numeric(tab$DenDF))
  expect_true(all(!is.na(tab$DenDF)))         # should have Satterthwaite approx

  # Reasonable values
  expect_true(all(tab$NumDF > 0))
  expect_true(all(tab$`F value` > 0 | is.na(tab$`F value`)))

  # Omega should be present for main effects / interactions
  expect_true(any(tab$`Omega2 (partial)` > 0, na.rm = TRUE))

  # No error when printing/summarizing
  expect_no_error({
    print(tab)
    summary(tab)
  })
})


# Optional extra test: graceful behavior on singular fit / tiny random variance
test_that("bbmake_anova_table does not crash on boundary/singular lmer fit", {
  skip_if_not_installed("lmerTest")

  mods <- setup_models()
  dat <- mods$data

  # Force very weak clustering → variance → ~0
  dat$y_weak <- dat$Total_Correct_Pokes + rnorm(nrow(dat), sd = 0.01)

  mod_singular <- suppressWarnings(
    lmer(y_weak ~ Sex + Diet + Stim + (1 | ID), data = dat)
  )

  # Should still run (even with warning)
  expect_no_error({
    tab <- bbmake_anova_table(mod_singular)
  })

  expect_true(is.data.frame(tab))
  expect_true(nrow(tab) > 0)
})
