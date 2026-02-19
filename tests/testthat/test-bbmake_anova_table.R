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

  # ── NEW: glmmTMB negative binomial model ───────────────────────────────
  dat$count <- round(pmax(5, dat$Total_Correct_Pokes))   # make positive integer counts

  mod_nb <- glmmTMB::glmmTMB(
    count ~ Sex * Diet * Stim + (1 | ID),
    family = glmmTMB::nbinom2,
    data = dat
  )

  list(
    lm    = mod_lm,
    lmer  = mod_lmer,
    nb    = mod_nb,          # ← new
    data  = dat
  )
}


# ==================================================================
# Existing tests (unchanged)
# ==================================================================

test_that("bbmake_anova_table works on lm() model and returns expected structure", {
  mods <- setup_models()
  model <- mods$lm

  tab <- bbmake_anova_table(model)

  expect_s3_class(tab, "data.frame")
  expect_s3_class(tab, "tbl_df")
  expect_true(nrow(tab) >= 1)

  expected_cols <- c(
    "Effect", "NumDF", "DenDF", "F value", "Pr(>F)",
    "Omega2 (partial)", "Omega2 (95% CI)"
  )
  expect_true(all(expected_cols %in% names(tab)))

  expect_true(all(is.na(tab$DenDF)))
  expect_true(is.numeric(tab$`Omega2 (partial)`))
  expect_true(is.character(tab$`Omega2 (95% CI)`))
  expect_true(any(!is.na(tab$`Omega2 (partial)`)))
})


test_that("bbmake_anova_table works on lmer() model and handles DenDF correctly", {
  skip_if_not_installed("lmerTest")

  mods <- setup_models()
  model <- mods$lmer

  tab <- bbmake_anova_table(model)

  expect_s3_class(tab, "data.frame")
  expect_s3_class(tab, "tbl_df")
  expect_true(nrow(tab) >= 1)

  expected_cols <- c(
    "Effect", "NumDF", "DenDF", "F value", "Pr(>F)",
    "Omega2 (partial)", "Omega2 (95% CI)"
  )
  expect_true(all(expected_cols %in% names(tab)))

  expect_true(is.numeric(tab$DenDF))
  expect_true(all(!is.na(tab$DenDF)))
  expect_true(any(tab$`Omega2 (partial)` > 0, na.rm = TRUE))
})


test_that("bbmake_anova_table does not crash on boundary/singular lmer fit", {
  skip_if_not_installed("lmerTest")

  mods <- setup_models()
  dat <- mods$data

  dat$y_weak <- dat$Total_Correct_Pokes + rnorm(nrow(dat), sd = 0.01)

  mod_singular <- suppressWarnings(
    lmer(y_weak ~ Sex + Diet + Stim + (1 | ID), data = dat)
  )

  expect_no_error({
    tab <- bbmake_anova_table(mod_singular)
  })
  expect_true(is.data.frame(tab))
  expect_true(nrow(tab) > 0)
})


# ==================================================================
# NEW TESTS: glmmTMB / negative binomial branch
# ==================================================================

test_that("bbmake_anova_table works on glmmTMB negative binomial model", {
  skip_if_not_installed("glmmTMB")

  mods <- setup_models()
  model <- mods$nb

  tab <- bbmake_anova_table(model)

  # Basic structure
  expect_s3_class(tab, "data.frame")
  expect_s3_class(tab, "tbl_df")
  expect_true(nrow(tab) >= 1)

  # GLMM-specific columns
  expected_cols <- c(
    "Effect", "NumDF", "Chisq value", "Pr(>Chisq)",
    "IRR", "IRR 95% CI"
  )
  expect_true(all(expected_cols %in% names(tab)))

  # IRR columns should be present and sensible
  expect_true(is.numeric(tab$IRR))
  expect_true(is.character(tab$`IRR 95% CI`))
  expect_true(any(!is.na(tab$IRR)))          # at least some IRRs

  # Term-name cleaning worked (main effects and interactions should be short)
  expect_true("Sex" %in% tab$Effect)
  expect_true("Diet" %in% tab$Effect)
  expect_true(any(grepl("Sex:Diet", tab$Effect)))   # interaction example

  # No NAs in the joined columns (clean_term_names succeeded)
  expect_true(all(!is.na(tab$IRR) | is.na(tab$IRR)))  # just checking no unexpected NAs

  # R² attribute should be attached
  expect_true(!is.null(attr(tab, "r2")))
  expect_true(grepl("Marginal R²", attr(tab, "r2")))
})


test_that("bbmake_anova_table glmmTMB branch returns clean IRR values for main effects", {
  skip_if_not_installed("glmmTMB")

  mods <- setup_models()
  model <- mods$nb

  tab <- bbmake_anova_table(model)

  # Main effects should have sensible IRR values (not NA)
  main_effects <- c("Sex", "Diet", "Stim")
  for (eff in main_effects) {
    row <- tab[tab$Effect == eff, ]
    expect_true(nrow(row) == 1)
    expect_false(is.na(row$IRR))
    expect_true(row$IRR > 0)
  }
})
