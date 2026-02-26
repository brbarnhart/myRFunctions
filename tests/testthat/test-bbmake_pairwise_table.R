# tests/testthat/test-bbmake_pairwise_table.R

library(testthat)
library(lme4)
library(lmerTest)
library(glmmTMB)
library(emmeans)
library(dplyr)
library(tibble)

# ── Helper: create test data and models (warnings suppressed) ─────────────────
setup_pairwise_models <- function() {
  set.seed(42)

  dat <- expand.grid(
    ID      = factor(1:20),
    Cohort  = factor(1:2),
    Sex     = factor(c("M", "F")),
    Diet    = factor(c("Ctrl", "HF")),
    Stim    = factor(c("Low", "High")),
    Satiety = factor(c("Low", "High"))
  ) |>
    as_tibble() |>
    mutate(
      Total_Correct_Pokes = 15 +
        5 * (Sex == "M") +
        8 * (Diet == "HF") +
        6 * (Stim == "High") +
        4 * (Satiety == "High") +
        rnorm(n(), sd = 4)
    )

  dat$count <- round(pmax(1, dat$Total_Correct_Pokes))

  # Suppress the common Hessian warning from small glmmTMB data
  mod_lmer <- suppressWarnings(
    lmer(Total_Correct_Pokes ~ Sex * Diet * Stim * Satiety + (1 | Cohort/ID),
         data = dat)
  )

  mod_nb <- suppressWarnings(
    glmmTMB(count ~ Sex * Diet * Stim * Satiety + (1 | Cohort/ID),
            family = nbinom2, data = dat)
  )

  list(lmer = mod_lmer, nb = mod_nb, data = dat)
}

# ==================================================================
# Gaussian models
# ==================================================================

test_that("bbmake_pairwise_table works on Gaussian (lmer) model - Mean Difference + Cohen's d", {
  skip_if_not_installed("lmerTest")
  skip_if_not_installed("emmeans")

  mods <- setup_pairwise_models()
  emm <- emmeans(mods$lmer, ~ Sex | Diet, infer = c(TRUE, TRUE))
  pw  <- pairs(emm, reverse = TRUE)

  tab <- bbmake_pairwise_table(pw, model = mods$lmer)

  # Column presence
  required_cols <- c("contrast", "Mean Difference", "SE", "Cohen's d", "d 95% CI")
  missing_cols  <- setdiff(required_cols, names(tab))
  expect_true(length(missing_cols) == 0,
              info = paste("Missing columns:", paste(missing_cols, collapse = ", ")))

  # Type checks with helpful messages
  expect_true(is.numeric(tab$`Mean Difference`),
              info = "Mean Difference column should be numeric")
  expect_true(is.character(tab$`d 95% CI`),
              info = "d 95% CI column should be character strings")

  # CI string format
  bad_ci <- tab$`d 95% CI`[!grepl("^\\[", tab$`d 95% CI`)]
  expect_true(length(bad_ci) == 0,
              info = paste("Bad CI strings (should start with '['):",
                           paste(head(bad_ci, 5), collapse = "; ")))

  # No NA p-values
  na_p <- sum(is.na(tab$p.value))
  expect_true(na_p == 0,
              info = paste(na_p, "NA p.value(s) found"))
})

# ==================================================================
# glmmTMB count models - type = "response"
# ==================================================================

test_that("bbmake_pairwise_table works on glmmTMB nbinom2 with type = 'response' (Rate Ratio branch)", {
  skip_if_not_installed("glmmTMB")
  skip_if_not_installed("emmeans")

  mods <- setup_pairwise_models()
  emm <- emmeans(mods$nb, ~ Sex | Diet, type = "response")
  pw  <- pairs(emm, reverse = TRUE, infer = TRUE)

  tab <- bbmake_pairwise_table(pw)

  # Column presence
  required_cols <- c("contrast", "Rate Ratio", "RR 95% CI", "% Change")
  missing_cols  <- setdiff(required_cols, names(tab))
  expect_true(length(missing_cols) == 0,
              info = paste("Missing columns:", paste(missing_cols, collapse = ", ")))

  # Type checks with helpful messages
  expect_true(is.numeric(tab$`Rate Ratio`),
              info = "Rate Ratio column should be numeric")
  expect_true(is.character(tab$`RR 95% CI`),
              info = "RR 95% CI column should be character strings")

  # CI string format
  bad_ci <- tab$`RR 95% CI`[!grepl("^\\[", tab$`RR 95% CI`)]
  expect_true(length(bad_ci) == 0,
              info = paste("Bad RR 95% CI strings:", paste(head(bad_ci, 5), collapse = "; ")))

  # % Change format
  bad_pct <- tab$`% Change`[!grepl("%$", tab$`% Change`)]
  expect_true(length(bad_pct) == 0,
              info = paste("Bad % Change strings:", paste(head(bad_pct, 5), collapse = "; ")))
})

test_that("bbmake_pairwise_table works on glmmTMB nbinom2 with type = 'response' + by grouping", {
  skip_if_not_installed("glmmTMB")
  skip_if_not_installed("emmeans")

  mods <- setup_pairwise_models()
  emm <- emmeans(mods$nb, ~ Stim * Satiety | Sex * Diet, type = "response")
  pw  <- pairs(emm, by = c("Sex", "Diet"), reverse = TRUE, infer = TRUE)

  tab <- bbmake_pairwise_table(pw)

  expect_s3_class(tab, "tbl_df")
  expect_true(nrow(tab) > 0,
              info = paste("Table was empty (0 rows)"))
  expect_true(all(c("Rate Ratio", "RR 95% CI") %in% names(tab)))
})

# ==================================================================
# Integration & polishing
# ==================================================================

test_that("bbmake_pairwise_table output works with bbnice_pairwise_table and has clean polishing", {
  skip_if_not_installed("glmmTMB")
  skip_if_not_installed("rempsyc")
  skip_if_not_installed("flextable")

  mods <- setup_pairwise_models()
  emm <- emmeans(mods$nb, ~ Sex | Diet, type = "response")
  pw  <- pairs(emm, reverse = TRUE, infer = TRUE)

  tab <- bbmake_pairwise_table(pw)

  # polishing checks
  expect_true(all(str_trim(tab$contrast) == tab$contrast),
              info = "contrast column should have no leading/trailing whitespace")

  unsorted_p <- na.omit(tab$p.value)
  expect_true(!is.unsorted(unsorted_p),
              info = paste("p.value not sorted ascending. First few p-values:",
                           paste(head(unsorted_p, 6), collapse = ", ")))

  # nice table should not error
  expect_no_error({
    ft <- bbnice_pairwise_table(tab, title = "Test Pairwise")
    expect_s3_class(ft, "flextable")
  })
})
