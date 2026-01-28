# tests/testthat/test-bbnice_anova_table.R

library(testthat)
library(lme4)
library(lmerTest)
library(effectsize)
library(dplyr)
library(tibble)
library(flextable)
library(rempsyc)

setup_models <- function() {
  set.seed(42)

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

  mod_lm   <- lm(Total_Correct_Pokes ~ Sex * Diet * Stim, data = dat)
  mod_lmer <- lmer(Total_Correct_Pokes ~ Sex * Diet * Stim + (1 | ID), data = dat)

  list(lm = mod_lm, lmer = mod_lmer, data = dat)
}


test_that("bbnice_anova_table produces valid flextable and keeps p numeric", {
  skip_if_not_installed("rempsyc")
  skip_if_not_installed("flextable")
  skip_if_not_installed("lmerTest")

  mods <- setup_models()

  # Real lm model
  tab_lm <- bbmake_anova_table(mods$lm)
  expect_true("Pr(>F)" %in% names(tab_lm))
  expect_true(is.numeric(tab_lm$`Pr(>F)`))

  ft_lm <- bbnice_anova_table(tab_lm, title = "Test LM")
  expect_s3_class(ft_lm, "flextable")
  expect_true(nrow(ft_lm$body$dataset) >= 1)
  expect_true(ncol(ft_lm$body$dataset) >= 5)

  # Real lmer model
  tab_lmer <- bbmake_anova_table(mods$lmer)
  expect_true("Pr(>F)" %in% names(tab_lmer))
  expect_true(is.numeric(tab_lmer$`Pr(>F)`))

  ft_lmer <- bbnice_anova_table(tab_lmer, title = "Test LMER")
  expect_s3_class(ft_lmer, "flextable")
  expect_true(nrow(ft_lmer$body$dataset) >= 1)
})


test_that("bbnice_anova_table adds stars when p-values are small enough", {
  skip_if_not_installed("rempsyc")

  # Forced small-p table → must produce stars
  small_p_tab <- tibble(
    Effect             = c("Main", "Interaction", "Small"),
    NumDF              = c(1, 2, 1),
    DenDF              = c(NA, NA, NA),
    `F value`          = c(6.8, 11.2, 4.9),
    `Pr(>F)`           = c(0.012, 0.0008, 0.032),          # <0.05, <0.001, <0.05
    `Omega2 (partial)` = c(0.18, 0.29, 0.11),
    `Omega2 (95% CI)`  = c("[0.05, 0.33]", "[0.12, 0.47]", "[0.01, 0.23]")
  )

  ft <- bbnice_anova_table(small_p_tab, title = "Small P Test")

  expect_s3_class(ft, "flextable")

  body_text <- as.character(unlist(ft$body$dataset))
  print(ft)
  expect_true(any(grepl("\\*+", body_text)))           # Should have at least one *
})


test_that("bbnice_anova_table errors when p-column is character", {
  skip_if_not_installed("rempsyc")

  bad_table <- tibble(
    Effect             = c("A", "B"),
    NumDF              = c(1, 1),
    DenDF              = c(NA, NA),
    `F value`          = c(5.2, 3.1),
    `Pr(>F)`           = c("<0.001", "0.042"),            # ← character
    `Omega2 (partial)` = c(0.22, 0.09),
    `Omega2 (95% CI)`  = c("[0.08, 0.38]", "[0.00, 0.20]")
  )

  expect_error(
    bbnice_anova_table(bad_table, title = "Bad P"),
    regexp = "p value must be numeric"
  )
})


test_that("bbnice_anova_table handles missing Pr(>F) column gracefully", {
  skip_if_not_installed("rempsyc")

  no_p_table <- tibble(
    Effect             = c("Main", "Other"),
    NumDF              = c(1, 2),
    DenDF              = c(NA, NA),
    `F value`          = c(7.5, 4.2),
    `Omega2 (partial)` = c(0.19, 0.13),
    `Omega2 (95% CI)`  = c("[0.06, 0.34]", "[0.02, 0.26]")
  )

  expect_no_error({
    ft <- bbnice_anova_table(no_p_table, title = "No P")
    expect_s3_class(ft, "flextable")
  })
})
