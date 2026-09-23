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

  required_cols <- c(
    "contrast", "Mean Difference", "SE", "lower.CL", "upper.CL",
    "p.value", "d"
  )
  missing_cols <- setdiff(required_cols, names(tab))
  expect_true(length(missing_cols) == 0,
              info = paste("Missing columns:", paste(missing_cols, collapse = ", ")))
  expect_true(any(c("t.ratio", "z.ratio") %in% names(tab)))
  expect_true(is.numeric(tab$`Mean Difference`))
  expect_true(is.numeric(tab$SE))
  expect_true(is.numeric(tab$lower.CL))
  expect_true(is.numeric(tab$d))
  expect_equal(sum(is.na(tab$p.value)), 0L)
})

test_that("bbmake_pairwise_table forwards cross.adjust to summary()", {
  skip_if_not_installed("emmeans")

  set.seed(1)
  dat <- expand.grid(
    Sex  = factor(c("F", "M")),
    Stim = factor(c("Cont", "Opto")),
    Diet = factor(c("Chow", "HFD")),
    id   = 1:8
  )
  dat$y <- rpois(nrow(dat), lambda = 5)
  mod <- glm(y ~ Sex * Stim * Diet, data = dat, family = poisson)
  emm <- emmeans(mod, ~ Stim | Diet * Sex, type = "response")
  pw  <- pairs(emm, reverse = TRUE, adjust = "none")

  none <- bbmake_pairwise_table(pw, cross.adjust = "none")
  holm <- bbmake_pairwise_table(pw, cross.adjust = "holm")
  expected <- as.data.frame(summary(pw, infer = c(TRUE, TRUE), by = NULL, adjust = "holm"))

  expect_true(all(holm$p.value >= none$p.value - 1e-12))
  expect_true(any(holm$p.value > none$p.value + 1e-12))
  expect_equal(sort(holm$p.value), sort(expected$p.value))
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

  required_cols <- c(
    "contrast", "IRR", "SE", "lower.CL", "upper.CL", "p.value"
  )
  missing_cols <- setdiff(required_cols, names(tab))
  expect_true(length(missing_cols) == 0,
              info = paste("Missing columns:", paste(missing_cols, collapse = ", ")))
  expect_true(any(c("z.ratio", "t.ratio") %in% names(tab)))
  expect_true(is.numeric(tab$IRR))
  expect_true(is.numeric(tab$SE))
  expect_true(is.numeric(tab$upper.CL))
  expect_true(is.numeric(tab$lower.CL))
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
  expect_true(all(c("IRR", "SE", "lower.CL", "upper.CL", "p.value") %in% names(tab)))
  expect_true(any(c("z.ratio", "t.ratio") %in% names(tab)))
})

# ==================================================================
# Integration & polishing
# ==================================================================

test_that("bbmake_pairwise_table output works with bbnice_pairwise_table and has clean polishing", {
  skip_if_not_installed("glmmTMB")
  skip_if_not_installed("rempsyc")
  skip_if_not_installed("flextable")

  mods <- setup_pairwise_models()
  emm <- emmeans(mods$nb, ~ Sex * Diet, type = "response")
  pw  <- pairs(emm, reverse = TRUE, infer = TRUE)

  tab <- bbmake_pairwise_table(pw)

  unsorted_p <- na.omit(tab$p.value)
  expect_true(!is.unsorted(unsorted_p),
              info = paste("p.value not sorted ascending. First few p-values:",
                           paste(head(unsorted_p, 6), collapse = ", ")))

  ft <- bbnice_pairwise_table(tab)
  stat <- if ("z.ratio" %in% names(tab)) "z" else "t"
  expect_s3_class(ft, "flextable")
  expect_equal(ft$col_keys, c("Contrast", stat, "IRR", "95% CI", "p"))
})
