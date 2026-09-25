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

test_that("default adjustment is Holm across every comparison", {
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
  pw  <- pairs(emm, reverse = TRUE, adjust = "tukey")
  expected <- as.data.frame(summary(
    pw, infer = c(TRUE, TRUE), by = NULL, adjust = "holm"
  ))

  expect_message(
    tab <- bbmake_pairwise_table(pw),
    "Holm across 4 comparisons"
  )
  expect_equal(sort(tab$p.value), sort(expected$p.value))
  expect_silent(bbmake_pairwise_table(pw, cross.adjust = "none"))
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

# ==================================================================
# Interaction contrasts
# ==================================================================

interaction_poisson <- function() {
  set.seed(1)
  dat <- expand.grid(
    Sex = factor(c("F", "M")),
    Stim = factor(c("A", "B", "C")),
    Diet = factor(c("C", "H")),
    id = 1:6
  )
  dat$y <- rpois(nrow(dat), lambda = 5)
  mod <- glm(y ~ Sex * Stim * Diet, data = dat, family = poisson)
  emm <- emmeans(mod, ~ Stim | Diet * Sex, type = "response")
  list(
    mod = mod,
    pw = contrast(emm, interaction = "pairwise", by = "Sex")
  )
}

test_that("simple pairs() tables keep the existing column order", {
  skip_if_not_installed("emmeans")

  set.seed(1)
  dat <- expand.grid(
    Sex = factor(c("F", "M")),
    Stim = factor(c("A", "B")),
    id = 1:8
  )
  dat$y <- rpois(nrow(dat), lambda = 5)
  mod <- glm(y ~ Sex * Stim, data = dat, family = poisson)
  pw <- pairs(emmeans(mod, ~ Stim | Sex, type = "response"), reverse = TRUE)

  tab <- bbmake_pairwise_table(pw, cross.adjust = "none")
  expect_equal(
    names(tab),
    c(
      "Sex", "contrast", "df", "z.ratio", "IRR",
      "lower.CL", "upper.CL", "SE", "p.value"
    )
  )
})

test_that("interaction contrasts keep one column per factor", {
  skip_if_not_installed("emmeans")

  built <- interaction_poisson()
  expected <- as.data.frame(summary(
    built$pw,
    infer = c(TRUE, TRUE),
    by = NULL,
    adjust = "holm"
  ))

  expect_message(
    tab <- bbmake_pairwise_table(built$pw),
    "Holm across 6 comparisons"
  )
  expect_equal(
    names(tab),
    c(
      "Sex", "Stim_pairwise", "Diet_pairwise", "df", "z.ratio", "IRR",
      "lower.CL", "upper.CL", "SE", "p.value"
    )
  )
  expect_false("contrast" %in% names(tab))
  expect_false("null" %in% names(tab))
  expect_equal(sort(tab$p.value), sort(expected$p.value))
  expect_equal(as.character(unique(tab$Sex)), c("F", "M"))
})

test_that("interaction tables pipe into bbnice_pairwise_table", {
  skip_if_not_installed("emmeans")
  skip_if_not_installed("flextable")
  skip_if_not_installed("rempsyc")

  built <- interaction_poisson()
  ft <- suppressMessages(bbmake_pairwise_table(built$pw)) |>
    bbnice_pairwise_table()

  expect_equal(
    ft$col_keys,
    c("Sex", "Stim Pairwise", "Diet Pairwise", "z", "IRR", "95% CI", "p")
  )
  spans <- as.numeric(ft$body$spans$columns[, match("Sex", ft$col_keys)])
  expect_equal(spans, c(3, 0, 0, 3, 0, 0))
})

test_that("revpairwise interaction columns are kept", {
  skip_if_not_installed("emmeans")

  set.seed(1)
  dat <- expand.grid(
    Sex = factor(c("F", "M")),
    Stim = factor(c("A", "B")),
    Diet = factor(c("C", "H")),
    id = 1:6
  )
  dat$y <- rpois(nrow(dat), lambda = 5)
  mod <- glm(y ~ Sex * Stim * Diet, data = dat, family = poisson)
  pw <- emmeans(mod, ~ Stim | Diet * Sex, type = "response") |>
    contrast(interaction = "revpairwise", by = "Sex")

  tab <- bbmake_pairwise_table(pw, cross.adjust = "none")
  expect_true(all(c("Sex", "Stim_revpairwise", "Diet_revpairwise") %in% names(tab)))
  expect_false("contrast" %in% names(tab))
})

test_that("interaction contrasts do not attach Cohen's d", {
  skip_if_not_installed("emmeans")

  set.seed(1)
  dat <- expand.grid(
    Sex = factor(c("F", "M")),
    Stim = factor(c("A", "B")),
    Diet = factor(c("C", "H")),
    id = 1:8
  )
  dat$y <- rnorm(nrow(dat))
  mod <- lm(y ~ Sex * Stim * Diet, data = dat)
  pw <- emmeans(mod, ~ Stim | Diet * Sex) |>
    contrast(interaction = "pairwise", by = "Sex")

  expect_no_warning(
    tab <- bbmake_pairwise_table(pw, model = mod, cross.adjust = "none")
  )
  expect_true("Mean Difference" %in% names(tab))
  expect_false("d" %in% names(tab))
  expect_false("contrast" %in% names(tab))
})

test_that("a contrast grid without pairwise columns still requires contrast", {
  skip_if_not_installed("emmeans")

  set.seed(1)
  dat <- expand.grid(
    Sex = factor(c("F", "M")),
    Stim = factor(c("A", "B")),
    id = 1:8
  )
  dat$y <- rpois(nrow(dat), lambda = 5)
  mod <- glm(y ~ Sex * Stim, data = dat, family = poisson)
  pw <- emmeans(mod, ~ Stim | Sex, type = "response") |>
    contrast(interaction = "consec")

  expect_error(bbmake_pairwise_table(pw), "contrast")
})
