# tests/testthat/test-bb_pairwise_labels.R

library(testthat)
library(emmeans)
library(lme4)
library(lmerTest)
library(dplyr)
library(tibble)

setup_pairwise <- function() {
  set.seed(42)
  dat <- expand.grid(
    subj  = factor(1:12),
    Stim  = factor(c("Low", "Med", "High"), levels = c("Low", "Med", "High")),
    Group = factor(c("Control", "Treatment"))
  ) |>
    as_tibble() |>
    mutate(
      y = 10 +
        3 * (Stim == "Med") +
        6 * (Stim == "High") +
        4 * (Group == "Treatment") +
        2.5 * (Stim == "High" & Group == "Treatment") +
        rnorm(n(), sd = 2)
    )

  mod <- suppressWarnings(lmer(y ~ Stim * Group + (1 | subj), data = dat))
  emm <- emmeans(mod, ~ Stim | Group)
  pw  <- pairs(emm, adjust = "tukey")
  list(mod = mod, emm = emm, pw = pw)
}

test_that("bb_pairwise_labels returns group1/group2, p.signif, y.position", {
  skip_if_not_installed("lmerTest")
  skip_if_not_installed("emmeans")

  s <- setup_pairwise()
  sig <- bb_pairwise_labels(s$pw, emm = s$emm, model = s$mod)

  expect_s3_class(sig, "tbl_df")
  expect_true(all(c("group1", "group2", "p.signif", "y.position",
                    "effect_annotation", "p.value") %in% names(sig)))
  expect_true("Group" %in% names(sig))
  expect_true(all(sig$p.signif %in% c("ns", "*", "**", "***")))
  expect_true(all(is.finite(sig$y.position)))
  expect_equal(attr(sig, "by.vars"), "Group")

  # group labels come from contrast strings
  expect_true(all(sig$group1 %in% c("Low", "Med", "High")))
  expect_true(all(sig$group2 %in% c("Low", "Med", "High")))
})

test_that("bb_pairwise_labels stacks y.position within panels", {
  skip_if_not_installed("lmerTest")
  skip_if_not_installed("emmeans")

  s <- setup_pairwise()
  sig <- bb_pairwise_labels(s$pw, emm = s$emm, step = 0.1)

  # Multiple contrasts per Group should not all share the exact same y
  n_unique <- sig |>
    group_by(Group) |>
    summarise(n_y = n_distinct(y.position), .groups = "drop")
  expect_true(all(n_unique$n_y > 1))
})

test_that("bb_pairwise_labels works with prebuilt pw_table (Cohen's d column)", {
  skip_if_not_installed("emmeans")

  s <- setup_pairwise()
  pw_table <- summary(s$pw) |>
    as_tibble() |>
    select(contrast, Group, p.value) |>
    mutate(`Cohen's d` = 0.5)

  sig <- bb_pairwise_labels(pw_table = pw_table, emm = s$emm)
  expect_true(all(grepl("Cohen's d", sig$effect_annotation)))
  expect_true(all(c("group1", "group2") %in% names(sig)))
})

test_that("bb_pairwise_labels hide.ns drops non-significant rows", {
  skip_if_not_installed("lmerTest")
  skip_if_not_installed("emmeans")

  s <- setup_pairwise()
  sig_all <- bb_pairwise_labels(s$pw, emm = s$emm, hide.ns = FALSE)
  sig_sig <- bb_pairwise_labels(s$pw, emm = s$emm, hide.ns = TRUE)

  expect_true(nrow(sig_sig) <= nrow(sig_all))
  expect_true(all(sig_sig$p.signif != "ns"))
})

test_that("bb_pairwise_labels errors without pw or pw_table", {
  expect_error(bb_pairwise_labels(), "pw")
})
