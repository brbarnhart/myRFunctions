# tests/testthat/test-bbmake_pairwise_sensitivity_table.R

library(testthat)
library(dplyr)
library(tibble)
library(emmeans)

# ── Oracle: the original hand-rolled pipeline ────────────────────────────────
get_pairs_manual <- function(
  model,
  model_name,
  adjust = "tukey",
  cross.adjust = "none"
) {
  pw <- emmeans(model, ~ Stim | Diet * Sex, type = "response") |>
    pairs(by = c("Sex", "Diet"), reverse = TRUE, adjust = adjust)
  sm <- if (identical(cross.adjust, "none")) {
    summary(pw, infer = c(TRUE, TRUE), adjust = adjust)
  } else {
    summary(pw, infer = c(TRUE, TRUE), by = NULL, adjust = cross.adjust)
  }
  as.data.frame(sm) |>
    mutate(Model = model_name) |>
    select(Sex, Diet, contrast, ratio, SE, p.value, Model)
}

setup_sensitivity_models <- function() {
  set.seed(1)
  dat <- expand.grid(
    Sex  = factor(c("F", "M")),
    Stim = factor(c("A", "B")),
    Diet = factor(c("C", "H")),
    id   = factor(1:8)
  )
  dat$y <- rpois(nrow(dat), lambda = 5)

  list(
    data = dat,
    Full = glm(y ~ Sex * Stim * Diet, data = dat, family = poisson),
    Reduced = glm(
      y ~ Sex * Stim * Diet,
      data = dat[dat$id != "1", ],
      family = poisson
    )
  )
}

# ==================================================================
# Core: named list matches the copy-pasted pipeline
# ==================================================================

test_that("named list of glm models matches the imap_dfr pairs pipeline", {
  s <- setup_sensitivity_models()
  mods <- list(Full = s$Full, `Outliers removed` = s$Reduced)

  expected <- dplyr::bind_rows(
    get_pairs_manual(s$Full, "Full", adjust = "none", cross.adjust = "holm"),
    get_pairs_manual(
      s$Reduced, "Outliers removed",
      adjust = "none", cross.adjust = "holm"
    )
  ) |>
    mutate(contrast = stringr::str_trim(contrast)) |>
    arrange(Sex, Diet, contrast, Model)

  expect_message(
    tab <- bbmake_pairwise_sensitivity_table(
      mods,
      ~ Stim | Diet * Sex,
      by = c("Sex", "Diet")
    ),
    "Holm across 4 comparisons within each model"
  )

  expect_s3_class(tab, "tbl_df")
  expect_equal(
    names(tab)[1:4],
    c("Sex", "Diet", "contrast", "Model")
  )
  expect_true(names(tab)[5] %in% c("z.ratio", "t.ratio"))
  expect_true(all(c("IRR", "SE", "lower.CL", "upper.CL", "p.value") %in% names(tab)))
  expect_equal(as.character(tab$contrast), expected$contrast)
  expect_equal(as.character(tab$Model), expected$Model)
  expect_equal(tab$IRR, expected$ratio)
  expect_equal(tab$SE, expected$SE)
  expect_equal(tab$p.value, expected$p.value)
  expect_equal(tab$Sex, expected$Sex)
  expect_equal(tab$Diet, expected$Diet)
})

test_that("Model column preserves list names and order", {
  s <- setup_sensitivity_models()
  mods <- list(`Outliers removed` = s$Reduced, Full = s$Full)
  tab <- bbmake_pairwise_sensitivity_table(
    mods,
    ~ Stim | Diet * Sex,
    by = c("Sex", "Diet")
  )

  expect_equal(levels(tab$Model), c("Outliers removed", "Full"))
  expect_true(all(c("Outliers removed", "Full") %in% as.character(tab$Model)))
})

test_that("each model contributes the same contrasts", {
  s <- setup_sensitivity_models()
  mods <- list(Full = s$Full, Reduced = s$Reduced)
  tab <- bbmake_pairwise_sensitivity_table(
    mods,
    ~ Stim | Diet * Sex,
    by = c("Sex", "Diet")
  )

  n_per_model <- table(tab$Model)
  expect_equal(unname(n_per_model[["Full"]]), unname(n_per_model[["Reduced"]]))
  expect_equal(n_per_model[["Full"]] * 2L, nrow(tab))
})

# ==================================================================
# Specs / by defaults
# ==================================================================

test_that("by defaults to emmeans by-variables from specs", {
  s <- setup_sensitivity_models()
  mods <- list(Full = s$Full, Reduced = s$Reduced)

  tab <- bbmake_pairwise_sensitivity_table(mods, ~ Stim | Diet * Sex)

  expect_true(all(c("Diet", "Sex") %in% names(tab)))
  expect_true("contrast" %in% names(tab))
})

test_that("two-way specs work without a Diet by-variable", {
  s <- setup_sensitivity_models()
  mods <- list(
    Full = glm(y ~ Sex * Stim, data = s$data, family = poisson),
    Reduced = glm(
      y ~ Sex * Stim,
      data = s$data[s$data$id != "1", ],
      family = poisson
    )
  )

  tab <- bbmake_pairwise_sensitivity_table(mods, ~ Stim | Sex, by = "Sex")

  expect_equal(names(tab)[1:3], c("Sex", "contrast", "Model"))
  expect_true(all(c("IRR", "SE", "p.value") %in% names(tab)))
  expect_false("Diet" %in% names(tab))
})

# ==================================================================
# Gaussian estimate column
# ==================================================================

test_that("gaussian models keep Mean Difference instead of IRR", {
  s <- setup_sensitivity_models()
  set.seed(1)
  dat <- s$data
  dat$yg <- rnorm(nrow(dat))
  mods <- list(
    Full = glm(yg ~ Sex * Stim * Diet, data = dat),
    Reduced = glm(yg ~ Sex * Stim * Diet, data = dat[dat$id != "1", ])
  )

  tab <- bbmake_pairwise_sensitivity_table(
    mods,
    ~ Stim | Diet * Sex,
    by = c("Sex", "Diet")
  )

  expect_true("Mean Difference" %in% names(tab))
  expect_false("IRR" %in% names(tab))
  expect_false("d" %in% names(tab))
  expect_true(is.numeric(tab$`Mean Difference`))
  expect_true(all(c("SE", "lower.CL", "upper.CL") %in% names(tab)))
})

# ==================================================================
# input checks
# ==================================================================

test_that("invalid inputs error clearly", {
  s <- setup_sensitivity_models()
  expect_error(bbmake_pairwise_sensitivity_table(), "`models` is missing")
  expect_error(
    bbmake_pairwise_sensitivity_table(s$Full, ~ Stim | Sex),
    "named list of models"
  )
  expect_error(
    bbmake_pairwise_sensitivity_table(list(s$Full, s$Reduced), ~ Stim | Sex),
    "named list"
  )
  expect_error(
    bbmake_pairwise_sensitivity_table(list(Full = s$Full, s$Reduced), ~ Stim | Sex),
    "named list"
  )
  expect_error(
    bbmake_pairwise_sensitivity_table(list(Full = s$Full)),
    "`specs` is missing"
  )
})

# ==================================================================
# adjust (multiple comparisons)
# ==================================================================

test_that("adjust is forwarded to contrast/pairs", {
  set.seed(1)
  dat <- expand.grid(
    Stim = factor(c("A", "B", "C")),
    id = 1:20
  )
  dat$y <- rpois(nrow(dat), lambda = 5)
  mods <- list(
    Full = glm(y ~ Stim, data = dat, family = poisson),
    Reduced = glm(y ~ Stim, data = dat[-1, ], family = poisson)
  )

  none <- bbmake_pairwise_sensitivity_table(
    mods, ~ Stim, adjust = "none", cross.adjust = "none"
  )
  bonf <- bbmake_pairwise_sensitivity_table(
    mods, ~ Stim, adjust = "bonferroni", cross.adjust = "none"
  )
  tukey <- bbmake_pairwise_sensitivity_table(
    mods, ~ Stim, adjust = "tukey", cross.adjust = "none"
  )
  expect_message(
    default <- bbmake_pairwise_sensitivity_table(mods, ~ Stim),
    "Holm across 3 comparisons within each model"
  )

  expect_false(isTRUE(all.equal(default$p.value, tukey$p.value)))
  expect_true(all(default$p.value >= none$p.value - 1e-12))
  expect_true(all(bonf$p.value >= none$p.value - 1e-12))
  expect_true(any(bonf$p.value > none$p.value + 1e-12))

  emm <- emmeans(mods$Full, ~ Stim, type = "response")
  expected <- as.data.frame(
    contrast(emm, method = "revpairwise", adjust = "bonferroni")
  )
  got <- bonf[as.character(bonf$Model) == "Full", ]
  expect_equal(got$p.value, expected$p.value)
  expect_equal(got$IRR, expected$ratio)
})

test_that("adjust matches the pairs() oracle on the factorial glm", {
  s <- setup_sensitivity_models()
  mods <- list(Full = s$Full, Reduced = s$Reduced)

  expected <- dplyr::bind_rows(
    get_pairs_manual(s$Full, "Full", adjust = "bonferroni"),
    get_pairs_manual(s$Reduced, "Reduced", adjust = "bonferroni")
  ) |>
    mutate(contrast = stringr::str_trim(contrast)) |>
    arrange(Sex, Diet, contrast, Model)

  tab <- bbmake_pairwise_sensitivity_table(
    mods,
    ~ Stim | Diet * Sex,
    by = c("Sex", "Diet"),
    adjust = "bonferroni",
    cross.adjust = "none"
  )

  expect_equal(tab$p.value, expected$p.value)
})

test_that("cross.adjust is forwarded across by-groups", {
  s <- setup_sensitivity_models()
  mods <- list(Full = s$Full, Reduced = s$Reduced)

  none <- bbmake_pairwise_sensitivity_table(
    mods, ~ Stim | Diet * Sex,
    by = c("Sex", "Diet"),
    adjust = "none",
    cross.adjust = "none"
  )
  cross <- bbmake_pairwise_sensitivity_table(
    mods, ~ Stim | Diet * Sex,
    by = c("Sex", "Diet"),
    adjust = "none",
    cross.adjust = "bonferroni"
  )

  expect_message(
    default <- bbmake_pairwise_sensitivity_table(
      mods, ~ Stim | Diet * Sex,
      by = c("Sex", "Diet")
    ),
    "Holm across 4 comparisons within each model"
  )
  expect_true(all(default$p.value >= none$p.value - 1e-12))
  expect_true(any(default$p.value > none$p.value + 1e-12))
  expect_true(all(cross$p.value >= none$p.value - 1e-12))
  expect_true(any(cross$p.value > none$p.value + 1e-12))

  expected <- dplyr::bind_rows(
    get_pairs_manual(s$Full, "Full", adjust = "none", cross.adjust = "bonferroni"),
    get_pairs_manual(s$Reduced, "Reduced", adjust = "none", cross.adjust = "bonferroni")
  ) |>
    mutate(contrast = stringr::str_trim(contrast)) |>
    arrange(Sex, Diet, contrast, Model)

  expect_equal(cross$p.value, expected$p.value)
  expect_equal(cross$IRR, expected$ratio)
})

# ==================================================================
# glmmTMB (the actual use case)
# ==================================================================

test_that("glmmTMB nbinom2 models return IRR contrasts", {
  suppressWarnings(skip_if_not_installed("glmmTMB"))

  s <- setup_sensitivity_models()
  dat <- s$data
  mods <- list(
    Full = suppressWarnings(glmmTMB::glmmTMB(
      y ~ Sex * Stim * Diet + (1 | id),
      data = dat,
      family = glmmTMB::nbinom2()
    )),
    Reduced = suppressWarnings(glmmTMB::glmmTMB(
      y ~ Sex * Stim * Diet + (1 | id),
      data = dat[dat$id != "1", ],
      family = glmmTMB::nbinom2()
    ))
  )

  tab <- suppressWarnings(bbmake_pairwise_sensitivity_table(
    mods,
    ~ Stim | Diet * Sex,
    by = c("Sex", "Diet")
  ))

  expect_s3_class(tab, "tbl_df")
  expect_equal(names(tab)[1:4], c("Sex", "Diet", "contrast", "Model"))
  expect_true(all(c("IRR", "SE", "lower.CL", "upper.CL", "p.value") %in% names(tab)))
  expect_true(is.numeric(tab$IRR))
  expect_true(all(tab$IRR > 0))
  expect_equal(levels(tab$Model), c("Full", "Reduced"))
})
