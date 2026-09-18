# tests/testthat/test-bbmake_lrt_sensitivity_table.R

library(testthat)
library(dplyr)
library(tibble)

setup_lrt_sensitivity_models <- function() {
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
# Core: stacked bbmake_lrt_table
# ==================================================================

test_that("named list matches bind_rows of bbmake_lrt_table", {
  s <- setup_lrt_sensitivity_models()
  mods <- list(Full = s$Full, `Outliers removed` = s$Reduced)

  expected <- dplyr::bind_rows(
    dplyr::mutate(bbmake_lrt_table(s$Full), Model = "Full"),
    dplyr::mutate(bbmake_lrt_table(s$Reduced), Model = "Outliers removed")
  )
  expected$Term <- factor(expected$Term, levels = unique(expected$Term))
  expected$Model <- factor(expected$Model, levels = names(mods))
  expected <- dplyr::select(expected, Term, Model, Df, LRT, p)
  expected <- dplyr::arrange(expected, Term, Model)

  tab <- bbmake_lrt_sensitivity_table(mods)

  expect_s3_class(tab, "tbl_df")
  expect_equal(names(tab), c("Term", "Model", "Df", "LRT", "p"))
  expect_equal(as.character(tab$Term), as.character(expected$Term))
  expect_equal(as.character(tab$Model), as.character(expected$Model))
  expect_equal(tab$Df, expected$Df)
  expect_equal(tab$LRT, expected$LRT)
  expect_equal(tab$p, expected$p)
})

test_that("Model column preserves list names and order", {
  s <- setup_lrt_sensitivity_models()
  mods <- list(`Outliers removed` = s$Reduced, Full = s$Full)
  tab <- bbmake_lrt_sensitivity_table(mods)

  expect_equal(levels(tab$Model), c("Outliers removed", "Full"))
  expect_true(all(c("Outliers removed", "Full") %in% as.character(tab$Model)))
})

test_that("each model contributes the same terms in Type II order", {
  s <- setup_lrt_sensitivity_models()
  mods <- list(Full = s$Full, Reduced = s$Reduced)
  tab <- bbmake_lrt_sensitivity_table(mods)
  single <- bbmake_lrt_table(s$Full)

  n_per_model <- table(tab$Model)
  expect_equal(unname(n_per_model[["Full"]]), unname(n_per_model[["Reduced"]]))
  expect_equal(n_per_model[["Full"]] * 2L, nrow(tab))
  expect_equal(
    as.character(tab$Term[as.character(tab$Model) == "Full"]),
    single$Term
  )
  expect_equal(as.character(tab$Term[[1]]), "Sex:Stim:Diet")
})

test_that("rows are grouped by term so models sit next to each other", {
  s <- setup_lrt_sensitivity_models()
  mods <- list(Full = s$Full, Reduced = s$Reduced)
  tab <- bbmake_lrt_sensitivity_table(mods)

  first_term <- as.character(tab$Term[1:2])
  first_model <- as.character(tab$Model[1:2])
  expect_equal(first_term, c("Sex:Stim:Diet", "Sex:Stim:Diet"))
  expect_equal(first_model, c("Full", "Reduced"))
})

# ==================================================================
# digits and input checks
# ==================================================================

test_that("digits rounds LRT and p only", {
  s <- setup_lrt_sensitivity_models()
  mods <- list(Full = s$Full, Reduced = s$Reduced)
  raw <- bbmake_lrt_sensitivity_table(mods)
  rnd <- bbmake_lrt_sensitivity_table(mods, digits = 3)

  expect_equal(rnd$LRT, round(raw$LRT, 3))
  expect_equal(rnd$p, round(raw$p, 3))
  expect_equal(as.character(rnd$Term), as.character(raw$Term))
})

test_that("invalid inputs error clearly", {
  s <- setup_lrt_sensitivity_models()
  expect_error(bbmake_lrt_sensitivity_table(), "`models` is missing")
  expect_error(
    bbmake_lrt_sensitivity_table(s$Full),
    "named list of models"
  )
  expect_error(
    bbmake_lrt_sensitivity_table(list(s$Full, s$Reduced)),
    "named list"
  )
  expect_error(
    bbmake_lrt_sensitivity_table(list(Full = s$Full), digits = -1),
    "`digits`"
  )
})

test_that("note attribute records stacked Type II LRTs", {
  s <- setup_lrt_sensitivity_models()
  mods <- list(Full = s$Full, Reduced = s$Reduced)
  tab <- bbmake_lrt_sensitivity_table(mods)
  expect_match(attr(tab, "note"), "Type II")
  expect_match(attr(tab, "note"), "stacked")
})

# ==================================================================
# glmmTMB (the actual use case)
# ==================================================================

test_that("glmmTMB two-way with random intercept stacks LRTs", {
  suppressWarnings(skip_if_not_installed("glmmTMB"))

  set.seed(1)
  dat <- expand.grid(
    Sex  = factor(c("F", "M")),
    Stim = factor(c("A", "B")),
    id   = factor(1:15)
  )
  dat$yb <- pmin(pmax(rbeta(nrow(dat), 2, 3), 1e-3), 1 - 1e-3)

  full <- suppressWarnings(glmmTMB::glmmTMB(
    yb ~ Sex * Stim + (1 | id),
    data = dat,
    family = glmmTMB::beta_family()
  ))
  reduced <- suppressWarnings(glmmTMB::glmmTMB(
    yb ~ Sex * Stim + (1 | id),
    data = dat[dat$id != "1", ],
    family = glmmTMB::beta_family()
  ))
  mods <- list(Full = full, Reduced = reduced)

  expected <- suppressWarnings(dplyr::bind_rows(
    dplyr::mutate(bbmake_lrt_table(full), Model = "Full"),
    dplyr::mutate(bbmake_lrt_table(reduced), Model = "Reduced")
  ))
  expected$Term <- factor(expected$Term, levels = unique(expected$Term))
  expected$Model <- factor(expected$Model, levels = names(mods))
  expected <- dplyr::arrange(expected, Term, Model)

  tab <- suppressWarnings(bbmake_lrt_sensitivity_table(mods))

  expect_equal(as.character(tab$Term), as.character(expected$Term))
  expect_equal(tab$LRT, expected$LRT)
  expect_equal(tab$p, expected$p)
  expect_false(any(grepl("|", as.character(tab$Term), fixed = TRUE)))
  expect_equal(nlevels(tab$Model), 2L)
})
