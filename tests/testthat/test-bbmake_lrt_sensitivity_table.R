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

# ==================================================================
# keep: forward a blocking term only where the model contains it
# ==================================================================

capture_sensitivity <- function(expr) {
  msgs <- character()
  tab <- withCallingHandlers(
    expr,
    message = function(m) {
      msgs <<- c(msgs, conditionMessage(m))
      invokeRestart("muffleMessage")
    }
  )
  list(tab = tab, message = paste(msgs, collapse = "\n"))
}

setup_block_sensitivity <- function() {
  set.seed(1)
  dat <- expand.grid(
    Sex  = factor(c("F", "M")),
    Stim = factor(c("A", "B")),
    Diet = factor(c("C", "H")),
    id   = factor(1:6)
  )
  dat$y <- rpois(nrow(dat), lambda = 5)
  dat$Block <- dat$id
  dat <- dat[!(dat$Block %in% c("1", "2") & dat$Sex == "F"), ]
  list(
    both = glm(
      y ~ Block + Diet + Sex * Stim, data = dat, family = poisson
    ),
    block_only = glm(
      y ~ Block + Sex * Stim, data = dat, family = poisson
    ),
    neither = glm(y ~ Sex * Stim, data = dat, family = poisson)
  )
}

lrt_rows_for <- function(tab, model) {
  rows <- tab[as.character(tab$Model) == model, c("Term", "Df", "LRT", "p")]
  rows$Term <- as.character(rows$Term)
  rows[order(rows$Term), ]
}

expect_same_lrt <- function(got, expected) {
  expected <- expected[order(expected$Term), ]
  expect_equal(got$Term, expected$Term)
  expect_equal(got$Df, expected$Df)
  expect_equal(got$LRT, expected$LRT)
  expect_equal(got$p, expected$p)
}

test_that("models without an additive term stay quiet", {
  s <- setup_lrt_sensitivity_models()
  mods <- list(Full = s$Full, Reduced = s$Reduced)
  expect_silent(bbmake_lrt_sensitivity_table(mods))
})

test_that("keep is applied only in models that contain the term", {
  s <- setup_block_sensitivity()
  mods <- list(
    Both = s$both,
    `Block only` = s$block_only,
    Neither = s$neither
  )
  out <- capture_sensitivity(
    bbmake_lrt_sensitivity_table(mods, keep = c("Block", "Diet"))
  )

  both <- suppressMessages(bbmake_lrt_table(s$both, keep = c("Block", "Diet")))
  block_only <- suppressMessages(bbmake_lrt_table(s$block_only, keep = "Block"))
  block_only_dropped <- suppressMessages(bbmake_lrt_table(s$block_only))
  neither <- bbmake_lrt_table(s$neither)

  expect_same_lrt(lrt_rows_for(out$tab, "Both"), both)
  expect_same_lrt(lrt_rows_for(out$tab, "Block only"), block_only)
  expect_same_lrt(lrt_rows_for(out$tab, "Neither"), neither)
  expect_false(isTRUE(all.equal(
    block_only$LRT[block_only$Term == "Sex"],
    block_only_dropped$LRT[block_only_dropped$Term == "Sex"]
  )))
  expect_false("Diet" %in% as.character(out$tab$Term[out$tab$Model == "Block only"]))
  expect_false(any(c("Block", "Diet") %in% as.character(
    out$tab$Term[out$tab$Model == "Neither"]
  )))
  expect_equal(
    as.character(out$tab$Model[out$tab$Term == "Diet"]),
    "Both"
  )

  expect_match(out$message, "`Both`:", fixed = TRUE)
  expect_match(out$message, "`Block` and `Diet` stay", fixed = TRUE)
  expect_match(out$message, "`Block only`:", fixed = TRUE)
  expect_match(
    out$message,
    "`Neither` does not contain `Block`, so that model is fit without it.",
    fixed = TRUE
  )
  expect_match(
    out$message,
    "`Block only` and `Neither` do not contain `Diet`, so those models are fit without it.",
    fixed = TRUE
  )
  expect_match(
    attr(out$tab, "note"),
    "Block, Diet retained in models that contain them",
    fixed = TRUE
  )
  expect_match(attr(out$tab, "note"), "stacked", fixed = TRUE)
})

test_that("omitting keep still names the model that drops a blocking factor", {
  s <- setup_block_sensitivity()
  out <- capture_sensitivity(
    bbmake_lrt_sensitivity_table(list(`Block only` = s$block_only))
  )

  expect_match(out$message, "`Block only`:", fixed = TRUE)
  expect_match(out$message, "Later LRTs are fit without `Block`", fixed = TRUE)
  expect_match(out$message, 'Pass keep = "Block"', fixed = TRUE)
  expect_false(grepl("retained", attr(out$tab, "note"), fixed = TRUE))
})

test_that("invalid keep fails before the models are stacked", {
  s <- setup_block_sensitivity()
  mods <- list(Both = s$both, Neither = s$neither)
  expect_error(
    bbmake_lrt_sensitivity_table(mods, keep = "Cohort"),
    "not in any model"
  )
  expect_error(
    bbmake_lrt_sensitivity_table(mods, keep = c("Block", "Cohort")),
    "Cohort"
  )
  expect_error(
    bbmake_lrt_sensitivity_table(mods, keep = 1),
    "`keep` must be"
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
