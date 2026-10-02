# tests/testthat/test-bbmake_lrt_table.R

library(testthat)
library(dplyr)
library(tibble)

# ── Oracle: the original hand-rolled pipeline ────────────────────────────────
lrt_drop1_manual <- function(model) {
  drop1(model, test = "Chisq") |>
    as.data.frame() |>
    tibble::rownames_to_column("Term") |>
    dplyr::filter(Term != "<none>")
}

setup_poisson <- function() {
  set.seed(1)
  dat <- expand.grid(
    Sex  = factor(c("F", "M")),
    Stim = factor(c("A", "B")),
    Diet = factor(c("C", "H")),
    id   = factor(1:6)
  )
  dat$y <- rpois(nrow(dat), lambda = 5)
  dat
}

# ==================================================================
# Core: two-way model matches the copy-pasted pipeline
# ==================================================================

test_that("two-way glm matches sequential drop1 bind_rows pipeline", {
  dat <- setup_poisson()
  mod2 <- glm(y ~ Sex * Stim, data = dat, family = poisson)
  mod1 <- update(mod2, . ~ . - Sex:Stim)

  expected <- dplyr::bind_rows(
    lrt_drop1_manual(mod2),
    lrt_drop1_manual(mod1)
  ) |>
    dplyr::select(Term, Df, LRT, `Pr(>Chi)`) |>
    dplyr::rename(p = `Pr(>Chi)`)

  tab <- bbmake_lrt_table(mod2)

  expect_s3_class(tab, "tbl_df")
  expect_equal(names(tab), c("Term", "Df", "LRT", "p"))
  expect_equal(tab$Term, expected$Term)
  expect_equal(tab$Df, as.integer(expected$Df))
  expect_equal(tab$LRT, expected$LRT)
  expect_equal(tab$p, expected$p)
  expect_false("<none>" %in% tab$Term)
})

test_that("two-way table tests the interaction before main effects", {
  dat <- setup_poisson()
  mod <- glm(y ~ Sex * Stim, data = dat, family = poisson)
  tab <- bbmake_lrt_table(mod)

  expect_equal(tab$Term, c("Sex:Stim", "Sex", "Stim"))
})

# ==================================================================
# Three-way: one row per term, highest order first
# ==================================================================

test_that("three-way glm walks interaction hierarchy then main effects", {
  dat <- setup_poisson()
  mod3 <- glm(y ~ Sex * Stim * Diet, data = dat, family = poisson)
  mod2 <- update(mod3, . ~ . - Sex:Stim:Diet)
  mod1 <- update(mod2, . ~ . - Sex:Stim - Sex:Diet - Stim:Diet)

  expected <- dplyr::bind_rows(
    lrt_drop1_manual(mod3),
    lrt_drop1_manual(mod2),
    lrt_drop1_manual(mod1)
  )

  tab <- bbmake_lrt_table(mod3)

  expect_equal(tab$Term, expected$Term)
  expect_equal(tab$LRT, expected$LRT)
  expect_equal(tab$p, expected$`Pr(>Chi)`)
  expect_equal(tab$Term[[1]], "Sex:Stim:Diet")
  expect_true(all(c("Sex:Stim", "Sex:Diet", "Stim:Diet") %in% tab$Term))
  expect_true(all(c("Sex", "Stim", "Diet") %in% tab$Term))
})

# ==================================================================
# Edge cases
# ==================================================================

test_that("main-effects-only model is a single drop1 pass", {
  dat <- setup_poisson()
  mod <- glm(y ~ Sex + Stim, data = dat, family = poisson)
  expected <- lrt_drop1_manual(mod)
  tab <- bbmake_lrt_table(mod)

  expect_equal(sort(tab$Term), sort(expected$Term))
  expect_equal(nrow(tab), 2L)
})

test_that("intercept-only model returns an empty table with the right columns", {
  dat <- setup_poisson()
  mod <- glm(y ~ 1, data = dat, family = poisson)
  tab <- bbmake_lrt_table(mod)

  expect_s3_class(tab, "tbl_df")
  expect_equal(names(tab), c("Term", "Df", "LRT", "p"))
  expect_equal(nrow(tab), 0L)
})

test_that("gaussian glm maps scaled deviance onto LRT", {
  dat <- setup_poisson()
  set.seed(1)
  dat$yg <- rnorm(nrow(dat))
  mod <- glm(yg ~ Sex * Stim, data = dat)

  d1 <- as.data.frame(drop1(mod, test = "Chisq"))
  tab <- bbmake_lrt_table(mod)

  expect_equal(tab$Term[[1]], "Sex:Stim")
  expect_equal(tab$LRT[[1]], d1["Sex:Stim", "scaled dev."])
  expect_equal(tab$p[[1]], d1["Sex:Stim", "Pr(>Chi)"])
})

test_that("digits rounds LRT and p only", {
  dat <- setup_poisson()
  mod <- glm(y ~ Sex * Stim, data = dat, family = poisson)
  raw <- bbmake_lrt_table(mod)
  rnd <- bbmake_lrt_table(mod, digits = 3)

  expect_equal(rnd$LRT, round(raw$LRT, 3))
  expect_equal(rnd$p, round(raw$p, 3))
  expect_equal(rnd$Term, raw$Term)
})

test_that("invalid digits and missing model error clearly", {
  expect_error(bbmake_lrt_table(), "`model` is missing")
  dat <- setup_poisson()
  mod <- glm(y ~ Sex, data = dat, family = poisson)
  expect_error(bbmake_lrt_table(mod, digits = -1), "`digits`")
  expect_error(bbmake_lrt_table(mod, digits = c(1, 2)), "`digits`")
})

test_that("note attribute records Type II LRTs", {
  dat <- setup_poisson()
  mod <- glm(y ~ Sex + Stim, data = dat, family = poisson)
  tab <- bbmake_lrt_table(mod)
  expect_match(attr(tab, "note"), "Type II")
  expect_match(attr(tab, "note"), "likelihood-ratio")
})

# ==================================================================
# glmmTMB (the actual use case)
# ==================================================================

# ==================================================================
# keep: blocking / additive fixed effects stay in later LRTs
# ==================================================================

lrt_bind_manual <- function(models) {
  pieces <- lapply(models, lrt_drop1_manual)
  pieces <- pieces[vapply(pieces, nrow, integer(1)) > 0L]
  dplyr::bind_rows(pieces) |>
    dplyr::select(Term, Df, LRT, `Pr(>Chi)`) |>
    dplyr::rename(p = `Pr(>Chi)`)
}

lrt_drop_update <- function(model, terms) {
  f <- stats::as.formula(
    paste("~ .", paste(paste0("-", terms), collapse = " ")),
    env = environment(stats::formula(model))
  )
  cl <- stats::update(model, f, evaluate = FALSE)
  eval(cl, envir = environment(stats::formula(model)))
}

capture_lrt <- function(expr) {
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

test_that("interaction formulas do not message about additive terms", {
  dat <- setup_poisson()
  expect_silent(bbmake_lrt_table(glm(
    y ~ Sex * Stim, data = dat, family = poisson
  )))
  expect_silent(bbmake_lrt_table(glm(
    y ~ Sex * Stim * Diet, data = dat, family = poisson
  )))
  expect_silent(bbmake_lrt_table(glm(
    y ~ id + Sex + Stim, data = dat, family = poisson
  )))
})

test_that("an additive term is dropped with the interaction unless kept", {
  dat <- setup_poisson()
  dat$Block <- dat$id
  mod <- glm(y ~ Block + Sex * Stim, data = dat, family = poisson)

  offered <- lrt_drop1_manual(mod)
  reduced <- lrt_drop_update(mod, offered$Term)
  expected <- lrt_bind_manual(list(mod, reduced))

  out <- capture_lrt(bbmake_lrt_table(mod))

  expect_equal(out$tab$Term, expected$Term)
  expect_equal(out$tab$Df, as.integer(expected$Df))
  expect_equal(out$tab$LRT, expected$LRT)
  expect_equal(out$tab$p, expected$p)
  expect_match(
    out$message,
    "`Block` is removed together with `Sex:Stim`",
    fixed = TRUE
  )
  expect_match(out$message, "Later LRTs are fit without `Block`", fixed = TRUE)
  expect_match(out$message, 'Pass keep = "Block"', fixed = TRUE)
  expect_false(grepl("retained", attr(out$tab, "note"), fixed = TRUE))
})

test_that("keep retains a blocking factor in later LRTs and says so", {
  dat <- setup_poisson()
  dat$Block <- dat$id
  # Drop a Sex x Block corner so the block is not orthogonal to Sex.
  dat <- dat[!(dat$Block %in% c("1", "2") & dat$Sex == "F"), ]
  mod <- glm(y ~ Block + Sex * Stim, data = dat, family = poisson)

  offered <- lrt_drop1_manual(mod)
  drop_now <- setdiff(offered$Term, "Block")
  reduced <- lrt_drop_update(mod, drop_now)
  offered_reduced <- lrt_drop1_manual(reduced)
  offered_reduced <- offered_reduced[offered_reduced$Term != "Block", ]
  expected <- dplyr::bind_rows(
    dplyr::select(offered, Term, Df, LRT, `Pr(>Chi)`),
    dplyr::select(offered_reduced, Term, Df, LRT, `Pr(>Chi)`)
  ) |>
    dplyr::rename(p = `Pr(>Chi)`)

  dropped <- capture_lrt(bbmake_lrt_table(mod))
  out <- capture_lrt(bbmake_lrt_table(mod, keep = "Block"))

  expect_equal(out$tab$Term, expected$Term)
  expect_equal(out$tab$Df, as.integer(expected$Df))
  expect_equal(out$tab$LRT, expected$LRT)
  expect_equal(out$tab$p, expected$p)
  expect_false(any(duplicated(out$tab$Term)))
  expect_equal(
    out$tab$LRT[out$tab$Term == "Sex:Stim"],
    dropped$tab$LRT[dropped$tab$Term == "Sex:Stim"]
  )
  expect_false(isTRUE(all.equal(
    out$tab$LRT[out$tab$Term == "Sex"],
    dropped$tab$LRT[dropped$tab$Term == "Sex"]
  )))
  expect_match(out$message, "`Block` stays in every reduced model", fixed = TRUE)
  expect_match(out$message, "Later LRTs are fit with `Block`", fixed = TRUE)
  expect_false(grepl("is removed", out$message, fixed = TRUE))
  expect_match(attr(out$tab, "note"), "Block retained in every reduced model", fixed = TRUE)
  expect_match(attr(out$tab, "note"), "Type II", fixed = TRUE)
})

test_that("keep can retain one additive term while another is still dropped", {
  dat <- setup_poisson()
  dat$Block <- dat$id
  dat$Cohort <- dat$Diet
  dat <- dat[!(dat$Block %in% c("1", "2") & dat$Sex == "F"), ]
  mod <- glm(y ~ Block + Cohort + Sex * Stim, data = dat, family = poisson)

  offered <- lrt_drop1_manual(mod)
  reduced <- lrt_drop_update(mod, setdiff(offered$Term, "Block"))
  offered_reduced <- lrt_drop1_manual(reduced)
  offered_reduced <- offered_reduced[offered_reduced$Term != "Block", ]
  expected <- dplyr::bind_rows(
    dplyr::select(offered, Term, Df, LRT, `Pr(>Chi)`),
    dplyr::select(offered_reduced, Term, Df, LRT, `Pr(>Chi)`)
  ) |>
    dplyr::rename(p = `Pr(>Chi)`)

  out <- capture_lrt(bbmake_lrt_table(mod, keep = "Block"))

  expect_equal(out$tab$Term, expected$Term)
  expect_equal(out$tab$LRT, expected$LRT)
  expect_match(out$message, "`Block` stays in every reduced model", fixed = TRUE)
  expect_match(
    out$message,
    "`Cohort` is removed together with `Sex:Stim`",
    fixed = TRUE
  )
  expect_match(out$message, "Add `Cohort` to `keep`", fixed = TRUE)
  expect_match(out$message, "Later LRTs are fit without `Cohort`", fixed = TRUE)
})

test_that("several kept terms stay in and are not described as removed", {
  dat <- setup_poisson()
  dat$Block <- dat$id
  dat$Cohort <- dat$Diet
  dat <- dat[!(dat$Block %in% c("1", "2") & dat$Sex == "F"), ]
  mod <- glm(y ~ Block + Cohort + Sex * Stim, data = dat, family = poisson)

  offered <- lrt_drop1_manual(mod)
  reduced <- lrt_drop_update(mod, setdiff(offered$Term, c("Block", "Cohort")))
  offered_reduced <- lrt_drop1_manual(reduced)
  offered_reduced <- offered_reduced[
    !offered_reduced$Term %in% c("Block", "Cohort"),
  ]
  expected <- dplyr::bind_rows(
    dplyr::select(offered, Term, Df, LRT, `Pr(>Chi)`),
    dplyr::select(offered_reduced, Term, Df, LRT, `Pr(>Chi)`)
  ) |>
    dplyr::rename(p = `Pr(>Chi)`)

  out <- capture_lrt(bbmake_lrt_table(mod, keep = c("Block", "Cohort")))

  expect_equal(out$tab$Term, expected$Term)
  expect_equal(out$tab$LRT, expected$LRT)
  expect_false(any(duplicated(out$tab$Term)))
  expect_match(
    out$message,
    "`Block` and `Cohort` stay in every reduced model",
    fixed = TRUE
  )
  expect_match(
    out$message,
    "Later LRTs are fit with `Block` and `Cohort`",
    fixed = TRUE
  )
  expect_false(grepl("removed", out$message, fixed = TRUE))
  expect_match(
    attr(out$tab, "note"),
    "Block, Cohort retained in every reduced model",
    fixed = TRUE
  )
})

test_that("keeping a last-stage main effect does not repeat or change tests", {
  dat <- setup_poisson()
  mod <- glm(y ~ Sex * Stim, data = dat, family = poisson)
  base <- bbmake_lrt_table(mod)
  out <- capture_lrt(bbmake_lrt_table(mod, keep = "Sex"))

  expect_equal(out$tab$Term, base$Term)
  expect_equal(out$tab$Df, base$Df)
  expect_equal(out$tab$LRT, base$LRT)
  expect_equal(out$tab$p, base$p)
  expect_match(out$message, "`Sex` stays in every reduced model", fixed = TRUE)
  expect_match(attr(out$tab, "note"), "Sex retained in every reduced model", fixed = TRUE)
})

test_that("keeping an interaction leaves its margins untested", {
  dat <- setup_poisson()
  mod <- glm(y ~ Sex * Stim, data = dat, family = poisson)
  manual <- lrt_drop1_manual(mod)
  out <- capture_lrt(bbmake_lrt_table(mod, keep = "Sex:Stim"))

  expect_equal(out$tab$Term, manual$Term)
  expect_equal(out$tab$LRT, manual$LRT)
  expect_false(any(c("Sex", "Stim") %in% out$tab$Term))
  expect_match(out$message, "`Sex:Stim` stays in every reduced model", fixed = TRUE)
})

test_that("unknown or invalid keep fails before fitting reduced models", {
  dat <- setup_poisson()
  mod <- glm(y ~ Sex * Stim, data = dat, family = poisson)

  expect_error(bbmake_lrt_table(mod, keep = "Block"), "not in the model")
  expect_error(bbmake_lrt_table(mod, keep = "Block"), "Sex:Stim", fixed = TRUE)
  expect_error(bbmake_lrt_table(mod, keep = c("Sex", "Block")), "Block")
  expect_error(bbmake_lrt_table(mod, keep = 1), "`keep` must be")
  expect_error(bbmake_lrt_table(mod, keep = NA_character_), "`keep` must be")
  expect_error(bbmake_lrt_table(mod, keep = ""), "`keep` must be")
})

test_that("glmmTMB two-way with random intercept matches manual pipeline", {
  suppressWarnings(skip_if_not_installed("glmmTMB"))

  set.seed(1)
  dat <- expand.grid(
    Sex  = factor(c("F", "M")),
    Stim = factor(c("A", "B")),
    id   = factor(1:15)
  )
  dat$yb <- pmin(pmax(rbeta(nrow(dat), 2, 3), 1e-3), 1 - 1e-3)

  mod2 <- suppressWarnings(glmmTMB::glmmTMB(
    yb ~ Sex * Stim + (1 | id),
    data = dat,
    family = glmmTMB::beta_family()
  ))
  mod1 <- suppressWarnings(update(mod2, . ~ . - Sex:Stim))

  expected <- suppressWarnings(dplyr::bind_rows(
    lrt_drop1_manual(mod2),
    lrt_drop1_manual(mod1)
  ))

  tab <- suppressWarnings(bbmake_lrt_table(mod2))

  expect_equal(tab$Term, expected$Term)
  expect_equal(tab$LRT, expected$LRT)
  expect_equal(tab$p, expected$`Pr(>Chi)`)
  expect_false(any(grepl("|", tab$Term, fixed = TRUE)))
  expect_false(any(tab$Term %in% c("1", "(Intercept)")))
  expect_equal(nrow(tab), 3L)
})
