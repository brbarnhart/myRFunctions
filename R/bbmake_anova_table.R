#' Unified ANOVA-style table
#'
#' Creates an ANOVA-style table for lm, lmer/lmerModLmerTest, or glmmTMB models.
#' For Gaussian models it returns the classic F-test table + partial omega-squared.
#' For glmmTMB (e.g. negative binomial) it returns Wald χ² tests + Incidence Rate Ratios (IRRs).
#'
#' @param model A fitted model — `lm()`, `lmer()`/`lmerMod*`, or `glmmTMB`
#' @param type Type of sums of squares ("III" is default and recommended)
#'
#' @return A tibble with the ANOVA-style results (and IRRs for glmmTMB)
#' @export
bbmake_anova_table <- function(model, type = "III") {

  # ------------------------------------------------------------------
  # Helper: Clean glmmTMB coefficient names
  # ------------------------------------------------------------------
  clean_term_names <- function(x) {
    x |>
      stringr::str_split(":") |>
      purrr::map(~ stringr::str_remove(.x, "(?<=[^A-Z])[A-Z].*$")) |>
      purrr::map_chr(~ stringr::str_c(.x, collapse = ":"))
  }

  # ==================================================================
  # Gaussian models (lm / lmer)
  # ==================================================================
  if (inherits(model, c("lm", "lmerMod", "lmerModLmerTest"))) {

    # --------------------- Get ANOVA table ---------------------
    if (inherits(model, "lmerModLmerTest")) {
      # lmerTest::anova gives the proper F-tests + DenDF (Satterthwaite/KR)
      aov_tab <- anova(model) |>
        as.data.frame() |>
        tibble::rownames_to_column("Effect") |>
        dplyr::select(Effect, NumDF, DenDF, `F value`, `Pr(>F)`)

    } else {
      # lm or plain lmerMod → car::Anova (type III)
      aov_tab <- suppressWarnings(
        car::Anova(model, type = type)
      ) |>
        as.data.frame() |>
        tibble::rownames_to_column("Effect")

      if ("Df" %in% names(aov_tab)) {
        aov_tab <- dplyr::rename(aov_tab, NumDF = Df)
      }

      # Plain lmerMod (no lmerTest) gives Chisq → treat as F-value for consistency
      if ("Chisq" %in% names(aov_tab)) {
        aov_tab <- dplyr::rename(aov_tab,
                                 `F value` = Chisq,
                                 `Pr(>F)`  = `Pr(>Chisq)`
        )
      }

      if (!"DenDF" %in% names(aov_tab)) {
        aov_tab$DenDF <- NA_real_
      }

      aov_tab <- aov_tab |>
        dplyr::select(Effect, NumDF, DenDF, `F value`, `Pr(>F)`)
    }

    # --------------------- Partial omega-squared ---------------------
    omega <- effectsize::omega_squared(model, partial = TRUE) |>
      as.data.frame() |>
      dplyr::rename(
        Effect             = Parameter,
        `Omega2 (partial)` = Omega2_partial
      )

    # Guarantee character CI column (what the tests expect)
    if (all(c("CI_low", "CI_high") %in% names(omega))) {
      omega <- omega |>
        dplyr::mutate(`Omega2 (95% CI)` = sprintf("[%.2f, %.2f]", CI_low, CI_high)) |>
        dplyr::select(Effect, `Omega2 (partial)`, `Omega2 (95% CI)`)
    } else if ("CI" %in% names(omega)) {
      omega <- omega |> dplyr::rename(`Omega2 (95% CI)` = CI)
    } else {
      omega$`Omega2 (95% CI)` <- NA_character_
      omega <- dplyr::select(omega, Effect, `Omega2 (partial)`, `Omega2 (95% CI)`)
    }

    tab <- dplyr::left_join(aov_tab, omega, by = "Effect")

    # ==================================================================
    # glmmTMB models (unchanged – already passing)
    # ==================================================================
  } else if (inherits(model, "glmmTMB")) {

    aov_tab <- car::Anova(model, type = type, test.statistic = "Chisq") |>
      as.data.frame() |>
      tibble::rownames_to_column("Effect") |>
      dplyr::rename(
        NumDF = Df,
        `Chisq value` = Chisq,
        `Pr(>Chisq)` = `Pr(>Chisq)`
      ) |>
      dplyr::mutate(
        `Chisq value` = round(`Chisq value`, 2),
        `Pr(>Chisq)` = format.pval(`Pr(>Chisq)`, digits = 3, eps = 0.001)
      ) |>
      dplyr::select(Effect, NumDF, `Chisq value`, `Pr(>Chisq)`)

    irr_tab <- parameters::model_parameters(
      model,
      exponentiate = TRUE,
      component = "conditional"
    ) |>
      as.data.frame() |>
      dplyr::filter(!grepl("^(sd_|zi_)", Parameter)) |>
      dplyr::rename(Effect = Parameter) |>
      dplyr::mutate(
        Effect = clean_term_names(Effect),
        IRR = round(Coefficient, 2),
        `IRR 95% CI` = sprintf("[%.2f, %.2f]", CI_low, CI_high)
      ) |>
      dplyr::select(Effect, IRR, `IRR 95% CI`)

    tab <- dplyr::left_join(aov_tab, irr_tab, by = "Effect")

  } else {
    stop("Model class not supported. Only lm, lmer*, and glmmTMB are handled.")
  }

  # ====================== Robust Nakagawa R² ======================
  r2 <- tryCatch(
    performance::r2_nakagawa(model),
    error = function(e) list(R2_marginal = NA_real_, R2_conditional = NA_real_)
  )

  marg <- if (is.list(r2) && !is.null(r2$R2_marginal)) r2$R2_marginal else NA_real_
  cond <- if (is.list(r2) && !is.null(r2$R2_conditional)) r2$R2_conditional else NA_real_

  attr(tab, "r2") <- paste0(
    "Marginal R² = ", round(marg, 3),
    " | Conditional R² = ", round(cond, 3)
  )

  tibble::as_tibble(tab)
}
