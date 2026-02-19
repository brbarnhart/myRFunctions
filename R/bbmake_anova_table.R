#' Unified ANOVA-style table
#'
#' Creates an ANOVA-style table for lm, lmer/lmerModLmerTest, or glmmTMB models.
#' For Gaussian models it returns the classic F-test table + partial omega-squared.
#' For glmmTMB (e.g. negative binomial) it returns Wald χ² tests + Incidence Rate Ratios (IRRs).
#'
#' @param model A fitted model — `lm()`, `lmer()`/`lmerMod*`, or `glmmTMB`
#' @param type For glmmTMB only: Type of sums of squares ("III" is default and recommended)
#'
#' @return A tibble with the ANOVA-style results (and IRRs for glmmTMB)
#' @export
#'
#' @examples
#' # Gaussian models (your old workflow)
#' bbmake_anova_table(mod_lm)
#' bbmake_anova_table(mod_lmer)
#'
#' # Negative binomial GLMM
#' bbmake_anova_table(model_nb)
bbmake_anova_table <- function(model, type = "III") {

  # ------------------------------------------------------------------
  # Helper: Clean glmmTMB coefficient names (SexMale → Sex, etc.)
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

    tab <- parameters::model_parameters(
      model,
      effects = "fixed",
      ci_method = "satterthwaite"
    ) |>
      tibble::as_tibble() |>
      dplyr::rename(
        Effect = Parameter,
        `F value` = F,
        `Pr(>F)` = p
      ) |>
      dplyr::select(Effect, `F value`, `Pr(>F)`, dplyr::everything())

    # Partial omega-squared (your favourite)
    omega <- effectsize::omega_squared(model, partial = TRUE) |>
      as.data.frame() |>
      dplyr::rename(
        Effect = Parameter,
        `Omega2 (partial)` = Omega2_partial,
        `Omega2 (95% CI)` = CI
      )

    tab <- dplyr::left_join(tab, omega, by = "Effect")

    # ==================================================================
    # glmmTMB models (negative binomial, Poisson, etc.)
    # ==================================================================
  } else if (inherits(model, "glmmTMB")) {

    # Wald χ² tests (Type III)
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

    # Incidence Rate Ratios (cleaned names)
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

  # Add overall model fit (Nakagawa R²)
  r2 <- performance::r2_nakagawa(model)
  attr(tab, "r2") <- paste0(
    "Marginal R² = ", round(r2$R2_marginal, 3),
    " | Conditional R² = ", round(r2$R2_conditional, 3)
  )

  tab
}
