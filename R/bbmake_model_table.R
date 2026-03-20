#' Unified ANOVA-style table
#'
#' Creates an ANOVA-style table for lm, lmer/lmerModLmerTest, or glmmTMB models.
#' For Gaussian models it returns the classic F-test table + partial omega-squared.
#' For glmmTMB (e.g. negative binomial) it returns Wald-z tests + Incidence Rate Ratios (IRRs).
#'
#' @param model A fitted model — `lm()`, `lmer()`/`lmerMod*`, or `glmmTMB`
#' @param type Type of sums of squares ("III" is default and recommended)
#'
#' @return A tibble with the ANOVA-style results (and IRRs for glmmTMB)
#' @export
bbmake_model_table <- function(model, type = NULL) {

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
        dplyr::rename(
          `F` = `F value`,
          p = `Pr(>F)`
        ) |>
        dplyr::select(Effect, NumDF, DenDF, `F`, p) |>
        dplyr:: filter(!Effect %in% c("(Intercept)", "Residuals"))

    } else {
      # lm or plain lmerMod → car::Anova (type II)
      if (is.null(type)) {type <- "II"}
      aov_tab <- suppressWarnings(
        car::Anova(model, type = type, test.statistic = "F")
      ) |>
        as.data.frame() |>
        tibble::rownames_to_column("Effect") |>
        rename(
          `F` = `F value`,
          p = `Pr(>F)`
        )

      aov_tab <- aov_tab |>
        dplyr::select(Effect, any_of(c("DF", "NumDF", "DenDF")), `F`, p) |>
        dplyr:: filter(!Effect %in% c("(Intercept)", "Residuals"))
    }

    # --------------------- Partial omega-squared ---------------------
    omega <- effectsize::omega_squared(model, partial = TRUE) |>
      as.data.frame() |>
      dplyr::rename(
        Effect             = Parameter,
        `Omega2 (partial)` = Omega2_partial
      ) |>
      select(Effect, `Omega2 (partial)`, CI_low, CI_high)

    tab <- dplyr::left_join(aov_tab, omega, by = "Effect")

    # ==================================================================
    # glmmTMB models
    # ==================================================================
  } else if (inherits(model, "glmmTMB")) {

    if (!is.null(type)) {

      # ── Term-wise Wald χ² table (ANOVA-style) ────────────────────────────────

      # Make sure car is available
      if (!requireNamespace("car", quietly = TRUE)) {
        stop("Package 'car' is required for glmm_anova_type. Please install it.")
      }

      # Run Anova – type can be "II", "III", 2 or 3 (car accepts both)
      aov_tab <- car::Anova(model,
                            type      = type,
                            component = "cond",          # main (mean) model
                            test.statistic = "Chisq")

      tab <- aov_tab |>
        as.data.frame() |>
        tibble::rownames_to_column("Effect") |>
        dplyr::as_tibble() |>
        dplyr::rename(
          Df       = Df,
          `Wald χ²` = Chisq,
          `p`      = `Pr(>Chisq)`
        ) |>
        dplyr::mutate(Effect = clean_term_names(Effect)) |>
        dplyr::select(Effect, Df, `Wald χ²`, p)

      # Add on IRR information for a pseudo effect size
      IRR_tab <- parameters::model_parameters(model, exponentiate = TRUE) |>
        dplyr::as_tibble() |>
        dplyr::rename(
          Effect     = Parameter,
          IRR        = Coefficient
        ) |>
        tidyr::drop_na(z) |>
        dplyr::filter(Effect != "(Intercept)") |>
        dplyr::mutate(Effect = clean_term_names(Effect)) |>
        dplyr::select(Effect, IRR, CI_low, CI_high)

      tab <- left_join(tab, IRR_tab, by = "Effect")

      # Optional: add note about test type
      attr(tab, "note") <- sprintf("Type %s Wald χ² tests (conditional component)",
                                   if (is.numeric(type)) type else toupper(type))

    } else {

      # ── Original per-coefficient style with IRRs (default) ────────────────────

      tab <- parameters::model_parameters(model, exponentiate = TRUE) |>
        dplyr::as_tibble() |>
        dplyr::rename(
          Effect     = Parameter,
          IRR        = Coefficient
        ) |>
        tidyr::drop_na(z) |>
        dplyr::filter(Effect != "(Intercept)") |>
        dplyr::mutate(Effect = clean_term_names(Effect)) |>
        dplyr::select(Effect, z, p, IRR, CI_low, CI_high)

    }

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
