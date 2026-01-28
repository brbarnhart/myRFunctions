bbmake_anova_table <- function(model) {

  # Get ANOVA table
  aov_obj <- anova(model)
  aov_tab <- as.data.frame(aov_obj) |>
    tibble::rownames_to_column(var = "Effect")

  # Normalize column names (lm vs lmerTest/lme4 differences)
  names(aov_tab) <- gsub("^Df$", "NumDF", names(aov_tab))          # lm → Df → NumDF
  names(aov_tab) <- gsub("^NumDF$", "NumDF", names(aov_tab))       # already correct
  names(aov_tab) <- gsub("^denDF$|^DenDF$", "DenDF", names(aov_tab), ignore.case = TRUE)

  # Ensure DenDF exists (lm doesn't have it)
  if (!"DenDF" %in% names(aov_tab)) {
    aov_tab$DenDF <- NA_real_
  }

  # Make sure NumDF exists (should after rename)
  if (!"NumDF" %in% names(aov_tab)) {
    stop("ANOVA table does not contain a numeric DF column (tried 'Df' / 'NumDF')")
  }

  # Compute partial omega-squared + CIs
  omega_tab <- effectsize::omega_squared(model, partial = TRUE) |>
    as.data.frame() |>
    dplyr::rename(Effect = Parameter) |>
    dplyr::mutate(
      `Omega2 (partial)` = round(Omega2_partial, 2),
      `Omega2 (95% CI)` = dplyr::if_else(
        is.na(CI_low) | is.na(CI_high),
        "—",
        sprintf("[%.2f, %.2f]", CI_low, CI_high)
      )
    ) |>
    dplyr::select(Effect, `Omega2 (partial)`, `Omega2 (95% CI)`)

  # Join and select (now safe)
  combined_table <- dplyr::left_join(aov_tab, omega_tab, by = "Effect") |>
    dplyr::select(
      Effect,
      NumDF,
      DenDF,
      `F value`,
      `Pr(>F)`,
      `Omega2 (partial)`,
      `Omega2 (95% CI)`
    ) |>
    # Optional: return as tibble for consistency
    tibble::as_tibble()

  # Format p-values (safe even if some are NA)
  if ("Pr(>F)" %in% names(combined_table)) {
    combined_table$`Pr(>F)` <- format.pval(combined_table$`Pr(>F)`, digits = 3, eps = 0.001)
  }

  combined_table
}
