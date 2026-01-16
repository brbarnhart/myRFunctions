#' Create APA-style mixed-effects ANOVA table with partial omega-squared and CIs
#'
#' @param model A fitted lmer model (from lmerTest)
#' @param title Character string for the table title (default: "Mixed-Effects ANOVA")
#'
#' @return A flextable object (from nice_table)
#' @export
#'
#' @examples
#' stat_table <- make_anova_table(model, title = "mPFC Delta AUC Analysis")
make_anova_table <- function(model, title = "Mixed-Effects ANOVA") {

  aov_tab <- anova(model) |>
    as.data.frame() |>
    tibble::rownames_to_column(var = "Effect")

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

  combined_table <- dplyr::left_join(aov_tab, omega_tab, by = "Effect") |>
    dplyr::select(Effect, NumDF, DenDF, `F value`, `Pr(>F)`,
                  `Omega2 (partial)`, `Omega2 (95% CI)`)

  rempsyc::nice_table(
    combined_table,
    col.format.p = 5,
    title = c("", title),
    note = c(" * p < .05, ** p < .01, *** p < .001")
    ) |>
    flextable::autofit()
}
