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

  aov_tab <- lmerTest::anova(model) |>
    as.data.frame() |>
    tibble::rownames_to_column(var = "Effect") |>
    dplyr::mutate(
      Sig. = dplyr::case_when(...),
      `Pr(>F)` = stats::format.pval(`Pr(>F)`, digits = 3, eps = .001)
    )

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
    dplyr::select(...)

  rempsyc::nice_table(
    combined_table,
    title = c("", title),
    note = c("Note. . p < .1, * p < .05, ** p < .01, *** p < .001")
    ) |>
    flextable::autofit()
}
