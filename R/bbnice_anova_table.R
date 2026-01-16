#' Create APA-style mixed-effects ANOVA table with partial omega-squared and CIs
#'
#' @param anova_table the anova table produced by make_anova_table()
#' @param title Character string for the table title (default: "ANOVA table")
#'
#' @return A flextable object (from nice_table)
#' @export
#'
#' @examples
#' stat_table <- bbnice_anova_table(anova_table, title = "mPFC Delta AUC Analysis")
bbnice_anova_table <- function(anova_table, title = "ANOVA Table") {

  rempsyc::nice_table(
    anova_table,
    col.format.p = 5,
    title = c("", title),
    note = c(" * p < .05, ** p < .01, *** p < .001")
    ) |>
    flextable::autofit()

}
