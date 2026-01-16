#' Create APA-style pairwise comparison table with Cohen's d and CIs
#'
#' @param pw An emmGrid object (output from pairs() from emmeans)
#' @param title Character string for the table title (default: "Pairwise Comparisons Table")
#'
#' @return A data.frame object
#' @export
#'
#' @examples
#' pw_table <- bbnice_pairwise_table(pw)
bbnice_pairwise_table <- function(pw, title = "Pairwise Comparisons Table") {

  rempsyc::nice_table(
    pw,
    title = c("", title),
    note = c(" * p < .05, ** p < .01, *** p < .001")
  ) |>
    flextable::autofit()

}
