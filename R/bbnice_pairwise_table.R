#' Create APA-style pairwise comparison table
#'
#' Works with raw `emmGrid` OR the output of `bbmake_pairwise_table()`.
#'
#' @param pw_table `emmGrid` or tibble from `bbmake_pairwise_table()`
#' @param title Table title
#'
#' @return A flextable object
#' @export
bbnice_pairwise_table <- function(pw_table, title = "Pairwise Comparisons Table") {

  # If user passed raw emmGrid, auto-process it
  if (inherits(pw_table, "emmGrid")) {
    pw_table <- bbmake_pairwise_table(pw_table)
  }

  rempsyc::nice_table(
    pw_table,
    col.format.p = 5,
    title = c("", title),
    note = c(" * p < .05, ** p < .01, *** p < .001")
  ) |>
    flextable::autofit()
}
