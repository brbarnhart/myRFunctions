#' Create APA-style pairwise comparison table
#'
#' Works with raw `emmGrid` OR the output of `bbmake_pairwise_table()`.
#' Automatically detects the p-value column (works for both Gaussian and GLMM tables).
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

  # Dynamically find the p-value column (robust to different table widths)
  p_col <- which(names(pw_table) %in% c("p.value", "p", "Pr(>|z|)", "Pr(>|t|)"))[1]

  # Fallback: any column whose name contains "p." or ends with ".p"
  if (is.na(p_col) || length(p_col) == 0) {
    p_col <- which(grepl("p\\.|\\.p$|^p$", names(pw_table), ignore.case = TRUE))[1]
  }

  # If still not found, let nice_table guess (or set to NULL)
  if (is.na(p_col) || length(p_col) == 0) {
    p_col <- NULL
  }

  rempsyc::nice_table(
    pw_table,
    col.format.p = p_col,          # ← now dynamic!
    title = c("", title),
    note = c(" * p < .05, ** p < .01, *** p < .001")
  ) |>
    flextable::autofit()
}
