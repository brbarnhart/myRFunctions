#' Create pairwise comparison table with Cohen's d and CIs
#'
#' @param pw An emmGrid object (output from pairs() from emmeans)
#'
#' @return A data.frame object
#' @export
#'
#' @examples
#' pw_table <- bbmake_pairwise_table(pw)
bbmake_pairwise_table <- function(pw) {

  pw_summary <- pw |>
    summary() |>
    as_tibble() |>
    rename(`mean diference` = estimate) |>
    select(contrast:df, t.ratio, p.value) |>
    rowid_to_column()

  cohen_d <- eff_size(
    pw,
    sigma = sigma(model),
    edf = df.residual(model),
    method = "identity"
  ) |>
    summary(infer = TRUE) |>
    mutate(contrast = str_remove_all(contrast, "[()]")) |>
    select(effect.size, lower.CL, upper.CL) |>
    rename(`Cohen's d` = effect.size) |>
    rowid_to_column()

  pw_table <- left_join(pw_summary, cohen_d, by = "rowid") |>
    select(!rowid)

  pw_table

}
