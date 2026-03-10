#' Create pairwise comparison plot (now handles both Gaussian and GLMM)
#'
#' Automatically shows the correct effect-size annotation in the top-right corner.
#'
#' @param emm emmGrid object from emmeans
#' @param pw_table Output from bbmake_pairwise_table()
#' @param formula Formula passed to emmip()
#' @param grouping_vars variables to group by for significance markers
#' @param aesthetic Default aes(x = xvar, y = yvar)
#' @param y.adjust Vertical nudge for significance stars
#'
#' @return ggplot2 object
#' @export
bbmake_pairwise_plot <- function(
    emm,
    pw_table,
    formula,
    grouping_vars = NULL,
    aesthetic = aes(x = xvar, y = yvar),
    y.adjust = 0
) {

  emm_df <- as.data.frame(emm)

  # Detect the mean column (emmean on link scale, response on response scale)
  if ("response" %in% names(emm_df)) {
    mean_col <- "response"
  } else if ("emmean" %in% names(emm_df)) {
    mean_col <- "emmean"
  } else {
    stop("No mean column (emmean or response) found in as.data.frame(emm)")
  }

  # Also handle lower.CL vs asymp.LCL (glmmTMB often uses asymp.*)
  lower_col <- ifelse("lower.CL" %in% names(emm_df), "lower.CL", "asymp.LCL")
  upper_col <- ifelse("upper.CL" %in% names(emm_df), "upper.CL", "asymp.UCL")

  emm_df <- as.data.frame(emm) |>
    rename(
      emmean = all_of(mean_col),
      lower  = all_of(lower_col),
      upper  = all_of(upper_col)
    )

  if (is.null(grouping_vars) || length(grouping_vars) == 0) {
    # No grouping → single global max
    emm_df <- emm_df |>
      mutate(
        y.position = max(upper, na.rm = TRUE),
        y.position = y.position + y.adjust
      )
  } else {
    emm_df <- emm_df |>
      group_by(across(all_of(grouping_vars))) |>
      mutate(
        y.position = max(upper, na.rm = TRUE),
        y.position = y.position + y.adjust) |>
      ungroup()
  }

  # Significance markers (works for any number of contrasts)
  pw_contrasts <- pw_table |>
    as.data.frame() |>
    mutate(
      p.signif  = case_when(
        p.value < 0.001 ~ "***",
        p.value < 0.01  ~ "**",
        p.value < 0.05  ~ "*",
        p.value < 0.1   ~ ".",
        TRUE            ~ "ns"
      )
    )

  if (attr(pw_table, "estName") == "estimate") {
    pw_contrasts <- pw_contrasts |>
      mutate(
        group1    = str_trim(str_split_i(contrast, " - ", 1)),
        group2    = str_trim(str_split_i(contrast, " - ", 2))
      )
  } else if (attr(pw_table, "estName") == "ratio") {
    pw_contrasts <- pw_contrasts |>
      mutate(
        group1    = str_trim(str_split_i(contrast, " / ", 1)),
        group2    = str_trim(str_split_i(contrast, " / ", 2))
      )
  }

  if (!is.null(grouping_vars) && length(grouping_vars) > 0) {
    y_pos_lookup <- emm_df |>
      select(all_of(grouping_vars), y.position) |>
      distinct()

    pw_contrasts <- pw_contrasts |>
      left_join(y_pos_lookup, by = grouping_vars)
  } else {
    # No grouping → just one value
    pw_contrasts <- pw_contrasts |>
      mutate(
        y.position = max(emm_df$upper, na.rm = TRUE),
        y.position = y.position + y.adjust
      )
  }

  plot_data <- emmip(emm, formula, CIs = TRUE, plotit = FALSE)

  # Dynamic top-right annotation (only first contrast's value is shown – cleanest look)
  if ("Cohen's d" %in% colnames(pw_table)) {
    pw_table <- pw_table |> mutate (
      effect_annotation = sprintf("Cohen's d = %.2f", abs(`Cohen's d`))
    )
  } else if ("IRR" %in% colnames(pw_table)) {
    pw_table <- pw_table |> mutate (
      effect_annotation = sprintf("IRR = %.2f", IRR)
    )
  } else if ("Odds Ratio" %in% colnames(pw_table)) {
    pw_table <- pw_table |> mutate (
      effect_annotation = sprintf("Odds Ratio = %.2f", `Odds Ratio`)
    )
  } else {
    pw_table <- pw_table |> mutate (
      effect_annotation = ""
    )
  }

  ggplot(plot_data, aesthetic) +
    geom_hline(yintercept = 0, linetype = "solid", color = "gray50", linewidth = 0.6) +
    geom_line(linewidth = 0.8) +
    geom_point(size = 4.5) +
    geom_errorbar(aes(ymin = LCL, ymax = UCL), width = 0.25) +
    theme_get() +
    geom_text(
      data = pw_table,
      aes(label = effect_annotation),
      x = Inf, y = Inf,
      hjust = 1.05, vjust = 2.2,
      size = 4.2, inherit.aes = FALSE
    ) +
    stat_pvalue_manual(
      data = pw_contrasts,
      label = "p.signif",
      xmin = "group1", xmax = "group2",
      y.position = "y.position",
      hide.ns = FALSE,
      tip.length = 0.01
    ) +
    scale_y_continuous(expand = expansion(mult = c(0.0, 0.18)))
}
