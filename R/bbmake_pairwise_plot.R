#' Create pairwise comparison plot (now handles both Gaussian and GLMM)
#'
#' Automatically shows the correct effect-size annotation in the top-right corner.
#'
#' @param emm emmGrid object from emmeans
#' @param pw_table Output from bbmake_pairwise_table()
#' @param formula Formula passed to emmip()
#' @param aesthetic Default aes(x = xvar, y = yvar)
#' @param y.adjust Vertical nudge for significance stars
#'
#' @return ggplot2 object
#' @export
bbmake_pairwise_plot <- function(
    emm,
    pw_table,
    formula,
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

  emm_df <- emm_df |>
    rename(
      emmean = all_of(mean_col),
      lower  = all_of(lower_col),
      upper  = all_of(upper_col)
    )

  # Significance markers (works for any number of contrasts)
  pw_contrasts <- pw_table |>
    as.data.frame() |>
    mutate(
      group1    = str_trim(str_split_i(contrast, " - ", 1)),
      group2    = str_trim(str_split_i(contrast, " - ", 2)),
      p.signif  = case_when(
        p.value < 0.001 ~ "***",
        p.value < 0.01  ~ "**",
        p.value < 0.05  ~ "*",
        p.value < 0.1   ~ ".",
        TRUE            ~ "ns"
      ),
      y.position = max(emm_df$upper, na.rm = TRUE) + y.adjust
    )

  plot_data <- emmip(emm, formula, CIs = TRUE, plotit = FALSE)

  # Dynamic top-right annotation (only first contrast's value is shown – cleanest look)
  if ("Cohen's d" %in% colnames(pw_table)) {
    annot <- sprintf("Cohen's d = %.2f", abs(pw_table$`Cohen's d`[1]))
  } else if ("Rate Ratio" %in% colnames(pw_table)) {
    annot <- sprintf("Rate Ratio = %.2f\n%s",
                     pw_table$`Rate Ratio`[1],
                     pw_table$`RR 95% CI`[1])
  } else if ("Odds Ratio" %in% colnames(pw_table)) {
    annot <- sprintf("Odds Ratio = %.2f\n%s",
                     pw_table$`Odds Ratio`[1],
                     pw_table$`OR 95% CI`[1])
  } else if ("Estimate" %in% colnames(pw_table) || "Mean Difference" %in% colnames(pw_table)) {
    annot <- sprintf("Difference = %.2f", pw_table[[which(colnames(pw_table) %in% c("Estimate", "Mean Difference"))[1]]][1])
  } else {
    annot <- ""
  }

  ggplot(plot_data, aesthetic) +
    geom_hline(yintercept = 0, linetype = "solid", color = "gray50", linewidth = 0.6) +
    geom_line(linewidth = 0.8) +
    geom_point(size = 4.5) +
    geom_errorbar(aes(ymin = LCL, ymax = UCL), width = 0.25) +
    theme_get() +
    annotate("text",
             x = Inf, y = Inf,
             label = annot,
             hjust = 1.05, vjust = 2.2,
             size = 4.2, lineheight = 0.9) +
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
