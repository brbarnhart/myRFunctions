#' Create pairwise comparison table with Cohen's d and CIs
#'
#' @param emm An emmGrid object (output from emmeans())
#' @param pw_table A table_df object (output from bbmake_pairwise_table())
#' @param formula Formula object used to create emm
#' @param group1 "String" First group for significance markers for pairwise comparisons
#' @param group2 "String" Second group for significance markers for pairwise comparisons
#' @param y.adjust Float or list of floats representing the y adjustment for each significance marker
#'
#' @return A ggplot2 plot object
#' @export
bbmake_pairwise_plot <- function(emm, pw_table, formula, group1, group2, y.adjust = 0) {
  # Convert to data frame for plotting
  emm_df <- as.data.frame(emm) |>
    rename(emmean = emmean,
           lower = lower.CL,
           upper = upper.CL)

  # Prepare significance labels (diet differences within each sex)
  pw_contrasts <- pw_table |>
    as.data.frame() |>
    mutate(
      group1 = group1,
      group2 = group2,  # order matches your output
      p.signif = case_when(
        p.value < 0.001 ~ "***",
        p.value < 0.01  ~ "**",
        p.value < 0.05  ~ "*",
        p.value < 0.1   ~ ".",
        TRUE            ~ "ns"
      ),
      y.position = max(emm_df$upper) + y.adjust
    )

  plot_data <- emmip(emm,
                     formula,
                     CIs = TRUE,
                     plotit = FALSE)

  p <- plot_data |>
    ggplot(aes(x = xvar, y = yvar)) +
    geom_hline(yintercept = 0, linetype = "solid", color = "gray50", linewidth = 0.6) +
    geom_line() +
    geom_point(size = 4) +
    geom_errorbar(aes(ymin = LCL, ymax = UCL), width = 0.2) +
    theme_get() +
    geom_text(
      data = pw_contrasts,
      aes(label = paste("Cohen's d: ", sprintf("%.2f", abs(`Cohen's d`)))),
      x = Inf, y = Inf,
      hjust = 1.05, vjust = 1.5, size = 4, inherit.aes = FALSE
    ) +
    stat_pvalue_manual(data = pw_contrasts,
                       label = "p.signif",
                       xmin = "group1", xmax = "group2",
                       y.position = "y.position",
                       hide.ns = FALSE) +
    scale_y_continuous(expand = expansion(mult = c(0.0, 0.12)))

  p
}



