# Purpose: generate the figure showing LR perturbation-profile stability
# across GTEx normal-reference sizes for the study cohort.
# Input: an in-memory reference_size_resample_diagnostics data frame, or
#        results/06_Assess_SSN_Stability/ssn_reference_size_stability.qs2.
# Outputs: figures/06_Assess_SSN_Stability/reference_size_stability.{png,pdf}.
#
# Run from the repository root with Rscript scripts/06_Plot_SSN_Stability.R

library(qs2)
library(ggplot2)

if (!exists("reference_size_resample_diagnostics")) {
  reference_size_stability <- qs2::qs_read(
    "results/06_Assess_SSN_Stability/ssn_reference_size_stability.qs2"
  )
  reference_size_resample_diagnostics <- reference_size_stability$resample_diagnostics
}

dir.create("figures/06_Assess_SSN_Stability", recursive = TRUE, showWarnings = FALSE)

reference_size_levels <- c("50", "100", "250", "500", "898")
reference_size_resample_diagnostics$reference_size_label <- factor(
  reference_size_resample_diagnostics$reference_size,
  levels = as.integer(reference_size_levels),
  labels = reference_size_levels
)

reference_size_fill_colors <- c(
  "50" = "#D8E5F0",
  "100" = "#C2D7E8",
  "250" = "#A9C7DE",
  "500" = "#84ADD1",
  "898" = "#E39D98"
)
reference_size_line_colors <- c(
  "50" = "#6F99BC",
  "100" = "#628FB4",
  "250" = "#5685AD",
  "500" = "#497AA5",
  "898" = "#B85A54"
)

reference_size_plot <- ggplot(
  reference_size_resample_diagnostics,
  aes(
    x = reference_size_label,
    y = spearman_correlation_to_full_reference,
    fill = reference_size_label,
    colour = reference_size_label
  )
) +
  geom_hline(
    yintercept = 1,
    linetype = 2,
    linewidth = 0.5,
    colour = "#6E6E6E"
  ) +
  geom_violin(
    width = 0.82,
    trim = TRUE,
    scale = "width",
    adjust = 1,
    linewidth = 0.45,
    alpha = 0.92
  ) +
  geom_boxplot(
    width = 0.16,
    outlier.shape = NA,
    fill = "white",
    colour = "#333333",
    linewidth = 0.45
  ) +
  stat_summary(
    fun = median,
    geom = "point",
    shape = 21,
    size = 2.2,
    stroke = 0.45,
    fill = "white",
    colour = "#333333"
  ) +
  scale_fill_manual(values = reference_size_fill_colors, guide = "none") +
  scale_colour_manual(values = reference_size_line_colors, guide = "none") +
  scale_y_continuous(
    breaks = seq(0.86, 1.00, by = 0.02),
    labels = scales::label_number(accuracy = 0.01),
    expand = expansion(mult = c(0.01, 0.025))
  ) +
  coord_cartesian(ylim = c(0.86, 1.002), clip = "off") +
  labs(
    title = "Reference-size stability of LR perturbation profiles",
    subtitle = "1,000 bootstrap resamples per reference size",
    x = "GTEx normal reference size",
    y = expression("Spearman " * rho * " with full-reference profile")
  ) +
  theme_classic(base_size = 11, base_family = "sans", base_line_size = 0.5) +
  theme(
    axis.text = element_text(size = 10, colour = "black"),
    axis.title = element_text(size = 11, face = "bold", colour = "black"),
    axis.title.x = element_text(margin = margin(t = 7)),
    axis.title.y = element_text(margin = margin(r = 7)),
    axis.ticks = element_line(linewidth = 0.45, colour = "black"),
    plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
    plot.subtitle = element_text(size = 9.5, colour = "#404040", hjust = 0.5),
    plot.title.position = "plot",
    panel.grid = element_blank(),
    plot.margin = margin(t = 6, r = 8, b = 5, l = 6, unit = "pt")
  )

reference_size_plot

ggsave(
  "figures/06_Assess_SSN_Stability/reference_size_stability.png",
  reference_size_plot,
  width = 5,
  height = 4,
  dpi = 600,
  bg = "white"
)
ggsave(
  "figures/06_Assess_SSN_Stability/reference_size_stability.pdf",
  reference_size_plot,
  width = 5,
  height = 4,
  device = pdf,
  bg = "white"
)
