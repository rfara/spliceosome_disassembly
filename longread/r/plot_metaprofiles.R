#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(readr)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 5) {
  stop("Usage: plot_metaprofiles.R COVERAGE_TSV DELETION_TSV AUC_TSV COVERAGE_PDF DELETION_PDF")
}

coverage_path <- args[[1]]
deletion_path <- args[[2]]
auc_path <- args[[3]]
coverage_pdf <- args[[4]]
deletion_pdf <- args[[5]]

dir.create(dirname(coverage_pdf), recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(deletion_pdf), recursive = TRUE, showWarnings = FALSE)

protein_levels <- c("tfip11", "dhx35")
protein_colors <- c(tfip11 = "violetred3", dhx35 = "darkorchid4")

coverage <- read_tsv(coverage_path, show_col_types = FALSE) %>%
  mutate(protein = factor(protein, levels = protein_levels))

coverage_plot <- ggplot(
  coverage,
  aes(
    x = strand_corrected_position,
    y = normalised_coverage,
    color = protein
  )
) +
  geom_line() +
  geom_vline(xintercept = 0, linetype = "dashed") +
  theme_classic() +
  theme(
    axis.ticks = element_line(color = "black"),
    axis.text = element_text(color = "black")
  ) +
  labs(
    x = "position relative to 3' splice site (nt)",
    y = "normalised read coverage",
    color = ""
  ) +
  scale_color_manual(values = protein_colors) +
  coord_cartesian(ylim = c(0, 1), xlim = c(-70, 30))

ggsave(coverage_pdf, coverage_plot, width = 6, height = 3)

deletion <- read_tsv(deletion_path, show_col_types = FALSE) %>%
  mutate(protein = factor(protein, levels = protein_levels))
auc <- read_tsv(auc_path, show_col_types = FALSE)

ribbon <- deletion %>%
  inner_join(auc %>% select(protein, baseline), by = "protein") %>%
  filter(strand_corrected_position >= -4, strand_corrected_position <= 2)

deletion_plot <- ggplot(
  deletion,
  aes(
    x = strand_corrected_position,
    y = coverage_normalised,
    color = protein
  )
) +
  geom_line() +
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_vline(xintercept = c(-4, 2), alpha = 0.5, linetype = "dotted") +
  geom_ribbon(
    data = ribbon,
    aes(ymin = baseline, ymax = coverage_normalised, fill = protein),
    color = NA,
    alpha = 0.3
  ) +
  annotate(
    "curve", x = -20, y = 0.07, xend = -1, yend = 0.06,
    linewidth = 0.5, curvature = -0.25,
    arrow = arrow(length = grid::unit(0.25, "cm"))
  ) +
  annotate(
    "curve", x = -20, y = 0.055, xend = -1, yend = 0.041,
    linewidth = 0.5, curvature = -0.25,
    arrow = arrow(length = grid::unit(0.25, "cm"))
  ) +
  annotate(
    "text", x = -24, y = 0.07,
    label = round(auc$area_under_peak[auc$protein == "tfip11"], 3), size = 3
  ) +
  annotate(
    "text", x = -24, y = 0.055,
    label = round(auc$area_under_peak[auc$protein == "dhx35"], 3), size = 3
  ) +
  theme_classic() +
  theme(
    axis.ticks = element_line(color = "black"),
    axis.text = element_text(color = "black")
  ) +
  labs(
    x = "position relative to branchpoint (nt)",
    y = "deletion rate per read",
    color = "",
    fill = ""
  ) +
  scale_color_manual(values = protein_colors) +
  scale_fill_manual(values = protein_colors, guide = "none") +
  coord_cartesian(ylim = c(0.02, 0.08), xlim = c(-60, 10))

ggsave(deletion_pdf, deletion_plot, width = 6, height = 3)
