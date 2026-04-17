rm(list = ls())

library(tidyverse)
library(patchwork)

setwd("/Users/rupert.faraway/Documents/GitHub/spliceosome_disassembly/shortread/")


# Load --------------------------------------------------------------------

bp_meta <- read_tsv("snake/results/branchpoints/combined/anchored_enrichment_cutoff_sweep/none/metaprofile.by_condition.tsv") %>%
  mutate(condition = condition %>% fct_relevel("DIS"))

ss_meta <- read_tsv("snake/results/branchpoints/combined/anchored_enrichment_cutoff_sweep/none/three_prime_coverage.by_condition.tsv") %>%
  mutate(condition = condition %>% fct_relevel("DIS"))

# Plot --------------------------------------------------------------------

ggplot(bp_meta, aes(x = offset_nt, color = condition, fill = condition)) +
  geom_line(aes(y = mean_coverage_anchored_percent), alpha = 0.5, show.legend = F) +
  geom_ribbon(aes(ymax = mean_coverage_anchored_percent, ymin = 0), alpha = 0.5, linewidth = 0, show.legend = F) +
  geom_ribbon(aes(ymin = mean_coverage_anchored_percent - ci95_coverage_anchored_percent, 
                  ymax = mean_coverage_anchored_percent + ci95_coverage_anchored_percent),
              show.legend = F,
              alpha = 0.5, linewidth = 0) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  coord_cartesian(xlim = c(-40, 10)) +
  theme_classic() +
  theme(axis.ticks = element_line(color = "black"),
        axis.text = element_text(color = "black")) +
  labs(x = "distance from branchpoint (nt)",
       y = "average read coverage\n(percentage)",
       fill = "") +
  scale_color_manual(values = c("#D33B76", "black")) +
  scale_fill_manual(values = c("#D33B76", "black")) ->
  bp_meta_plot

ggplot(ss_meta, aes(x = offset_nt, color = condition, fill = condition)) +
  geom_line(aes(y = mean_coverage_spanning_percent), alpha = 0.5, show.legend = F) +
  geom_ribbon(aes(ymax = mean_coverage_spanning_percent, ymin = 0), alpha = 0.5, linewidth = 0,
              show.legend = F) +
  geom_ribbon(aes(ymin = mean_coverage_spanning_percent - ci95_coverage_spanning_percent, 
                  ymax = mean_coverage_spanning_percent + ci95_coverage_spanning_percent),
              alpha = 0.5, linewidth = 0) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  coord_cartesian(xlim = c(-10, 10)) +
  theme_classic() +
  theme(axis.ticks = element_line(color = "black"),
        axis.text = element_text(color = "black"),
        axis.line.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.y = element_blank(),
        axis.title.y = element_blank()) +
  guides(color = "none") +
  labs(x = "distance from 3'SS (nt)",
       y = "average read coverage\n(percentage)",
       fill = "")  +
  scale_color_manual(values = c("#D33B76", "black")) +
  scale_fill_manual(values = c("#D33B76", "black")) ->
  ss_meta_plot

(bp_meta_plot | ss_meta_plot) + plot_layout(widths = c(5, 2))


bp_feature_rects <- tibble(
  feature = "intron",
  xmin = -100,
  xmax = 100,
  ymin = 0.46,
  ymax = 0.54
)

bp_feature_points <- tibble(
  feature = "branchpoint",
  x = 0,
  y = 0.5
)

ss_feature_rects <- tibble(
  feature = c("intron", "exon"),
  xmin = c(-10, 0),
  xmax = c(0, 10),
  ymin = c(0.46, 0.35),
  ymax = c(0.54, 0.65)
)

bp_feature_plot <- ggplot() +
  geom_rect(
    data = bp_feature_rects,
    aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
    fill = "black",
    color = NA
  ) +
  geom_point(
    data = bp_feature_points,
    aes(x = x, y = y),
    shape = 21,
    size = 4,
    stroke = 0,
    fill = "black",
    color = "black"
  ) +
  coord_cartesian(xlim = c(-40, 10), ylim = c(0, 1)) +
  theme_classic() +
  theme(
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.line.x = element_blank(),
    
    axis.title.y = element_text(color = "transparent"),
    axis.text.y = element_text(color = "transparent"),
    axis.ticks.y = element_line(color = "transparent"),
    axis.line.y = element_line(color = "transparent"),
  ) +
  labs(y = "average read coverage\n(percentage)")

ss_feature_plot <- ggplot(ss_feature_rects) +
  geom_rect(
    aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
    fill = "black",
    color = NA
  ) +
  coord_cartesian(xlim = c(-10, 10), ylim = c(0, 1), expand = FALSE, clip = "off") +
  theme_classic() +
  theme(
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.line.x = element_blank(),
    
    axis.title.y = element_text(color = "transparent"),
    axis.text.y = element_text(color = "transparent"),
    axis.ticks.y = element_line(color = "transparent"),
    axis.line.y = element_line(color = "transparent"),
  ) +
  labs(y = "average read coverage\n(percentage)")


bp_panel <- (bp_meta_plot / bp_feature_plot) +
  plot_layout(heights = c(10, 1))

ss_panel <- (ss_meta_plot / ss_feature_plot) +
  plot_layout(heights = c(10, 1))

(bp_panel | ss_panel) +
  plot_layout(widths = c(5, 2)) ->
  combined_metaprofile_panel

ggsave("plots/combined_metaprofile_panel.pdf", combined_metaprofile_panel, width = 8, height = 2.5)

