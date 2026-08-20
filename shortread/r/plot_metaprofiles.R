rm(list = ls())

library(tidyverse)
library(patchwork)

cli_args <- commandArgs(trailingOnly = TRUE)

get_cli_arg <- function(name, default = NULL) {
  inline_prefix <- paste0("--", name, "=")
  inline_value <- cli_args[str_starts(cli_args, fixed(inline_prefix))]
  if (length(inline_value) > 0) {
    return(str_remove(inline_value[[length(inline_value)]], fixed(inline_prefix)))
  }

  flag_index <- match(paste0("--", name), cli_args)
  if (!is.na(flag_index) && flag_index < length(cli_args)) {
    return(cli_args[[flag_index + 1]])
  }

  default
}

required_arg <- function(name) {
  value <- get_cli_arg(name)
  if (is.null(value) || !nzchar(value)) {
    stop("Missing required argument --", name)
  }
  value
}

branchpoint_metaprofile_path <- required_arg("branchpoint-metaprofile")
splice_site_metaprofile_path <- required_arg("splice-site-metaprofile")
output_path <- required_arg("output")
splice_site_metaprofile_all_path <- get_cli_arg("splice-site-metaprofile-all")
splice_site_output_all_path <- get_cli_arg("output-splice-site-all")
dir.create(dirname(output_path), recursive = TRUE, showWarnings = FALSE)
if (!is.null(splice_site_output_all_path) && nzchar(splice_site_output_all_path)) {
  dir.create(dirname(splice_site_output_all_path), recursive = TRUE, showWarnings = FALSE)
}

# Load --------------------------------------------------------------------

bp_meta <- read_tsv(branchpoint_metaprofile_path) %>%
  mutate(condition = condition %>% fct_relevel("DIS"))

ss_meta <- read_tsv(splice_site_metaprofile_path) %>%
  mutate(condition = condition %>% fct_relevel("DIS"))

ss_meta_all <- NULL
if (!is.null(splice_site_metaprofile_all_path) && nzchar(splice_site_metaprofile_all_path)) {
  ss_meta_all <- read_tsv(splice_site_metaprofile_all_path) %>%
    mutate(condition = condition %>% fct_relevel("DIS"))
}

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

build_ss_meta_plot <- function(ss_data, hide_y_axis = FALSE) {
  plot <- ggplot(ss_data, aes(x = offset_nt, color = condition, fill = condition)) +
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
          axis.text = element_text(color = "black")) +
    guides(color = "none") +
    labs(x = "distance from 3'SS (nt)",
         y = "average read coverage\n(percentage)",
         fill = "") +
    scale_color_manual(values = c("#D33B76", "black")) +
    scale_fill_manual(values = c("#D33B76", "black"))

  if (hide_y_axis) {
    plot <- plot +
      theme(
        axis.line.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.y = element_blank(),
        axis.title.y = element_blank()
      )
  }

  plot
}

ss_meta_plot <- build_ss_meta_plot(ss_meta, hide_y_axis = TRUE)

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

ggsave(output_path, combined_metaprofile_panel, width = 8, height = 2.5)

if (!is.null(ss_meta_all) && !is.null(splice_site_output_all_path) && nzchar(splice_site_output_all_path)) {
  ss_all_panel <- (build_ss_meta_plot(ss_meta_all, hide_y_axis = FALSE) / ss_feature_plot) +
    plot_layout(heights = c(10, 1))
  ggsave(splice_site_output_all_path, ss_all_panel, width = 3.2, height = 2.5)
}
