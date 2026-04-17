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

metaprofile_path <- required_arg("metaprofile")
output_path <- required_arg("output")
dir.create(dirname(output_path), recursive = TRUE, showWarnings = FALSE)

bp_meta <- read_tsv(metaprofile_path) %>%
  filter(condition == "DIS")

ggplot(bp_meta, aes(x = offset_nt, color = condition, fill = condition)) +
  geom_line(aes(y = mean_coverage_percent), alpha = 0.5, show.legend = F) +
  geom_ribbon(aes(ymax = mean_coverage_percent, ymin = 0), alpha = 0.5, linewidth = 0, show.legend = F) +
  geom_ribbon(aes(ymin = mean_coverage_percent - ci95_coverage_percent, 
                  ymax = mean_coverage_percent + ci95_coverage_percent),
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
  scale_color_manual(values = c("#D33B76")) +
  scale_fill_manual(values = c("#D33B76")) ->
  bp_meta_plot

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

bp_panel <- (bp_meta_plot / bp_feature_plot) +
  plot_layout(heights = c(10, 1))

ggsave(output_path, bp_panel, width = 6, height = 2.5)
