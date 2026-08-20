rm(list = ls())

library(tidyverse)
library(patchwork)
library(scales)

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

input_dir <- required_arg("input-dir")
output_pdf <- required_arg("output-pdf")
output_png <- required_arg("output-png")
output_broad <- required_arg("output-broad")
output_subtypes <- required_arg("output-subtypes")
output_counts_pdf <- get_cli_arg(
  "output-counts-pdf",
  default = sub("\\.pdf$", "_mapped_pairs.pdf", output_pdf)
)
output_counts_png <- get_cli_arg(
  "output-counts-png",
  default = sub("\\.png$", "_mapped_pairs.png", output_png)
)
output_counts_tsv <- get_cli_arg(
  "output-counts-tsv",
  default = if (str_detect(output_broad, "\\.broad\\.tsv$")) {
    sub("\\.broad\\.tsv$", ".mapped_pairs.tsv", output_broad)
  } else {
    sub("\\.tsv$", "_mapped_pairs.tsv", output_broad)
  }
)

dir.create(dirname(output_pdf), recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(output_png), recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(output_broad), recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(output_subtypes), recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(output_counts_pdf), recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(output_counts_png), recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(output_counts_tsv), recursive = TRUE, showWarnings = FALSE)

sample_levels <- c("ILS_1", "ILS_2", "ILS_3", "DIS_1", "DIS_2", "DIS_3")
condition_levels <- c("ILS", "DIS")

raw_files <- list.files(input_dir, pattern = "\\.rna_content\\.tsv$", full.names = TRUE)
if (length(raw_files) == 0) {
  stop("No .rna_content.tsv files found in ", input_dir)
}

sample_df <- raw_files %>%
  set_names(basename(.) %>% str_remove("\\.rna_content\\.tsv$")) %>%
  map_dfr(
    \(path) read_tsv(path, show_col_types = FALSE),
    .id = "sample"
  ) %>%
  mutate(
    condition = str_extract(sample, "^[^_]+"),
    percent_total = fraction_of_total_input_pairs * 100
  )

main_category_map <- c(
  "rRNA" = "rRNA",
  "snRNA" = "snRNA",
  "protein_coding_intron" = "pre-mRNA",
  "protein_coding_exon" = "mRNA",
  "intergenic" = "intergenic"
)

other_subtype_map <- c(
  "mitochondrial" = "mitochondrial",
  "snoRNA" = "snoRNA",
  "scaRNA" = "scaRNA",
  "miRNA" = "miRNA",
  "tRNA" = "tRNA",
  "ribozyme" = "ribozyme",
  "vaultRNA" = "vaultRNA",
  "misc_RNA" = "misc_RNA",
  "lncRNA" = "lncRNA",
  "pseudogene" = "pseudogene",
  "other_annotated_genic" = "other annotated/ambiguous",
  "ambiguous" = "other annotated/ambiguous"
)

broad_levels <- c("rRNA", "snRNA", "pre-mRNA", "mRNA", "intergenic", "other RNA")
broad_palette <- c(
  "rRNA" = "#6C8EAD",
  "snRNA" = "#D35D47",
  "pre-mRNA" = "#4A7C59",
  "mRNA" = "#D8A23A",
  "intergenic" = "#9D7A5E",
  "other RNA" = "#8A8F98"
)

subtype_palette <- c(
  "ILS" = "#1F3C5B",
  "DIS" = "#B54736"
)

broad_df <- sample_df %>%
  filter(category != "genome_unmapped") %>%
  mutate(
    broad_category = case_when(
      category %in% names(main_category_map) ~ unname(main_category_map[category]),
      TRUE ~ "other RNA"
    )
  ) %>%
  group_by(sample, condition, broad_category) %>%
  summarise(percent_total = sum(percent_total), .groups = "drop") %>%
  group_by(sample, condition) %>%
  mutate(
    mapped_percent_of_input = sum(percent_total),
    percent_mapped = percent_total / mapped_percent_of_input * 100
  ) %>%
  ungroup() %>%
  mutate(
    sample = factor(sample, levels = sample_levels),
    condition = factor(condition, levels = condition_levels),
    broad_category = factor(broad_category, levels = broad_levels)
  ) %>%
  arrange(sample, broad_category)

subtype_df <- sample_df %>%
  filter(category %in% names(other_subtype_map)) %>%
  mutate(subtype = unname(other_subtype_map[category])) %>%
  group_by(sample, condition, subtype) %>%
  summarise(percent_total = sum(percent_total), .groups = "drop") %>%
  left_join(
    broad_df %>%
      distinct(sample, condition, mapped_percent_of_input),
    by = c("sample", "condition")
  ) %>%
  mutate(percent_mapped = percent_total / mapped_percent_of_input * 100) %>%
  group_by(subtype) %>%
  mutate(mean_percent = mean(percent_mapped)) %>%
  ungroup() %>%
  filter(mean_percent > 0) %>%
  mutate(
    sample = factor(sample, levels = sample_levels),
    condition = factor(condition, levels = condition_levels)
  )

subtype_levels <- subtype_df %>%
  distinct(subtype, mean_percent) %>%
  arrange(mean_percent) %>%
  pull(subtype)

subtype_df <- subtype_df %>%
  mutate(subtype = factor(subtype, levels = subtype_levels))

subtype_condition_df <- subtype_df %>%
  group_by(condition, subtype) %>%
  summarise(mean_percent = mean(percent_mapped), .groups = "drop") %>%
  mutate(
    condition = factor(condition, levels = condition_levels),
    subtype = factor(subtype, levels = subtype_levels)
  )

subtype_segments <- subtype_condition_df %>%
  pivot_wider(names_from = condition, values_from = mean_percent)

mapped_counts_df <- sample_df %>%
  filter(category != "genome_unmapped") %>%
  group_by(sample, condition) %>%
  summarise(mapped_pairs = sum(count), .groups = "drop") %>%
  mutate(
    sample = factor(sample, levels = sample_levels),
    condition = factor(condition, levels = condition_levels),
    mapped_pairs_millions = mapped_pairs / 1e6
  ) %>%
  arrange(sample)

write_tsv(broad_df, output_broad)
write_tsv(subtype_df, output_subtypes)
write_tsv(mapped_counts_df, output_counts_tsv)

label_df <- broad_df %>%
  group_by(sample, condition) %>%
  arrange(broad_category, .by_group = TRUE) %>%
  mutate(
    label_x = cumsum(percent_mapped) - percent_mapped / 2,
    label = label_number(accuracy = 0.1, suffix = "%")(percent_mapped)
  ) %>%
  ungroup() %>%
  filter(percent_mapped >= 4)

main_plot <- ggplot(
  broad_df,
  aes(x = percent_mapped, y = sample, fill = broad_category)
) +
  geom_col(
    width = 0.72,
    color = "white",
    linewidth = 0.3,
    position = position_stack(reverse = TRUE)
  ) +
  geom_text(
    data = label_df,
    aes(x = label_x, label = label),
    size = 3,
    color = "white",
    fontface = "bold"
  ) +
  facet_grid(condition ~ ., scales = "free_y", space = "free_y", switch = "y") +
  scale_fill_manual(values = broad_palette, breaks = broad_levels, drop = FALSE) +
  scale_x_continuous(
    limits = c(0, 100),
    breaks = seq(0, 100, 20),
    expand = expansion(mult = c(0, 0.01)),
    labels = label_number(suffix = "%")
  ) +
  labs(
    x = "Percent of mapped read pairs",
    y = NULL,
    fill = NULL
  ) +
  theme_minimal(base_size = 11) +
  theme(
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_line(color = "#E8E1D7", linewidth = 0.35),
    strip.placement = "outside",
    strip.text.y.left = element_text(angle = 0, face = "bold", color = "#31261E"),
    strip.background = element_rect(fill = "#F4EFE8", color = NA),
    axis.text = element_text(color = "#31261E"),
    axis.title.x = element_text(color = "#31261E"),
    legend.position = "bottom",
    legend.key.width = unit(0.9, "cm"),
    plot.background = element_rect(fill = "white", color = NA)
  )

subtype_plot <- ggplot(subtype_condition_df, aes(x = mean_percent, y = subtype, color = condition)) +
  geom_segment(
    data = subtype_segments,
    aes(x = ILS, xend = DIS, y = subtype, yend = subtype),
    inherit.aes = FALSE,
    linewidth = 1.15,
    color = "#D9D1C7"
  ) +
  geom_point(size = 3.2) +
  scale_color_manual(values = subtype_palette) +
  scale_x_continuous(
    expand = expansion(mult = c(0.02, 0.06)),
    labels = label_number(accuracy = 0.1, suffix = "%")
  ) +
  labs(
    x = "Mean percent of mapped read pairs",
    y = NULL,
    color = NULL
  ) +
  theme_minimal(base_size = 11) +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.major.x = element_line(color = "#E8E1D7", linewidth = 0.35),
    axis.text = element_text(color = "#31261E"),
    axis.title.x = element_text(color = "#31261E"),
    legend.position = "bottom",
    plot.background = element_rect(fill = "white", color = NA)
  )

combined_plot <- main_plot / subtype_plot +
  plot_layout(heights = c(3.3, 2.1), guides = "collect") +
  plot_annotation(
    title = "RNA content across ILS and DIS libraries",
    subtitle = paste(
      "Top: per-sample mapped composition.",
      "Bottom: minor RNA subtypes shown as condition means.",
      "Both panels exclude unmapped and residual unassigned reads from the denominator."
    ),
    theme = theme(
      plot.title = element_text(face = "bold", size = 15, color = "#241A13"),
      plot.subtitle = element_text(size = 10, color = "#5A4C40", margin = margin(b = 8))
    )
  ) &
  theme(legend.position = "bottom")

counts_plot <- ggplot(
  mapped_counts_df,
  aes(x = mapped_pairs, y = sample, fill = condition)
) +
  geom_col(width = 0.72, color = "white", linewidth = 0.3) +
  geom_text(
    aes(label = label_number(accuracy = 0.1, scale = 1e-6, suffix = "M")(mapped_pairs)),
    hjust = -0.08,
    size = 3.3,
    color = "#31261E",
    fontface = "bold"
  ) +
  facet_grid(condition ~ ., scales = "free_y", space = "free_y", switch = "y") +
  scale_fill_manual(values = subtype_palette, guide = "none") +
  scale_x_continuous(
    expand = expansion(mult = c(0, 0.14)),
    labels = label_number(accuracy = 0.1, scale = 1e-6, suffix = "M")
  ) +
  labs(
    title = "Mapped read pairs per sample",
    subtitle = "Counts are summed across all non-unmapped RNA categories.",
    x = "Mapped read pairs",
    y = NULL
  ) +
  theme_minimal(base_size = 11) +
  theme(
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_line(color = "#E8E1D7", linewidth = 0.35),
    strip.placement = "outside",
    strip.text.y.left = element_text(angle = 0, face = "bold", color = "#31261E"),
    strip.background = element_rect(fill = "#F4EFE8", color = NA),
    axis.text = element_text(color = "#31261E"),
    axis.title.x = element_text(color = "#31261E"),
    plot.title = element_text(face = "bold", size = 14, color = "#241A13"),
    plot.subtitle = element_text(size = 10, color = "#5A4C40", margin = margin(b = 8)),
    plot.background = element_rect(fill = "white", color = NA)
  )

ggsave(output_pdf, combined_plot, width = 9.4, height = 8.2)
ggsave(output_png, combined_plot, width = 9.4, height = 8.2, dpi = 300)
ggsave(output_counts_pdf, counts_plot, width = 8.2, height = 4.6)
ggsave(output_counts_png, counts_plot, width = 8.2, height = 4.6, dpi = 300)
