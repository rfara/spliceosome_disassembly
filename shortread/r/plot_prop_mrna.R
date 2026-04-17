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
summary_path <- required_arg("summary")
output_path <- required_arg("output")
dir.create(dirname(output_path), recursive = TRUE, showWarnings = FALSE)

condition_meta <- read_tsv(metaprofile_path)

ggplot(condition_meta, aes(x = offset_nt, color = condition, fill = condition)) +
  geom_line(aes(y = mean_coverage_percent_gene_reads)) +
  geom_ribbon(aes(ymax = mean_coverage_percent_gene_reads, ymin = 0),
              alpha = 0.5, linewidth = 0, show.legend = F) +
  geom_ribbon(aes(ymin = mean_coverage_percent_gene_reads - ci95_coverage_percent_gene_reads, 
                  ymax = mean_coverage_percent_gene_reads + ci95_coverage_percent_gene_reads),
              show.legend = F,
              alpha = 0.25, linewidth = 0) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  coord_cartesian(xlim = c(-150, 15)) +
  theme_classic() +
  theme(axis.ticks = element_line(color = "black"),
        axis.text = element_text(color = "black")) +
  guides(fill = "none") +
  labs(x = "distance from gene 3' end (nt)",
       y = "mean read coverage\n(percentage of genic reads)",
       fill = "",
       color = "") +
  scale_color_manual(values = c("black", "#D33B76"), breaks = c("ILS", "DIS")) +
  scale_fill_manual(values = c("black", "#D33B76"), breaks = c("ILS", "DIS")) ->
  metaprofile_3end

summary_df <- read_tsv(summary_path) %>%
  mutate(condition = fct_relevel(condition, "ILS"))

t.test(summary_df$mrna_percent_gene_reads ~ summary_df$condition) ->
  mrna_ttest

ggplot(summary_df,
       aes(x = condition, color = condition,
           y = mrna_percent_gene_reads)) +
  geom_point(position = position_dodge2(width = 0.7), shape = 21, size = 3, fill = "white") +
  geom_bar(aes(x = condition, y = mean_mrna_percent), 
           stat = "identity",
           fill = NA,
           show.legend = F,
           data = summary_df %>%
             group_by(condition) %>% 
             summarise(mean_mrna_percent = mean(mrna_percent_gene_reads))) +
  expand_limits(y = c(0, 100)) +
  annotate(geom = "text", 
           label = paste0("p = ", round(mrna_ttest$p.value, 5)),
           size = 3,
           x = 1.5,
           y = max(summary_df$mrna_percent_gene_reads) * 1.5) +
  theme_classic() +
  theme(axis.ticks = element_line(color = "black"),
        axis.text = element_text(color = "black")) +
  labs(x = "",
       y = "percentage of reads\nmapping to exons",
       color = "") +
  scale_color_manual(values = rev(c("#D33B76", "black"))) +
  scale_fill_manual(values = rev(c("#D33B76", "black"))) ->
  barplot_proportion

(metaprofile_3end | barplot_proportion) +
  plot_layout(widths = c(3, 1)) ->
  combined_plot

ggsave(output_path, combined_plot, width = 7, height = 2.5)
