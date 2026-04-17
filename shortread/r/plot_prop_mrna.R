rm(list = ls())

library(tidyverse)
library(patchwork)

condition_meta <- read_tsv("snake/results/premrna/combined/three_prime_exonic_coverage.by_condition.tsv")

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

summary_df <- read_tsv("snake/results/premrna/combined/summary.by_sample.tsv") %>%
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

ggsave("plots/proportion_exonic_combined_plot.pdf", combined_plot, width = 7, height = 2.5)

