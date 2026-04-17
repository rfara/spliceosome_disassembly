rm(list = ls())

library(tidyverse)

bp_reads <- read_tsv("snake/results/branchpoints/combined/no_filtering/summary.by_sample.tsv") %>%
  mutate(condition = fct_relevel(condition, "ILS"))


t.test(bp_reads$zero_or_plus_one_branchpoint_percent_anchored ~ bp_reads$condition) ->
  bp_ttest

ggplot(bp_reads, aes(x = condition, color = condition, y = zero_or_plus_one_branchpoint_percent_anchored)) +
  geom_point(position = position_dodge2(width = 0.7), shape = 21, size = 3, fill = "white") +
  geom_bar(aes(x = condition, y = mean_bp_percent), 
           stat = "identity",
           fill = NA,
           show.legend = F,
           data = bp_reads %>%
             group_by(condition) %>% 
             summarise(mean_bp_percent = mean(zero_or_plus_one_branchpoint_percent_anchored))) +
  expand_limits(y = c(0, 65)) +
  annotate(geom = "text", 
           label = paste0("p = ", round(bp_ttest$p.value, 5)),
           size = 3,
           x = 1.5,
           y = max(bp_reads$zero_or_plus_one_branchpoint_percent_anchored) * 1.1) +
  theme_classic() +
  theme(axis.ticks = element_line(color = "black"),
        axis.text = element_text(color = "black")) +
  labs(x = "",
       y = "percentage of reads\nstopping at BP",
       color = "") +
  scale_color_manual(values = rev(c("#D33B76", "black"))) +
  scale_fill_manual(values = rev(c("#D33B76", "black"))) ->
  barplot_proportion

ggsave("plots/proportion_reads_stop_at_bp.pdf", barplot_proportion, width = 3, height = 3)


