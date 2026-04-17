rm(list = ls())

library(tidyverse)

dis_ils <- read_tsv("snake/results/branchpoints/combined/anchored_percentile_overlap/dis_ils_percentile_summary.tsv")
ils_dis <- read_tsv("snake/results/branchpoints/combined/anchored_percentile_overlap/percentile_summary.tsv")

ggplot(ils_dis,
       aes(x = reference_percentile_cutoff, y = query_covered_percent)) +
  geom_line() +
  geom_ribbon(aes(ymax = query_covered_percent, ymin = 0), fill = "black", alpha = 0.5) +
  expand_limits(y = 0) +
  theme_classic() +
  theme(axis.ticks = element_line(color = "black"),
        axis.text = element_text(color = "black")) +
  labs(x = "intron count percentile in ILS",
       y = "percent of intron with DIS reads") ->
  percentile_ils_dis

ggplot(dis_ils,
       aes(x = reference_percentile_cutoff, y = query_covered_percent)) +
  geom_line() +
  geom_ribbon(aes(ymax = query_covered_percent, ymin = 0), fill = "black", alpha = 0.5) +
  expand_limits(y = 0) +
  theme_classic() +
  theme(axis.ticks = element_line(color = "black"),
        axis.text = element_text(color = "black")) +
  labs(x = "intron count percentile in DIS",
       y = "percent of intron with ILS reads") ->
  percentile_dis_ils

ggsave("plots/percentile_ils_dis.pdf", percentile_ils_dis, width = 3, height = 3)
ggsave("plots/percentile_dis_ils.pdf", percentile_dis_ils, width = 3, height = 3)

ils_dis %>%
  filter(reference_percentile_cutoff == 50) %>%
  select(query_covered_percent) %>%
  mutate(query_not_covered_percent = 100 - query_covered_percent) %>%
  pivot_longer(c(1:2)) %>%
  mutate(name = c("observed in DIS", "not observed in DIS") %>%
           fct_relevel("observed in DIS")) ->
  pie_chart_df

ggplot(pie_chart_df, aes(x = "", y = value, fill = name)) +
  geom_bar(stat = "identity", color = "white", linewidth = 1.5) +
  coord_polar("y", start = 0) +
  scale_fill_manual(values = c("#D33B76", "grey70")) +
  theme_void() +
  labs(fill = "",
       title = "top 50% of ILS introns") +
  geom_text(aes(label = paste0(round(value, 1), "%")), position = position_stack(vjust = 0.5)) +
  theme(legend.position = "bottom",
        plot.title = element_text(hjust = 0.5)) ->
  ils_dis_pie

ggsave("plots/ils_dis_pie.pdf", ils_dis_pie, width = 3, height = 3)
            