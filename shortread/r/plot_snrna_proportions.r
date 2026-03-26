rm(list = ls())

library(tidyverse)


setwd("/groups/plaschka/rupert.faraway/Vytautesome/spliceosome_disassembly/")

snrna_dir <- "shortread/snake/results/qc/snrna_counts/"

list.files(snrna_dir, full.names = T) %>%
  lapply(read_tsv) %>%
  set_names(
    list.files(snrna_dir) %>% str_remove(".snrna_counts.tsv")
  ) %>%
  bind_rows(.id = "sample") %>%
  mutate(ip = sample %>% word(1, sep = "_") %>% fct_relevel("ILS")) %>%
  mutate(snrna = str_remove(reference, "RN") %>%
           word(1, sep = "::") %>%
           word(1, sep = "-")) %>%
  mutate(snrna = case_when(str_detect(snrna, "U4ATAC") ~ "U4ATAC",
                           str_detect(snrna, "U6ATAC") ~ "U6ATAC",
                           str_detect(snrna, "U5") ~ "U5",
                           str_detect(snrna, "VU1") ~ "U1",
                           str_detect(snrna, "U6V") ~ "U6",
                           T ~ snrna)) ->
  snrna_counts


snrna_counts %>%
  group_by(ip, sample, snrna) %>%
  summarise(total_counts = sum(fragment_count)) %>%
  ungroup() %>%
  filter(!str_detect(snrna, "ENSG")) %>%
  group_by(sample) %>%
  mutate(prop_counts = total_counts / sum(total_counts)) %>%
  ungroup() ->
  snrna_summarised_counts
  
snrna_summarised_counts %>%
  filter(snrna %in% c("U1", "U2", "U4", "U5", "U6")) ->
  snrna_summarised_counts_canonical

snrna_summarised_counts_canonical %>%
  group_by(ip, snrna) %>%
  summarise(prop_counts = mean(prop_counts)) ->
  snrna_mean_counts_canonical

ggplot(snrna_summarised_counts_canonical,
       aes(x = snrna, y = prop_counts, color = factor(ip))) +
  geom_point(position = position_dodge2(width = 0.9)) +
  geom_bar(stat = "identity", data = snrna_mean_counts_canonical,
           position = position_dodge2(width = 0.9),
           fill = NA, show.legend = F) +
  labs(x = "", y = "proportion of total\nsnRNA counts", color = "") +
  theme_classic() +
  theme(axis.text = element_text(color = "black"),
        axis.ticks = element_line(color = "black")) +
  scale_color_manual(values = c("black", "#D33B76"))
  
snrna_summarised_counts_canonical %>%
  filter(ip == "ILS") %>%
  group_by(snrna) %>%
  summarise(mean_ils_prop = mean(prop_counts)) ->
  mean_ils_prop

snrna_summarised_counts_canonical %>%
  left_join(mean_ils_prop) %>%
  mutate(relative_to_ils = prop_counts / mean_ils_prop) ->
  relative_to_ils

ggplot(relative_to_ils,
       aes(x = snrna, y = relative_to_ils, color = factor(ip))) +
  geom_point(position = position_dodge2(width = 0.75)) +
  labs(x = "", y = "proportion of total\nsnRNA counts", color = "") +
  theme_classic() +
  scale_y_log10() +
  geom_hline(yintercept = 1, linetype = "dashed") +
  theme(axis.text = element_text(color = "black"),
        axis.ticks = element_line(color = "black")) +
  scale_color_manual(values = c("black", "#D33B76"))
