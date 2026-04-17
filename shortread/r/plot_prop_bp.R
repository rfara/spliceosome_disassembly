rm(list = ls())

library(tidyverse)

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

summary_path <- required_arg("summary")
output_path <- required_arg("output")
dir.create(dirname(output_path), recursive = TRUE, showWarnings = FALSE)

bp_reads <- read_tsv(summary_path) %>%
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

ggsave(output_path, barplot_proportion, width = 3, height = 3)

