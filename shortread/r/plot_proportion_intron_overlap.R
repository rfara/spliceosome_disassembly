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

ils_dis_summary_path <- required_arg("ils-dis-summary")
dis_ils_summary_path <- required_arg("dis-ils-summary")
output_ils_dis_path <- required_arg("output-ils-dis")
output_dis_ils_path <- required_arg("output-dis-ils")
output_pie_path <- required_arg("output-pie")
dir.create(dirname(output_ils_dis_path), recursive = TRUE, showWarnings = FALSE)

dis_ils <- read_tsv(dis_ils_summary_path)
ils_dis <- read_tsv(ils_dis_summary_path)

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

ggsave(output_ils_dis_path, percentile_ils_dis, width = 3, height = 3)
ggsave(output_dis_ils_path, percentile_dis_ils, width = 3, height = 3)

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

ggsave(output_pie_path, ils_dis_pie, width = 3, height = 3)
            
