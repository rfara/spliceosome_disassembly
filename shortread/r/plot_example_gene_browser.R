rm(list = ls())

library(tidyverse)
library(patchwork)
library(data.table)

# Settings ----------------------------------------------------------------

example_gene_name <- "HNRNPA1"
example_region <- "chr12:54282250-54282820" # e.g. "chr16:70251298-70290506"; overrides the gene-centred region.
annotation_gene_names <- NULL # In coordinate mode, NULL annotates all overlapping shared genes.
region_padding_nt <- 100
max_reads_per_condition <- 300
random_seed <- 1

condition_levels <- c("DIS", "ILS")
condition_colours <- c(DIS = "#D33B76", ILS = "black")
unique_alignment_expression <- "([NH] == 1 || (![NH] && mapq > 0)) && cigar !~ \"N\""

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

as_optional_string <- function(value) {
  if (is.null(value) || length(value) == 0) {
    return(NULL)
  }
  if (length(value) > 1) {
    return(value)
  }
  if (is.na(value) || !nzchar(value)) {
    return(NULL)
  }
  if (tolower(value) %in% c("null", "none", "na")) {
    return(NULL)
  }
  value
}

example_gene_name <- as_optional_string(get_cli_arg("gene", example_gene_name))
example_region <- as_optional_string(get_cli_arg("region", example_region))
annotation_gene_names <- as_optional_string(get_cli_arg("annotation-genes", annotation_gene_names))
if (!is.null(annotation_gene_names)) {
  annotation_gene_names <- str_split(annotation_gene_names, ",", simplify = TRUE) %>%
    as.character() %>%
    str_trim() %>%
    discard(~ .x == "")
}
max_reads_per_condition <- as.integer(get_cli_arg("max-reads", max_reads_per_condition))
random_seed <- as.integer(get_cli_arg("seed", random_seed))


# Paths -------------------------------------------------------------------

find_shortread_dir <- function() {
  script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  script_path <- if (length(script_arg) > 0) {
    normalizePath(sub("^--file=", "", script_arg[[1]]), mustWork = TRUE)
  } else {
    NA_character_
  }

  candidates <- c(
    if (!is.na(script_path)) file.path(dirname(script_path), "..") else NA_character_,
    getwd(),
    file.path(getwd(), "shortread"),
    file.path(getwd(), "..")
  )

  candidates <- unique(normalizePath(candidates[!is.na(candidates)], mustWork = FALSE))
  for (candidate in candidates) {
    if (
      dir.exists(file.path(candidate, "results")) &&
        dir.exists(file.path(candidate, "r"))
    ) {
      return(candidate)
    }
  }

  stop("Could not find the shortread directory. Run from the repo root, shortread/, or with Rscript.")
}

find_samtools <- function() {
  candidates <- c(
    Sys.getenv("SAMTOOLS", unset = NA_character_),
    unname(Sys.which("samtools"))
  )
  candidates <- unique(candidates[!is.na(candidates) & nzchar(candidates)])
  candidates <- candidates[file.exists(candidates)]

  if (length(candidates) == 0) {
    stop("Could not find samtools. Set the SAMTOOLS environment variable to the samtools executable.")
  }

  candidates[[1]]
}

shortread_dir <- find_shortread_dir()

samtools <- get_cli_arg("samtools", find_samtools())

gtf_path <- get_cli_arg(
  "gtf",
  file.path(shortread_dir, "..", "annotation", "gencode.v44.primary_assembly.basic.annotation.gtf.gz")
)
shared_genes_path <- get_cli_arg(
  "shared-genes",
  file.path(shortread_dir, "results", "processed_data", "premrna_mrna", "combined", "shared_genes.tsv")
)
sample_summary_path <- get_cli_arg(
  "sample-summary",
  file.path(shortread_dir, "results", "processed_data", "premrna_mrna", "combined", "summary.by_sample.tsv")
)
branchpoint_reference_path <- get_cli_arg(
  "branchpoint-reference",
  file.path(
    shortread_dir,
    "results",
    "processed_data",
    "branchpoint_metaprofiles",
    "reference",
    "mane_select_top_branchpoints.tsv"
  )
)
dedup_dir <- get_cli_arg("dedup-dir", file.path(shortread_dir, "results", "dedup"))
output_stem_arg <- as_optional_string(get_cli_arg("output-stem", NULL))
stats_output_path <- as_optional_string(get_cli_arg("stats-output", NULL))
coverage_data_output_path <- as_optional_string(get_cli_arg("coverage-data-output", NULL))
reads_data_output_path <- as_optional_string(get_cli_arg("reads-data-output", NULL))
annotations_data_output_path <- as_optional_string(get_cli_arg("annotations-data-output", NULL))


# Helpers -----------------------------------------------------------------

chromosome_label <- function(chrom) {
  paste("chromosome", str_remove(chrom, "^chr"))
}

sanitize_plot_id <- function(value) {
  str_replace_all(value, "[^A-Za-z0-9]+", "_") %>%
    str_replace_all("^_+|_+$", "")
}

parse_genomic_region <- function(region_string) {
  region_string <- str_remove_all(region_string, ",")
  region_match <- str_match(region_string, "^([^:]+):(\\d+)-(\\d+)$")

  if (is.na(region_match[1, 1])) {
    stop("Region must look like chr16:70251298-70290506")
  }

  chrom <- region_match[1, 2]
  if (!str_starts(chrom, "chr")) {
    chrom <- paste0("chr", chrom)
  }

  start <- as.integer(region_match[1, 3])
  end <- as.integer(region_match[1, 4])

  if (is.na(start) || is.na(end) || start < 1 || end < start) {
    stop("Region has invalid start/end coordinates: ", region_string)
  }

  list(
    chrom = chrom,
    start = start,
    end = end,
    label = paste0(chrom, ":", start, "-", end)
  )
}

gtf_attribute <- function(attributes, key) {
  pattern <- paste0(key, ' "([^"]+)"')
  str_match(attributes, pattern)[, 2]
}

read_transcript_exons <- function(gtf, transcript_id) {
  command <- paste(
    "gzip -dc", shQuote(gtf),
    "| grep -F", shQuote(paste0('transcript_id "', transcript_id, '"')),
    "| awk 'BEGIN{FS=\"\\t\"; OFS=\"\\t\"} $3 == \"exon\" {print}'"
  )

  exons <- data.table::fread(
    cmd = command,
    sep = "\t",
    header = FALSE,
    quote = "",
    fill = TRUE,
    col.names = c(
      "chrom", "source", "feature", "start", "end", "score",
      "strand", "frame", "attributes"
    )
  ) %>%
    as_tibble() %>%
    transmute(
      chrom,
      start = as.integer(start),
      end = as.integer(end),
      strand,
      exon_number = as.integer(gtf_attribute(attributes, "exon_number"))
    ) %>%
    arrange(start, end)

  if (nrow(exons) == 0) {
    stop("No exon rows found in the GTF for transcript ", transcript_id)
  }

  exons
}

genes_overlapping_region <- function(genes, region) {
  genes %>%
    mutate(
      transcript_start = as.integer(transcript_start),
      transcript_end = as.integer(transcript_end)
    ) %>%
    filter(
      chrom == region$chrom,
      transcript_start <= region$end,
      transcript_end >= region$start
    )
}

select_gene_centred_region <- function(genes, gene_name) {
  if (is.null(gene_name)) {
    stop("Set example_gene_name when example_region is NULL.")
  }

  gene <- genes %>%
    filter(gene_name == .env$gene_name) %>%
    slice(1) %>%
    mutate(
      transcript_start = as.integer(transcript_start),
      transcript_end = as.integer(transcript_end)
    )

  if (nrow(gene) == 0) {
    stop("Gene ", gene_name, " was not found in ", shared_genes_path)
  }

  region <- list(
    chrom = gene$chrom[[1]],
    start = max(1L, gene$transcript_start[[1]] - region_padding_nt),
    end = gene$transcript_end[[1]] + region_padding_nt
  )
  region$label <- paste0(region$chrom, ":", region$start, "-", region$end)

  list(region = region, genes = gene)
}

select_coordinate_region <- function(genes, region_string, selected_gene_names = NULL) {
  region <- parse_genomic_region(region_string)
  annotation_genes <- genes_overlapping_region(genes, region)

  if (!is.null(selected_gene_names) && length(selected_gene_names) > 0) {
    annotation_genes <- annotation_genes %>%
      filter(gene_name %in% selected_gene_names)
  }

  list(region = region, genes = annotation_genes)
}

read_annotation_exons <- function(annotation_genes) {
  if (nrow(annotation_genes) == 0) {
    return(tibble(
      chrom = character(),
      start = integer(),
      end = integer(),
      strand = character(),
      exon_number = integer(),
      gene_id = character(),
      gene_name = character(),
      transcript_id = character(),
      annotation_y = integer()
    ))
  }

  map_dfr(seq_len(nrow(annotation_genes)), function(index) {
    gene <- annotation_genes[index, ]
    read_transcript_exons(gtf_path, gene$transcript_id[[1]]) %>%
      mutate(
        gene_id = gene$gene_id[[1]],
        gene_name = gene$gene_name[[1]],
        transcript_id = gene$transcript_id[[1]],
        annotation_y = gene$annotation_y[[1]],
        .before = chrom
      )
  })
}

select_branchpoints <- function(branchpoints, annotation_genes, region) {
  branchpoints <- branchpoints %>%
    mutate(
      branchpoint_position = as.integer(branchpoint_position),
      intron_start = as.integer(intron_start),
      intron_end = as.integer(intron_end)
    ) %>%
    filter(
      chrom == region$chrom,
      branchpoint_position >= region$start,
      branchpoint_position <= region$end
    )

  if (nrow(annotation_genes) == 0) {
    return(branchpoints %>% mutate(annotation_y = NA_integer_))
  }

  branchpoints %>%
    inner_join(
      annotation_genes %>%
        select(gene_id, transcript_id, annotation_y),
      by = c("gene_id", "transcript_id")
    )
}

cigar_tokens <- function(cigar) {
  if (is.na(cigar) || cigar == "*") {
    return(tibble(length = integer(), op = character()))
  }

  matches <- str_match_all(cigar, "([0-9]+)([MIDNSHP=X])")[[1]]
  if (nrow(matches) == 0) {
    return(tibble(length = integer(), op = character()))
  }

  tibble(
    length = as.integer(matches[, 2]),
    op = matches[, 3]
  )
}

cigar_reference_blocks <- function(position, cigar) {
  tokens <- cigar_tokens(cigar)
  if (nrow(tokens) == 0) {
    return(tibble(start = integer(), end = integer()))
  }

  reference_position <- as.integer(position)
  blocks <- vector("list", nrow(tokens))
  block_count <- 0

  for (index in seq_len(nrow(tokens))) {
    token_length <- tokens$length[[index]]
    token_op <- tokens$op[[index]]

    if (token_op %in% c("M", "=", "X")) {
      block_count <- block_count + 1
      blocks[[block_count]] <- tibble(
        start = reference_position,
        end = reference_position + token_length - 1L
      )
      reference_position <- reference_position + token_length
    } else if (token_op %in% c("D", "N")) {
      reference_position <- reference_position + token_length
    }
  }

  bind_rows(blocks[seq_len(block_count)])
}

read_samtools_depth <- function(bam, region) {
  command <- paste(
    shQuote(samtools),
    "view -b -f 2 -F 3852 -e",
    shQuote(unique_alignment_expression),
    shQuote(bam),
    shQuote(region$label),
    "|",
    shQuote(samtools),
    "depth -"
  )

  stderr_file <- tempfile()
  on.exit(unlink(stderr_file), add = TRUE)

  lines <- system(
    paste(command, "2>", shQuote(stderr_file)),
    intern = TRUE
  )

  status <- attr(lines, "status")
  if (!is.null(status) && status != 0) {
    stderr_text <- readLines(stderr_file, warn = FALSE)
    stop("samtools depth failed for ", bam, ":\n", paste(stderr_text, collapse = "\n"))
  }

  positions <- tibble(position = seq.int(region$start, region$end))
  if (length(lines) == 0) {
    return(positions %>% mutate(coverage_count = 0))
  }

  data.table::fread(
    text = paste(lines, collapse = "\n"),
    sep = "\t",
    header = FALSE,
    col.names = c("chrom", "position", "coverage_count"),
    data.table = FALSE
  ) %>%
    as_tibble() %>%
    filter(
      chrom == region$chrom,
      position >= region$start,
      position <= region$end
    ) %>%
    group_by(position) %>%
    summarise(coverage_count = sum(coverage_count), .groups = "drop") %>%
    right_join(positions, by = "position") %>%
    mutate(coverage_count = replace_na(coverage_count, 0L)) %>%
    arrange(position)
}

read_samtools_region <- function(bam, region) {
  command <- paste(
    shQuote(samtools),
    "view -f 2 -F 3852 -e",
    shQuote(unique_alignment_expression),
    shQuote(bam),
    shQuote(region$label)
  )

  stderr_file <- tempfile()
  on.exit(unlink(stderr_file), add = TRUE)

  lines <- system(
    paste(command, "2>", shQuote(stderr_file)),
    intern = TRUE
  )

  status <- attr(lines, "status")
  if (!is.null(status) && status != 0) {
    stderr_text <- if (file.exists(stderr_file)) {
      readLines(stderr_file, warn = FALSE)
    } else {
      character()
    }
    stop("samtools view failed for ", bam, ":\n", paste(stderr_text, collapse = "\n"))
  }

  if (length(lines) == 0) {
    return(tibble())
  }

  alignments <- data.table::fread(
    text = paste(lines, collapse = "\n"),
    sep = "\t",
    header = FALSE,
    quote = "",
    fill = TRUE,
    data.table = FALSE
  )

  colnames(alignments)[seq_len(min(11, ncol(alignments)))] <- c(
    "qname", "flag", "rname", "pos", "mapq", "cigar",
    "rnext", "pnext", "tlen", "seq", "qual"
  )[seq_len(min(11, ncol(alignments)))]

  alignments <- alignments %>%
    as_tibble() %>%
    mutate(
      flag = as.integer(flag),
      pos = as.integer(pos),
      mapq = as.integer(mapq),
      pnext = as.integer(pnext),
      tlen = as.integer(tlen)
    )

  alignments %>%
    filter(rnext == "=" | rnext == rname) %>%
    transmute(
      read_name = qname,
      flag,
      pos,
      cigar
    )
}

extract_sample_region <- function(sample_row, region) {
  sample <- sample_row$sample
  condition <- sample_row$condition
  bam <- file.path(dedup_dir, paste0(sample, ".dedup.bam"))
  bai <- paste0(bam, ".bai")

  if (!file.exists(bam)) {
    stop("Missing BAM for sample ", sample, ": ", bam)
  }
  if (!file.exists(bai)) {
    stop("Missing BAM index for sample ", sample, ": ", bai)
  }

  coverage <- read_samtools_depth(bam, region) %>%
    mutate(
      sample = .env$sample,
      condition = .env$condition,
      .before = position
    )

  alignments <- read_samtools_region(bam, region) %>%
    mutate(
      sample = .env$sample,
      condition = .env$condition
    )

  if (nrow(alignments) == 0) {
    return(list(
      coverage = coverage,
      reads = tibble()
    ))
  }

  list(
    coverage = coverage,
    reads = alignments
  )
}

summarise_coverage <- function(coverage, sample_summary) {
  coverage %>%
    left_join(
      sample_summary %>% select(sample, intronic_fragments),
      by = "sample"
    ) %>%
    mutate(
      condition = factor(condition, levels = condition_levels),
      coverage_per_million_intronic_reads = coverage_count * 1e6 / intronic_fragments
    ) %>%
    group_by(condition, position) %>%
    summarise(
      sample_count = n_distinct(sample),
      mean_coverage = mean(coverage_per_million_intronic_reads),
      .groups = "drop"
    )
}

read_span_from_cigar <- function(position, cigar, region) {
  blocks <- cigar_reference_blocks(position, cigar)
  if (nrow(blocks) == 0) {
    return(tibble(
      xmin = integer(),
      xmax = integer()
    ))
  }

  tibble(
    xmin = max(min(blocks$start), region$start),
    xmax = min(max(blocks$end), region$end)
  ) %>%
    filter(xmin <= xmax)
}

pack_nonoverlapping_reads <- function(reads) {
  if (nrow(reads) == 0) {
    return(reads %>% mutate(pack_lane = integer()))
  }

  reads <- reads %>%
    arrange(xmin, xmax, read_name)

  lane_ends <- numeric()
  pack_lanes <- integer(nrow(reads))

  for (read_index in seq_len(nrow(reads))) {
    available_lanes <- which(lane_ends < reads$xmin[[read_index]])
    if (length(available_lanes) == 0) {
      available_lane <- length(lane_ends) + 1L
      lane_ends[[available_lane]] <- -Inf
    } else {
      available_lane <- available_lanes[[1]]
    }

    pack_lanes[[read_index]] <- available_lane
    lane_ends[[available_lane]] <- reads$xmax[[read_index]]
  }

  reads %>%
    mutate(pack_lane = pack_lanes)
}

subsample_reads <- function(reads, region, sample_summary) {
  set.seed(random_seed)

  sample_read_caps <- sample_summary %>%
    distinct(sample, condition) %>%
    count(condition, name = "samples_in_condition") %>%
    mutate(max_reads_per_sample = ceiling(max_reads_per_condition / samples_in_condition))

  sampled_reads <- reads %>%
    mutate(condition = factor(condition, levels = condition_levels)) %>%
    left_join(sample_read_caps, by = "condition") %>%
    group_by(condition, sample) %>%
    group_modify(~ slice_sample(.x, n = min(nrow(.x), .x$max_reads_per_sample[[1]]))) %>%
    ungroup() %>%
    mutate(span = map2(pos, cigar, read_span_from_cigar, region = region)) %>%
    unnest(span) %>%
    select(-samples_in_condition, -max_reads_per_sample) %>%
    arrange(condition, sample, xmin, xmax, read_name)

  packed_reads <- sampled_reads %>%
    group_by(condition, sample) %>%
    group_modify(~ pack_nonoverlapping_reads(.x)) %>%
    ungroup()

  sample_gap <- 2
  lane_counts <- sample_summary %>%
    distinct(sample, condition) %>%
    arrange(condition, sample) %>%
    left_join(
      packed_reads %>%
        group_by(condition, sample) %>%
        summarise(row_count = max(pack_lane), read_count = n(), .groups = "drop"),
      by = c("condition", "sample")
    ) %>%
    mutate(
      row_count = replace_na(row_count, 1L),
      read_count = replace_na(read_count, 0L),
      y_offset = lag(cumsum(row_count + sample_gap), default = 0),
      y_mid = y_offset + pmax(row_count, 1) / 2,
      separator_y = y_offset - sample_gap / 2
    )

  packed_reads <- packed_reads %>%
    left_join(lane_counts, by = c("condition", "sample")) %>%
    mutate(
      y = y_offset + pack_lane,
      ymin = y - 0.505,
      ymax = y + 0.505
    )

  list(reads = packed_reads, lanes = lane_counts)
}

browser_theme <- function(base_size = 8) {
  theme_classic(base_size = base_size) +
    theme(
      axis.ticks = element_line(color = "black"),
      axis.text = element_text(color = "black"),
      legend.title = element_blank(),
      plot.margin = margin(2, 5, 2, 5)
    )
}


# Load --------------------------------------------------------------------

shared_genes <- read_tsv(shared_genes_path, show_col_types = FALSE)
branchpoint_reference <- read_tsv(branchpoint_reference_path, show_col_types = FALSE)
sample_summary <- read_tsv(sample_summary_path, show_col_types = FALSE) %>%
  mutate(condition = factor(condition, levels = condition_levels)) %>%
  filter(condition %in% condition_levels)

region_selection <- if (is.null(example_region)) {
  select_gene_centred_region(shared_genes, example_gene_name)
} else {
  select_coordinate_region(shared_genes, example_region, annotation_gene_names)
}

plot_region <- region_selection$region
annotation_genes <- region_selection$genes %>%
  arrange(transcript_start, transcript_end, gene_name) %>%
  mutate(
    annotation_y = row_number(),
    label_x = (
      pmax(transcript_start, plot_region$start) +
        pmin(transcript_end, plot_region$end)
    ) / 2
  )

annotation_exons <- read_annotation_exons(annotation_genes)
branchpoints <- select_branchpoints(branchpoint_reference, annotation_genes, plot_region)

message("Extracting alignments from ", plot_region$label, " with ", samtools)
message("Annotating ", nrow(annotation_genes), " shared gene(s) in the region")
message("Marking ", nrow(branchpoints), " branchpoint annotation(s) in the region")

sample_tracks <- sample_summary %>%
  arrange(condition, sample) %>%
  split(.$sample) %>%
  map(extract_sample_region, region = plot_region)

coverage_by_sample <- sample_tracks %>%
  map("coverage") %>%
  bind_rows()

read_rects <- sample_tracks %>%
  map("reads") %>%
  bind_rows()

coverage_by_condition <- summarise_coverage(coverage_by_sample, sample_summary)
subsampled <- subsample_reads(read_rects, plot_region, sample_summary)

message("Extracted ", nrow(read_rects), " candidate read alignments")
message(
  "Drawing ",
  n_distinct(paste(subsampled$reads$sample, subsampled$reads$read_name)),
  " read alignments across ",
  sum(subsampled$lanes$row_count),
  " packed rows"
)


# Plot --------------------------------------------------------------------

coverage_plot <- ggplot(
  coverage_by_condition,
  aes(x = position, color = condition, fill = condition)
) +
  geom_vline(
    data = branchpoints,
    aes(xintercept = branchpoint_position),
    linetype = "dashed",
    color = "black",
    alpha = 0.5,
    linewidth = 0.25
  ) +
  geom_line(aes(y = mean_coverage), linewidth = 0.45) +
  geom_ribbon(aes(ymax = mean_coverage, ymin = 0),
              alpha = 0.5, linewidth = 0, show.legend = F) +
  scale_color_manual(values = condition_colours, drop = FALSE) +
  scale_fill_manual(values = condition_colours, drop = FALSE) +
  scale_x_continuous(
    limits = c(plot_region$start, plot_region$end),
    labels = scales::label_comma(),
    expand = expansion(mult = 0)
  ) +
  labs(
    x = NULL,
    y = "mean read coverage\n(per million intronic reads)"
  ) +
  browser_theme() +
  theme(
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.line.x = element_blank(),
    legend.position = "top"
  )

reads_plot <- ggplot(
  subsampled$reads,
  aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = condition)
) +
  geom_vline(
    data = branchpoints,
    aes(xintercept = branchpoint_position),
    linetype = "dashed",
    color = "black",
    alpha = 0.5,
    linewidth = 0.25
  ) +
  # geom_segment(
  #   data = subsampled$lanes %>% filter(separator_y > 0),
  #   aes(
  #     x = plot_region$start,
  #     xend = plot_region$end,
  #     y = separator_y,
  #     yend = separator_y
  #   ),
  #   inherit.aes = FALSE,
  #   color = "grey88",
  #   linewidth = 0.25
  # ) +
  geom_rect(color = NA, alpha = 0.8, show.legend = FALSE) +
  scale_fill_manual(values = condition_colours, drop = FALSE) +
  scale_y_continuous(
    breaks = subsampled$lanes$y_mid,
    labels = subsampled$lanes$sample,
    expand = expansion(mult = c(0.02, 0.02))
  ) +
  scale_x_continuous(
    limits = c(plot_region$start, plot_region$end),
    labels = scales::label_comma(),
    expand = expansion(mult = 0)
  ) +
  labs(
    x = chromosome_label(plot_region$chrom),
    y = "subsampled reads"
  ) +
  browser_theme()

annotation_rects <- bind_rows(
  annotation_genes %>%
    transmute(
      feature = "intron",
      gene_id,
      xmin = pmax(transcript_start, plot_region$start),
      xmax = pmin(transcript_end, plot_region$end),
      ymin = annotation_y - 0.04,
      ymax = annotation_y + 0.04
    ),
  annotation_exons %>%
    transmute(
      feature = "exon",
      gene_id,
      xmin = pmax(start, plot_region$start),
      xmax = pmin(end, plot_region$end),
      ymin = annotation_y - 0.18,
      ymax = annotation_y + 0.18
    )
) %>%
  filter(xmin <= xmax)

annotation_source_data <- bind_rows(
  annotation_rects %>%
    left_join(
      annotation_genes %>%
        select(gene_id, gene_name, transcript_id, strand),
      by = "gene_id"
    ) %>%
    transmute(
      feature_type = feature,
      gene_name,
      transcript_id,
      strand,
      x = NA_real_,
      xmin,
      xmax,
      y = NA_real_,
      ymin,
      ymax,
      label = NA_character_
    ),
  branchpoints %>%
    transmute(
      feature_type = "branchpoint",
      gene_name,
      transcript_id,
      strand,
      x = as.numeric(branchpoint_position),
      xmin = NA_real_,
      xmax = NA_real_,
      y = as.numeric(annotation_y),
      ymin = NA_real_,
      ymax = NA_real_,
      label = NA_character_
    ),
  if (nrow(annotation_genes) == 0) {
    tibble(
      feature_type = "gene_label",
      gene_name = NA_character_,
      transcript_id = NA_character_,
      strand = NA_character_,
      x = mean(c(plot_region$start, plot_region$end)),
      xmin = NA_real_,
      xmax = NA_real_,
      y = 1,
      ymin = NA_real_,
      ymax = NA_real_,
      label = "no shared genes in region"
    )
  } else {
    annotation_genes %>%
      transmute(
        feature_type = "gene_label",
        gene_name,
        transcript_id,
        strand,
        x = label_x,
        xmin = NA_real_,
        xmax = NA_real_,
        y = annotation_y + 0.35,
        ymin = NA_real_,
        ymax = NA_real_,
        label = paste0(gene_name, " (", strand, ")")
      )
  }
) %>%
  arrange(feature_type, gene_name, transcript_id, x, xmin, xmax)

annotation_ylim <- c(0.5, max(1, nrow(annotation_genes)) + 0.65)

annotation_plot <- ggplot(annotation_rects) +
  geom_vline(
    data = branchpoints,
    aes(xintercept = branchpoint_position),
    linetype = "dashed",
    color = "black",
    alpha = 0.5,
    linewidth = 0.25
  ) +
  geom_rect(
    aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
    fill = "black",
    color = NA
  ) +
  geom_point(
    data = branchpoints %>% filter(!is.na(annotation_y)),
    aes(x = branchpoint_position, y = annotation_y),
    size = 2.5
  ) +
  geom_text(
    data = annotation_genes,
    aes(
      x = label_x,
      y = annotation_y + 0.35,
      label = paste0(gene_name, " (", strand, ")")
    ),
    size = 3
  ) +
  geom_text(
    data = if (nrow(annotation_genes) == 0) {
      tibble(
        x = mean(c(plot_region$start, plot_region$end)),
        y = 1,
        label = "no shared genes in region"
      )
    } else {
      tibble(x = numeric(), y = numeric(), label = character())
    },
    aes(x = x, y = y, label = label),
    size = 3
  ) +
  scale_x_continuous(
    limits = c(plot_region$start, plot_region$end),
    labels = scales::label_comma(),
    expand = expansion(mult = 0)
  ) +
  coord_cartesian(ylim = annotation_ylim, clip = "off") +
  labs(x = NULL, y = NULL) +
  browser_theme() +
  theme(
    axis.title = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    axis.line = element_blank()
  )

example_gene_plot <- coverage_plot / reads_plot / annotation_plot +
  plot_layout(heights = c(2.2, 2.1, max(0.7, 0.45 * max(1, nrow(annotation_genes)))))

example_gene_plot_coverage <- coverage_plot / annotation_plot +
  plot_layout(heights = c(3, 1))

example_gene_plot_reads <- reads_plot / annotation_plot +
  plot_layout(heights = c(3, 1))

plot_id <- if (is.null(example_region)) {
  sanitize_plot_id(example_gene_name)
} else {
  sanitize_plot_id(plot_region$label)
}

output_stem <- if (is.null(output_stem_arg)) {
  file.path(shortread_dir, "results", "plots", paste0("example_gene_browser_", plot_id))
} else {
  output_stem_arg
}
dir.create(dirname(output_stem), recursive = TRUE, showWarnings = FALSE)

if (!is.null(stats_output_path)) {
  dir.create(dirname(stats_output_path), recursive = TRUE, showWarnings = FALSE)
  tibble(
    region = plot_region$label,
    requested_gene = if (is.null(example_gene_name)) NA_character_ else example_gene_name,
    annotated_shared_genes = nrow(annotation_genes),
    branchpoint_annotations = nrow(branchpoints),
    candidate_read_alignments = nrow(read_rects),
    drawn_distinct_reads = n_distinct(paste(subsampled$reads$sample, subsampled$reads$read_name)),
    packed_read_rows = sum(subsampled$lanes$row_count),
    max_reads_per_condition = max_reads_per_condition,
    random_seed = random_seed
  ) %>%
    write_tsv(stats_output_path)
  message("Wrote ", stats_output_path)
}

if (!is.null(coverage_data_output_path)) {
  dir.create(dirname(coverage_data_output_path), recursive = TRUE, showWarnings = FALSE)
  coverage_by_condition %>%
    transmute(
      condition = as.character(condition),
      position,
      mean_coverage_per_million_intronic_reads = mean_coverage
    ) %>%
    arrange(condition, position) %>%
    write_tsv(coverage_data_output_path)
  message("Wrote ", coverage_data_output_path)
}

if (!is.null(reads_data_output_path)) {
  dir.create(dirname(reads_data_output_path), recursive = TRUE, showWarnings = FALSE)
  subsampled$reads %>%
    transmute(
      condition = as.character(condition),
      sample,
      read_name,
      xmin,
      xmax,
      ymin,
      ymax
    ) %>%
    arrange(condition, sample, ymin, xmin, xmax, read_name) %>%
    write_tsv(reads_data_output_path)
  message("Wrote ", reads_data_output_path)
}

if (!is.null(annotations_data_output_path)) {
  dir.create(dirname(annotations_data_output_path), recursive = TRUE, showWarnings = FALSE)
  annotation_source_data %>%
    write_tsv(annotations_data_output_path)
  message("Wrote ", annotations_data_output_path)
}

ggsave(paste0(output_stem, ".pdf"), example_gene_plot, width = 9, height = 5)
ggsave(paste0(output_stem, ".png"), example_gene_plot, width = 9, height = 5, dpi = 300)

ggsave(paste0(output_stem, ".coverage_only.pdf"), example_gene_plot_coverage, width = 8, height = 3)
ggsave(paste0(output_stem, ".reads_only.pdf"), example_gene_plot_reads, width = 8, height = 3)


message("Wrote ", paste0(output_stem, ".pdf"))
message("Wrote ", paste0(output_stem, ".png"))
