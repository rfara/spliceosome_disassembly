#!/usr/bin/env python3

import argparse
import csv
from collections import defaultdict
from pathlib import Path


def parse_args():
    parser = argparse.ArgumentParser(
        description="Summarize read and feature counts used by each workflow plot."
    )
    parser.add_argument("--branchpoint-summary", required=True)
    parser.add_argument("--branchpoint-three-prime", required=True)
    parser.add_argument("--downstream-summary", required=True)
    parser.add_argument("--premrna-summary", required=True)
    parser.add_argument("--percentile-ils-dis", required=True)
    parser.add_argument("--percentile-dis-ils", required=True)
    parser.add_argument("--example-browser-stats", required=True)
    parser.add_argument("--combined-metaprofile-output", required=True)
    parser.add_argument("--proportion-reads-stop-at-bp-output", required=True)
    parser.add_argument("--downstream-exon-output", required=True)
    parser.add_argument("--premrna-output", required=True)
    parser.add_argument("--percentile-ils-dis-output", required=True)
    parser.add_argument("--percentile-dis-ils-output", required=True)
    parser.add_argument("--ils-dis-pie-output", required=True)
    parser.add_argument("--example-browser-outputs", required=True)
    parser.add_argument("--output", required=True)
    return parser.parse_args()


def read_tsv(path):
    with open(path) as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path, rows, fieldnames):
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def count_value(row, field):
    raw_value = row.get(field, 0)
    if raw_value in {"", None, "NA"}:
        return None
    return int(float(raw_value))


def format_per_sample_counts(sample_rows, field):
    parts = []
    for row in sample_rows:
        value = count_value(row, field)
        if value is None:
            continue
        parts.append(f"{row['sample']}={value}")
    return ";".join(parts)


def group_rows_by_condition(rows):
    grouped = defaultdict(list)
    for row in rows:
        grouped[row["condition"]].append(row)
    return grouped


def parse_example_output_paths(raw_value):
    outputs = [value for value in raw_value.split(";") if value]
    if len(outputs) != 4:
        raise ValueError(
            f"Expected four example browser outputs (pdf, png, coverage_only, reads_only); got {outputs}"
        )
    return {
        "combined_pdf": outputs[0],
        "combined_png": outputs[1],
        "coverage_only": outputs[2],
        "reads_only": outputs[3],
    }


def branchpoint_selection_label(analysis_min_reads):
    if analysis_min_reads is None or analysis_min_reads <= 0:
        return "all_anchor_positive_introns_union"
    return f"introns_included_in_analysis_with_min_{analysis_min_reads}_anchored_reads_in_every_sample"


def branchpoint_notes(analysis_min_reads, base_message):
    if analysis_min_reads is None or analysis_min_reads <= 0:
        return f"{base_message} All anchor-positive introns are included in analysis without shared-intron filtering."
    return f"{base_message} Introns included in analysis are defined by the minimum anchored-read threshold in every sample."


def build_condition_row(
    *,
    plot_id,
    output_files,
    panel,
    condition="",
    reference_condition="",
    query_condition="",
    feature_type="",
    feature_count="",
    secondary_feature_type="",
    secondary_feature_count="",
    read_type="",
    read_count="",
    per_sample_read_counts="",
    per_sample_feature_counts="",
    selection="",
    notes="",
):
    return {
        "plot_id": plot_id,
        "output_files": output_files,
        "panel": panel,
        "condition": condition,
        "reference_condition": reference_condition,
        "query_condition": query_condition,
        "feature_type": feature_type,
        "feature_count": feature_count,
        "secondary_feature_type": secondary_feature_type,
        "secondary_feature_count": secondary_feature_count,
        "read_type": read_type,
        "read_count": read_count,
        "per_sample_read_counts": per_sample_read_counts,
        "per_sample_feature_counts": per_sample_feature_counts,
        "selection": selection,
        "notes": notes,
    }


def percentile_row(rows, percentile):
    for row in rows:
        if int(row["reference_percentile_cutoff"]) == percentile:
            return row
    raise ValueError(f"Did not find percentile cutoff {percentile}")


def main():
    args = parse_args()

    branchpoint_summary_rows = read_tsv(args.branchpoint_summary)
    branchpoint_three_prime_rows = read_tsv(args.branchpoint_three_prime)
    downstream_summary_rows = read_tsv(args.downstream_summary)
    premrna_summary_rows = read_tsv(args.premrna_summary)
    percentile_ils_dis_rows = read_tsv(args.percentile_ils_dis)
    percentile_dis_ils_rows = read_tsv(args.percentile_dis_ils)
    example_browser_stats = read_tsv(args.example_browser_stats)

    if len(example_browser_stats) != 1:
        raise ValueError(f"Expected one row in {args.example_browser_stats}")

    example_outputs = parse_example_output_paths(args.example_browser_outputs)
    example_stats = example_browser_stats[0]

    fieldnames = [
        "plot_id",
        "output_files",
        "panel",
        "condition",
        "reference_condition",
        "query_condition",
        "feature_type",
        "feature_count",
        "secondary_feature_type",
        "secondary_feature_count",
        "read_type",
        "read_count",
        "per_sample_read_counts",
        "per_sample_feature_counts",
        "selection",
        "notes",
    ]

    rows = []

    branchpoint_by_condition = group_rows_by_condition(branchpoint_summary_rows)
    three_prime_at_zero_by_sample = {}
    for row in branchpoint_three_prime_rows:
        if int(row["offset_nt"]) != 0:
            continue
        three_prime_at_zero_by_sample[row["sample"]] = int(float(row["coverage_count"]))

    for condition, sample_rows in branchpoint_by_condition.items():
        analysis_introns = count_value(sample_rows[0], "introns_included_in_analysis")
        analysis_min_reads = count_value(
            sample_rows[0],
            "min_anchored_reads_all_samples_for_analysis_inclusion",
        )
        rows.append(
            build_condition_row(
                plot_id="combined_metaprofile_panel",
                output_files=args.combined_metaprofile_output,
                panel="branchpoint_coverage",
                condition=condition,
                feature_type="introns_included_in_analysis",
                feature_count=analysis_introns,
                read_type="anchored_fragments",
                read_count=sum(count_value(row, "anchored_fragments") for row in sample_rows),
                per_sample_read_counts=format_per_sample_counts(sample_rows, "anchored_fragments"),
                selection=branchpoint_selection_label(analysis_min_reads),
                notes=branchpoint_notes(
                    analysis_min_reads,
                    "Condition traces are replicate means on the intron set included in analysis.",
                ),
            )
        )
        rows.append(
            build_condition_row(
                plot_id="combined_metaprofile_panel",
                output_files=args.combined_metaprofile_output,
                panel="three_prime_splice_site_coverage",
                condition=condition,
                feature_type="introns_included_in_analysis",
                feature_count=analysis_introns,
                read_type="three_prime_spanning_fragments",
                read_count=sum(three_prime_at_zero_by_sample[row["sample"]] for row in sample_rows),
                per_sample_read_counts=";".join(
                    f"{row['sample']}={three_prime_at_zero_by_sample[row['sample']]}" for row in sample_rows
                ),
                selection=branchpoint_selection_label(analysis_min_reads),
                notes=branchpoint_notes(
                    analysis_min_reads,
                    "Read count equals aggregate 3'SS coverage_count at offset 0 on the same intron set included in analysis as the branchpoint panel.",
                ),
            )
        )
        rows.append(
            build_condition_row(
                plot_id="proportion_reads_stop_at_bp",
                output_files=args.proportion_reads_stop_at_bp_output,
                panel="main",
                condition=condition,
                feature_type="introns_included_in_analysis",
                feature_count=analysis_introns,
                read_type="anchored_fragments",
                read_count=sum(count_value(row, "anchored_fragments") for row in sample_rows),
                per_sample_read_counts=format_per_sample_counts(sample_rows, "anchored_fragments"),
                selection=branchpoint_selection_label(analysis_min_reads),
                notes=branchpoint_notes(
                    analysis_min_reads,
                    "Plot points are per-sample zero-or-plus-one branchpoint percentages on the intron set included in analysis.",
                ),
            )
        )

    downstream_by_condition = group_rows_by_condition(downstream_summary_rows)
    for condition, sample_rows in downstream_by_condition.items():
        if condition != "DIS":
            continue
        rows.append(
            build_condition_row(
                plot_id="dis_metaprofile_from_downstream_exon",
                output_files=args.downstream_exon_output,
                panel="main",
                condition=condition,
                feature_type="downstream_exon_spanning_introns",
                secondary_feature_type="reference_introns",
                secondary_feature_count=count_value(sample_rows[0], "reference_introns"),
                read_type="downstream_exon_spanning_fragments",
                read_count=sum(count_value(row, "downstream_exon_spanning_fragments") for row in sample_rows),
                per_sample_read_counts=format_per_sample_counts(sample_rows, "downstream_exon_spanning_fragments"),
                per_sample_feature_counts=format_per_sample_counts(sample_rows, "downstream_exon_spanning_introns"),
                selection=(
                    f"proper_pairs_with_fragment_end_at_least_"
                    f"{count_value(sample_rows[0], 'downstream_exon_min_offset_nt')}_nt_into_downstream_exon"
                ),
                notes="Condition trace is the DIS replicate mean; downstream-exon-spanning intron counts vary by sample.",
            )
        )

    premrna_by_condition = group_rows_by_condition(premrna_summary_rows)
    for condition, sample_rows in premrna_by_condition.items():
        shared_genes = count_value(sample_rows[0], "shared_genes")
        shared_min_reads = count_value(sample_rows[0], "shared_min_total_reads_all_samples")
        rows.append(
            build_condition_row(
                plot_id="proportion_exonic_combined_plot",
                output_files=args.premrna_output,
                panel="main",
                condition=condition,
                feature_type="shared_genes",
                feature_count=shared_genes,
                read_type="assigned_gene_fragments",
                read_count=sum(count_value(row, "assigned_gene_fragments") for row in sample_rows),
                per_sample_read_counts=format_per_sample_counts(sample_rows, "assigned_gene_fragments"),
                selection=f"shared_genes_with_min_{shared_min_reads}_total_gene_reads_in_every_sample",
                notes="Bar plot uses mRNA fraction among assigned gene fragments; metaprofile uses the same shared gene set.",
            )
        )

    ils_dis_full = percentile_row(percentile_ils_dis_rows, 1)
    dis_ils_full = percentile_row(percentile_dis_ils_rows, 1)
    ils_dis_p50 = percentile_row(percentile_ils_dis_rows, 50)

    rows.append(
        build_condition_row(
            plot_id="percentile_ils_dis",
            output_files=args.percentile_ils_dis_output,
            panel="main",
            reference_condition=ils_dis_full["reference_condition"],
            query_condition=ils_dis_full["query_condition"],
            feature_type="reference_anchor_positive_introns",
            feature_count=count_value(ils_dis_full, "total_reference_anchor_positive_introns"),
            read_type="reference_anchored_fragments",
            read_count=count_value(ils_dis_full, "reference_anchored_fragments_at_or_above_cutoff"),
            selection="full_ranked_reference_universe",
            notes="Plot varies by percentile cutoff; counts here describe the full ILS reference universe at the lowest cutoff.",
        )
    )
    rows.append(
        build_condition_row(
            plot_id="percentile_dis_ils",
            output_files=args.percentile_dis_ils_output,
            panel="main",
            reference_condition=dis_ils_full["reference_condition"],
            query_condition=dis_ils_full["query_condition"],
            feature_type="reference_anchor_positive_introns",
            feature_count=count_value(dis_ils_full, "total_reference_anchor_positive_introns"),
            read_type="reference_anchored_fragments",
            read_count=count_value(dis_ils_full, "reference_anchored_fragments_at_or_above_cutoff"),
            selection="full_ranked_reference_universe",
            notes="Plot varies by percentile cutoff; counts here describe the full DIS reference universe at the lowest cutoff.",
        )
    )
    rows.append(
        build_condition_row(
            plot_id="ils_dis_pie",
            output_files=args.ils_dis_pie_output,
            panel="main",
            reference_condition=ils_dis_p50["reference_condition"],
            query_condition=ils_dis_p50["query_condition"],
            feature_type="reference_introns_at_or_above_cutoff",
            feature_count=count_value(ils_dis_p50, "reference_introns_at_or_above_cutoff"),
            read_type="reference_anchored_fragments_at_or_above_cutoff",
            read_count=count_value(ils_dis_p50, "reference_anchored_fragments_at_or_above_cutoff"),
            selection="reference_percentile_cutoff_50",
            notes="Pie chart summarizes the top 50% of ILS introns by anchored fragment count.",
        )
    )

    region_label = example_stats["region"]
    requested_gene = example_stats["requested_gene"]
    shared_gene_count = count_value(example_stats, "annotated_shared_genes")
    branchpoint_count = count_value(example_stats, "branchpoint_annotations")
    candidate_read_alignments = count_value(example_stats, "candidate_read_alignments")
    drawn_distinct_reads = count_value(example_stats, "drawn_distinct_reads")
    packed_read_rows = count_value(example_stats, "packed_read_rows")
    max_reads_per_condition = count_value(example_stats, "max_reads_per_condition")
    random_seed = count_value(example_stats, "random_seed")
    example_note = (
        f"Region {region_label}; requested gene {requested_gene}; "
        f"max_reads_per_condition={max_reads_per_condition}; random_seed={random_seed}."
    )

    rows.append(
        build_condition_row(
            plot_id="example_gene_browser",
            output_files=";".join([example_outputs["combined_pdf"], example_outputs["combined_png"]]),
            panel="coverage",
            feature_type="annotated_shared_genes",
            feature_count=shared_gene_count,
            secondary_feature_type="branchpoint_annotations",
            secondary_feature_count=branchpoint_count,
            read_type="candidate_read_alignments",
            read_count=candidate_read_alignments,
            selection="region_overlap",
            notes=example_note,
        )
    )
    rows.append(
        build_condition_row(
            plot_id="example_gene_browser",
            output_files=";".join([example_outputs["combined_pdf"], example_outputs["combined_png"]]),
            panel="reads",
            feature_type="annotated_shared_genes",
            feature_count=shared_gene_count,
            secondary_feature_type="branchpoint_annotations",
            secondary_feature_count=branchpoint_count,
            read_type="drawn_distinct_reads",
            read_count=drawn_distinct_reads,
            selection="region_overlap_after_plot_subsampling",
            notes=f"{example_note} Packed read rows={packed_read_rows}.",
        )
    )
    rows.append(
        build_condition_row(
            plot_id="example_gene_browser_coverage_only",
            output_files=example_outputs["coverage_only"],
            panel="coverage",
            feature_type="annotated_shared_genes",
            feature_count=shared_gene_count,
            secondary_feature_type="branchpoint_annotations",
            secondary_feature_count=branchpoint_count,
            read_type="candidate_read_alignments",
            read_count=candidate_read_alignments,
            selection="region_overlap",
            notes=example_note,
        )
    )
    rows.append(
        build_condition_row(
            plot_id="example_gene_browser_reads_only",
            output_files=example_outputs["reads_only"],
            panel="reads",
            feature_type="annotated_shared_genes",
            feature_count=shared_gene_count,
            secondary_feature_type="branchpoint_annotations",
            secondary_feature_count=branchpoint_count,
            read_type="drawn_distinct_reads",
            read_count=drawn_distinct_reads,
            selection="region_overlap_after_plot_subsampling",
            notes=f"{example_note} Packed read rows={packed_read_rows}.",
        )
    )

    write_tsv(args.output, rows, fieldnames)

    print(f"Wrote {len(rows)} plot analysis rows to {args.output}")


if __name__ == "__main__":
    main()
