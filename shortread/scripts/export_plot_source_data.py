#!/usr/bin/env python3

import argparse
import csv
from pathlib import Path


def parse_args():
    parser = argparse.ArgumentParser(
        description="Write compact supplementary tables with the data underlying each plot."
    )
    parser.add_argument("--branchpoint-metaprofile", required=True)
    parser.add_argument("--branchpoint-three-prime", required=True)
    parser.add_argument("--branchpoint-summary", required=True)
    parser.add_argument("--downstream-metaprofile", required=True)
    parser.add_argument("--premrna-metaprofile", required=True)
    parser.add_argument("--premrna-summary", required=True)
    parser.add_argument("--percentile-ils-dis", required=True)
    parser.add_argument("--percentile-dis-ils", required=True)
    parser.add_argument("--example-browser-coverage", required=True)
    parser.add_argument("--example-browser-reads", required=True)
    parser.add_argument("--example-browser-annotations", required=True)
    parser.add_argument("--combined-metaprofile-output", required=True)
    parser.add_argument("--proportion-reads-stop-at-bp-output", required=True)
    parser.add_argument("--downstream-exon-output", required=True)
    parser.add_argument("--premrna-output", required=True)
    parser.add_argument("--percentile-ils-dis-output", required=True)
    parser.add_argument("--percentile-dis-ils-output", required=True)
    parser.add_argument("--ils-dis-pie-output", required=True)
    parser.add_argument("--example-browser-outputs", required=True)
    parser.add_argument("--output-branchpoint-coverage", required=True)
    parser.add_argument("--output-three-prime-coverage", required=True)
    parser.add_argument("--output-branchpoint-stop", required=True)
    parser.add_argument("--output-downstream-exon", required=True)
    parser.add_argument("--output-premrna-metaprofile", required=True)
    parser.add_argument("--output-premrna-exonic-fraction", required=True)
    parser.add_argument("--output-percentile-ils-dis", required=True)
    parser.add_argument("--output-percentile-dis-ils", required=True)
    parser.add_argument("--output-ils-dis-pie", required=True)
    parser.add_argument("--output-manifest", required=True)
    return parser.parse_args()


def read_tsv(path):
    with open(path) as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path, rows, fieldnames):
    output_path = Path(path)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with output_path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def select_columns(rows, columns):
    return [{column: row[column] for column in columns} for row in rows]


def parse_example_output_paths(raw_value):
    outputs = [value for value in raw_value.split(";") if value]
    if len(outputs) != 4:
        raise ValueError(
            f"Expected four example browser outputs (pdf, png, coverage_only, reads_only); got {outputs}"
        )
    return outputs


def percentile_row(rows, percentile):
    for row in rows:
        if int(row["reference_percentile_cutoff"]) == percentile:
            return row
    raise ValueError(f"Missing percentile cutoff {percentile}")


def main():
    args = parse_args()

    branchpoint_metaprofile_rows = read_tsv(args.branchpoint_metaprofile)
    branchpoint_three_prime_rows = read_tsv(args.branchpoint_three_prime)
    branchpoint_summary_rows = read_tsv(args.branchpoint_summary)
    downstream_metaprofile_rows = read_tsv(args.downstream_metaprofile)
    premrna_metaprofile_rows = read_tsv(args.premrna_metaprofile)
    premrna_summary_rows = read_tsv(args.premrna_summary)
    percentile_ils_dis_rows = read_tsv(args.percentile_ils_dis)
    percentile_dis_ils_rows = read_tsv(args.percentile_dis_ils)
    example_browser_outputs = parse_example_output_paths(args.example_browser_outputs)

    branchpoint_coverage_rows = select_columns(
        branchpoint_metaprofile_rows,
        [
            "condition",
            "offset_nt",
            "mean_coverage_anchored_percent",
            "ci95_coverage_anchored_percent",
        ],
    )
    write_tsv(
        args.output_branchpoint_coverage,
        branchpoint_coverage_rows,
        [
            "condition",
            "offset_nt",
            "mean_coverage_anchored_percent",
            "ci95_coverage_anchored_percent",
        ],
    )

    three_prime_coverage_rows = select_columns(
        branchpoint_three_prime_rows,
        [
            "condition",
            "offset_nt",
            "mean_coverage_spanning_percent",
            "ci95_coverage_spanning_percent",
        ],
    )
    write_tsv(
        args.output_three_prime_coverage,
        three_prime_coverage_rows,
        [
            "condition",
            "offset_nt",
            "mean_coverage_spanning_percent",
            "ci95_coverage_spanning_percent",
        ],
    )

    branchpoint_stop_rows = select_columns(
        branchpoint_summary_rows,
        [
            "sample",
            "condition",
            "zero_or_plus_one_branchpoint_percent_anchored",
        ],
    )
    write_tsv(
        args.output_branchpoint_stop,
        branchpoint_stop_rows,
        [
            "sample",
            "condition",
            "zero_or_plus_one_branchpoint_percent_anchored",
        ],
    )

    downstream_exon_rows = select_columns(
        downstream_metaprofile_rows,
        [
            "condition",
            "offset_nt",
            "mean_coverage_percent",
            "ci95_coverage_percent",
        ],
    )
    write_tsv(
        args.output_downstream_exon,
        downstream_exon_rows,
        [
            "condition",
            "offset_nt",
            "mean_coverage_percent",
            "ci95_coverage_percent",
        ],
    )

    premrna_metaprofile_source_rows = select_columns(
        premrna_metaprofile_rows,
        [
            "condition",
            "offset_nt",
            "mean_coverage_percent_gene_reads",
            "ci95_coverage_percent_gene_reads",
        ],
    )
    write_tsv(
        args.output_premrna_metaprofile,
        premrna_metaprofile_source_rows,
        [
            "condition",
            "offset_nt",
            "mean_coverage_percent_gene_reads",
            "ci95_coverage_percent_gene_reads",
        ],
    )

    premrna_fraction_rows = select_columns(
        premrna_summary_rows,
        [
            "sample",
            "condition",
            "mrna_percent_gene_reads",
        ],
    )
    write_tsv(
        args.output_premrna_exonic_fraction,
        premrna_fraction_rows,
        [
            "sample",
            "condition",
            "mrna_percent_gene_reads",
        ],
    )

    percentile_ils_dis_source_rows = select_columns(
        percentile_ils_dis_rows,
        [
            "reference_condition",
            "query_condition",
            "reference_percentile_cutoff",
            "query_covered_percent",
        ],
    )
    write_tsv(
        args.output_percentile_ils_dis,
        percentile_ils_dis_source_rows,
        [
            "reference_condition",
            "query_condition",
            "reference_percentile_cutoff",
            "query_covered_percent",
        ],
    )

    percentile_dis_ils_source_rows = select_columns(
        percentile_dis_ils_rows,
        [
            "reference_condition",
            "query_condition",
            "reference_percentile_cutoff",
            "query_covered_percent",
        ],
    )
    write_tsv(
        args.output_percentile_dis_ils,
        percentile_dis_ils_source_rows,
        [
            "reference_condition",
            "query_condition",
            "reference_percentile_cutoff",
            "query_covered_percent",
        ],
    )

    ils_dis_at_50 = percentile_row(percentile_ils_dis_rows, 50)
    pie_rows = [
        {
            "category": "observed in DIS",
            "percent": ils_dis_at_50["query_covered_percent"],
        },
        {
            "category": "not observed in DIS",
            "percent": str(100.0 - float(ils_dis_at_50["query_covered_percent"])),
        },
    ]
    write_tsv(args.output_ils_dis_pie, pie_rows, ["category", "percent"])

    manifest_rows = [
        {
            "plot_id": "combined_metaprofile_panel",
            "panel_or_view": "branchpoint_coverage",
            "plot_files": args.combined_metaprofile_output,
            "table_file": Path(args.output_branchpoint_coverage).name,
            "description": "Condition mean branchpoint-centred coverage percentages with 95% confidence intervals.",
        },
        {
            "plot_id": "combined_metaprofile_panel",
            "panel_or_view": "three_prime_splice_site_coverage",
            "plot_files": args.combined_metaprofile_output,
            "table_file": Path(args.output_three_prime_coverage).name,
            "description": "Condition mean 3' splice site coverage percentages with 95% confidence intervals.",
        },
        {
            "plot_id": "proportion_reads_stop_at_bp",
            "panel_or_view": "main",
            "plot_files": args.proportion_reads_stop_at_bp_output,
            "table_file": Path(args.output_branchpoint_stop).name,
            "description": "Per-sample percentages of anchored reads stopping at the branchpoint or one nucleotide downstream.",
        },
        {
            "plot_id": "dis_metaprofile_from_downstream_exon",
            "panel_or_view": "main",
            "plot_files": args.downstream_exon_output,
            "table_file": Path(args.output_downstream_exon).name,
            "description": "Condition mean DIS downstream-exon metaprofile coverage percentages with 95% confidence intervals.",
        },
        {
            "plot_id": "proportion_exonic_combined_plot",
            "panel_or_view": "metaprofile",
            "plot_files": args.premrna_output,
            "table_file": Path(args.output_premrna_metaprofile).name,
            "description": "Condition mean coverage near gene 3' ends, expressed as percentages of gene-assigned reads, with 95% confidence intervals.",
        },
        {
            "plot_id": "proportion_exonic_combined_plot",
            "panel_or_view": "barplot",
            "plot_files": args.premrna_output,
            "table_file": Path(args.output_premrna_exonic_fraction).name,
            "description": "Per-sample percentages of assigned gene fragments mapping to exons.",
        },
        {
            "plot_id": "percentile_ils_dis",
            "panel_or_view": "main",
            "plot_files": args.percentile_ils_dis_output,
            "table_file": Path(args.output_percentile_ils_dis).name,
            "description": "Percent of ILS-ranked introns that are observed in DIS across percentile cutoffs.",
        },
        {
            "plot_id": "percentile_dis_ils",
            "panel_or_view": "main",
            "plot_files": args.percentile_dis_ils_output,
            "table_file": Path(args.output_percentile_dis_ils).name,
            "description": "Percent of DIS-ranked introns that are observed in ILS across percentile cutoffs.",
        },
        {
            "plot_id": "ils_dis_pie",
            "panel_or_view": "main",
            "plot_files": args.ils_dis_pie_output,
            "table_file": Path(args.output_ils_dis_pie).name,
            "description": "Observed and not-observed proportions in DIS for the top 50% of ILS introns.",
        },
        {
            "plot_id": "example_gene_browser",
            "panel_or_view": "coverage",
            "plot_files": ";".join(example_browser_outputs),
            "table_file": Path(args.example_browser_coverage).name,
            "description": "Condition mean coverage track for the plotted genomic region.",
        },
        {
            "plot_id": "example_gene_browser",
            "panel_or_view": "reads",
            "plot_files": ";".join(example_browser_outputs),
            "table_file": Path(args.example_browser_reads).name,
            "description": "Exact subsampled read rectangles drawn in the gene browser read track.",
        },
        {
            "plot_id": "example_gene_browser",
            "panel_or_view": "annotations",
            "plot_files": ";".join(example_browser_outputs),
            "table_file": Path(args.example_browser_annotations).name,
            "description": "Gene-model rectangles, branchpoint markers, and text labels drawn in the annotation track.",
        },
    ]
    write_tsv(
        args.output_manifest,
        manifest_rows,
        ["plot_id", "panel_or_view", "plot_files", "table_file", "description"],
    )

    print(f"Wrote {len(manifest_rows)} plot source-data manifest rows to {args.output_manifest}")


if __name__ == "__main__":
    main()
