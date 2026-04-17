#!/usr/bin/env python3

import argparse
import csv
import gzip
from collections import defaultdict
from pathlib import Path

import matplotlib
import numpy as np

matplotlib.use("Agg")
import matplotlib.pyplot as plt


CONDITION_COLORS = {
    "ILS": "#1f77b4",
    "DIS": "#d95f02",
}

METADATA_FIELDS = [
    "gene_id",
    "gene_name",
    "transcript_id",
    "intron_number",
    "chrom",
    "strand",
    "intron_start",
    "intron_end",
    "intron_length",
    "three_prime_ss",
    "branchpoint_position",
    "branchpoint_score",
    "branchpoint_to_3ss_nt",
    "branchpoint_type",
    "branchpoint_candidates",
]


def parse_args():
    parser = argparse.ArgumentParser(
        description=(
            "Measure how often introns observed with anchored reads in one condition "
            "also have anchored-read coverage in another condition across percentile cutoffs."
        )
    )
    parser.add_argument("--site-counts", action="append", dest="site_counts", required=True)
    parser.add_argument("--reference-condition", default="ILS")
    parser.add_argument("--query-condition", default="DIS")
    parser.add_argument("--percentile-start", type=int, default=1)
    parser.add_argument("--percentile-end", type=int, default=99)
    parser.add_argument("--output-summary", required=True)
    parser.add_argument("--output-intron-counts", required=True)
    parser.add_argument("--output-plot-png", required=True)
    parser.add_argument("--output-plot-pdf", required=True)
    return parser.parse_args()


def open_text(path, mode="rt"):
    if str(path).endswith(".gz"):
        return gzip.open(path, mode)
    return open(path, mode)


def count_value(row, field):
    raw_value = row.get(field, 0)
    if raw_value in {"", None}:
        return 0
    return int(float(raw_value))


def infer_sample_name_from_path(path):
    return Path(path).name.split(".")[0]


def read_site_counts(paths):
    condition_counts = defaultdict(lambda: defaultdict(int))
    condition_sample_counts = defaultdict(lambda: defaultdict(dict))
    condition_samples = defaultdict(list)
    metadata = {}

    for path in paths:
        path_sample = infer_sample_name_from_path(path)
        file_samples = set()
        file_conditions = set()

        with open_text(path) as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            for row in reader:
                sample = row.get("sample") or path_sample
                condition = row["condition"]
                intron_id = row["intron_id"]
                anchored_fragments = count_value(row, "anchored_fragments")

                file_samples.add(sample)
                file_conditions.add(condition)
                condition_counts[condition][intron_id] += anchored_fragments
                condition_sample_counts[condition][intron_id][sample] = anchored_fragments
                metadata.setdefault(intron_id, {field: row.get(field, "") for field in METADATA_FIELDS})

        if len(file_samples) != 1:
            raise ValueError(f"Expected one sample in {path}, found {sorted(file_samples)}")
        if len(file_conditions) != 1:
            raise ValueError(f"Expected one condition in {path}, found {sorted(file_conditions)}")

        sample = next(iter(file_samples))
        condition = next(iter(file_conditions))
        if sample not in condition_samples[condition]:
            condition_samples[condition].append(sample)

    return condition_counts, condition_sample_counts, condition_samples, metadata


def write_rows(path, rows, fieldnames):
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    with open_text(path, "wt") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def build_intron_rows(
    reference_condition,
    query_condition,
    reference_counts,
    query_counts,
    query_sample_counts,
    metadata,
):
    rows = []
    for intron_id, reference_count in reference_counts.items():
        if reference_count <= 0:
            continue
        query_count = query_counts.get(intron_id, 0)
        query_present_samples = sum(1 for value in query_sample_counts.get(intron_id, {}).values() if value > 0)
        row = {
            "intron_id": intron_id,
            **metadata.get(intron_id, {field: "" for field in METADATA_FIELDS}),
            "reference_condition": reference_condition,
            "query_condition": query_condition,
            "reference_anchored_fragments": reference_count,
            "query_anchored_fragments": query_count,
            "query_samples_with_anchored_fragments": query_present_samples,
            "query_has_any_anchored_fragments": int(query_present_samples > 0),
        }
        rows.append(row)

    rows.sort(
        key=lambda row: (
            -row["reference_anchored_fragments"],
            row.get("gene_name", ""),
            row["intron_id"],
        )
    )
    return rows


def build_percentile_rows(
    intron_rows,
    reference_condition,
    query_condition,
    percentile_start,
    percentile_end,
    query_to_reference_depth_ratio,
    random_thinning_fraction,
):
    reference_values = np.array(
        [row["reference_anchored_fragments"] for row in intron_rows if row["reference_anchored_fragments"] > 0],
        dtype=np.float64,
    )
    if reference_values.size == 0:
        raise ValueError(f"No introns had anchored fragments in {reference_condition}")
    if percentile_start < 0 or percentile_end > 100 or percentile_start > percentile_end:
        raise ValueError("Percentile range must satisfy 0 <= start <= end <= 100")

    rows = []
    for percentile in range(percentile_start, percentile_end + 1):
        cutoff = float(np.percentile(reference_values, percentile, method="linear"))
        selected = [row for row in intron_rows if row["reference_anchored_fragments"] >= cutoff]
        covered = [row for row in selected if row["query_has_any_anchored_fragments"]]
        selected_count = len(selected)
        covered_count = len(covered)
        selected_reference_fragments = sum(row["reference_anchored_fragments"] for row in selected)
        selected_query_fragments = sum(row["query_anchored_fragments"] for row in selected)
        random_thinning_expected_covered = sum(
            1.0 - ((1.0 - random_thinning_fraction) ** row["reference_anchored_fragments"])
            for row in selected
        )

        rows.append(
            {
                "reference_condition": reference_condition,
                "query_condition": query_condition,
                "reference_percentile_cutoff": percentile,
                "reference_anchored_fragment_cutoff": cutoff,
                "reference_introns_at_or_above_cutoff": selected_count,
                "query_introns_with_any_anchored_fragments": covered_count,
                "query_introns_without_anchored_fragments": selected_count - covered_count,
                "query_covered_percent": 0.0 if selected_count == 0 else covered_count * 100.0 / selected_count,
                "query_to_reference_anchored_fragment_ratio": query_to_reference_depth_ratio,
                "random_thinning_fraction": random_thinning_fraction,
                "random_thinning_expected_introns_with_anchored_fragments": random_thinning_expected_covered,
                "random_thinning_expected_covered_percent": 0.0
                if selected_count == 0
                else random_thinning_expected_covered * 100.0 / selected_count,
                "reference_anchored_fragments_at_or_above_cutoff": selected_reference_fragments,
                "query_anchored_fragments_in_reference_introns_at_or_above_cutoff": selected_query_fragments,
                "total_reference_anchor_positive_introns": int(reference_values.size),
            }
        )
    return rows


def plot_percentile_summary(rows, reference_condition, query_condition, output_png, output_pdf):
    Path(output_png).parent.mkdir(parents=True, exist_ok=True)
    Path(output_pdf).parent.mkdir(parents=True, exist_ok=True)
    x_values = [int(row["reference_percentile_cutoff"]) for row in rows]
    y_values = [float(row["query_covered_percent"]) for row in rows]
    random_thinning_values = [float(row["random_thinning_expected_covered_percent"]) for row in rows]
    random_thinning_fraction = float(rows[0]["random_thinning_fraction"])
    query_to_reference_depth_ratio = float(rows[0]["query_to_reference_anchored_fragment_ratio"])
    include_depth_model = query_to_reference_depth_ratio <= 1.0

    fig, axis = plt.subplots(figsize=(7.2, 4.8))
    color = CONDITION_COLORS.get(query_condition, "#4c4c4c")
    axis.plot(x_values, y_values, color=color, linewidth=2.5, label=f"Observed {query_condition}")
    axis.scatter(x_values, y_values, color=color, s=14, linewidths=0)
    if include_depth_model:
        axis.plot(
            x_values,
            random_thinning_values,
            color="#4c4c4c",
            linewidth=2.2,
            linestyle="--",
            label=f"{reference_condition} downsampled to {query_condition} depth",
        )
    axis.set_xlim(min(x_values), max(x_values))
    axis.set_ylim(0, 100)
    axis.set_xlabel(f"{reference_condition} anchored-count percentile cutoff")
    axis.set_ylabel(f"{query_condition} introns with anchored coverage (%)")
    axis.set_title(f"{query_condition} coverage of {reference_condition}-anchored introns")
    axis.grid(axis="y", color="#d9d9d9", linewidth=0.8)
    axis.spines["top"].set_visible(False)
    axis.spines["right"].set_visible(False)
    axis.legend(loc="center left", bbox_to_anchor=(1.02, 0.5), frameon=False, fontsize=9)

    annotation = f"Pooled counts; {query_condition} covered: >=1 anchored fragment in any replicate"
    if include_depth_model:
        annotation = (
            f"{annotation}\n"
            f"Downsampling q={random_thinning_fraction:.4f} "
            f"({query_condition}/{reference_condition} anchored fragments)"
        )
    axis.text(
        0.02,
        0.04,
        annotation,
        transform=axis.transAxes,
        ha="left",
        va="bottom",
        fontsize=9,
        color="#333333",
    )

    fig.tight_layout()
    fig.savefig(output_png, dpi=300, bbox_inches="tight")
    fig.savefig(output_pdf, bbox_inches="tight")
    plt.close(fig)


def main():
    args = parse_args()
    condition_counts, condition_sample_counts, condition_samples, metadata = read_site_counts(args.site_counts)

    if args.reference_condition not in condition_counts:
        raise ValueError(f"Reference condition not found: {args.reference_condition}")
    if args.query_condition not in condition_counts:
        raise ValueError(f"Query condition not found: {args.query_condition}")

    reference_counts = condition_counts[args.reference_condition]
    query_counts = condition_counts[args.query_condition]
    query_sample_counts = condition_sample_counts[args.query_condition]
    reference_anchored_fragments = sum(reference_counts.values())
    query_anchored_fragments = sum(query_counts.values())
    query_to_reference_depth_ratio = 0.0 if reference_anchored_fragments == 0 else (
        query_anchored_fragments / reference_anchored_fragments
    )
    random_thinning_fraction = min(1.0, query_to_reference_depth_ratio)

    intron_rows = build_intron_rows(
        args.reference_condition,
        args.query_condition,
        reference_counts,
        query_counts,
        query_sample_counts,
        metadata,
    )
    percentile_rows = build_percentile_rows(
        intron_rows,
        args.reference_condition,
        args.query_condition,
        args.percentile_start,
        args.percentile_end,
        query_to_reference_depth_ratio,
        random_thinning_fraction,
    )

    intron_fieldnames = [
        "intron_id",
        *METADATA_FIELDS,
        "reference_condition",
        "query_condition",
        "reference_anchored_fragments",
        "query_anchored_fragments",
        "query_samples_with_anchored_fragments",
        "query_has_any_anchored_fragments",
    ]
    summary_fieldnames = [
        "reference_condition",
        "query_condition",
        "reference_percentile_cutoff",
        "reference_anchored_fragment_cutoff",
        "reference_introns_at_or_above_cutoff",
        "query_introns_with_any_anchored_fragments",
        "query_introns_without_anchored_fragments",
        "query_covered_percent",
        "query_to_reference_anchored_fragment_ratio",
        "random_thinning_fraction",
        "random_thinning_expected_introns_with_anchored_fragments",
        "random_thinning_expected_covered_percent",
        "reference_anchored_fragments_at_or_above_cutoff",
        "query_anchored_fragments_in_reference_introns_at_or_above_cutoff",
        "total_reference_anchor_positive_introns",
    ]

    write_rows(args.output_intron_counts, intron_rows, intron_fieldnames)
    write_rows(args.output_summary, percentile_rows, summary_fieldnames)
    plot_percentile_summary(
        percentile_rows,
        args.reference_condition,
        args.query_condition,
        args.output_plot_png,
        args.output_plot_pdf,
    )

    print(f"Reference condition: {args.reference_condition}")
    print(f"Query condition: {args.query_condition}")
    print(f"Reference samples: {', '.join(condition_samples[args.reference_condition])}")
    print(f"Query samples: {', '.join(condition_samples[args.query_condition])}")
    print(f"Reference anchor-positive introns: {len(intron_rows)}")
    print(f"Query/reference anchored-fragment ratio: {query_to_reference_depth_ratio}")
    print(f"Depth-model per-fragment detection probability: {random_thinning_fraction}")
    print(f"Percentile cutoffs written: {len(percentile_rows)}")


if __name__ == "__main__":
    main()
