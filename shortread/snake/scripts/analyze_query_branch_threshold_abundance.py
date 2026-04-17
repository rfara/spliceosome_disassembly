#!/usr/bin/env python3

import argparse
import csv
import gzip
import math
from collections import defaultdict
from pathlib import Path

import matplotlib
import numpy as np

matplotlib.use("Agg")
import matplotlib.pyplot as plt


GROUP_COLORS = ["#c7c7c7", "#d95f02"]
GROUP_LABELS = {
    "other": "Other shared introns",
    "threshold_high": "High-branch DIS introns",
}
METRIC_LABELS = {
    "control_anchored_fragments": "ILS anchored fragments",
    "query_anchored_fragments": "DIS anchored fragments",
    "control_global_anchored_fraction": "ILS anchored fraction of pooled intron-end reads",
    "query_global_anchored_fraction": "DIS anchored fraction of pooled intron-end reads",
    "global_anchored_log2_fold_change": "Global anchored abundance log2(DIS/ILS)",
    "control_gene_share": "ILS within-gene anchored share",
    "query_gene_share": "DIS within-gene anchored share",
    "gene_share_log2_fold_change": "Within-gene anchored share log2(DIS/ILS)",
    "gene_share_delta": "Within-gene anchored share delta (DIS - ILS)",
}
GROUP_COMPARISON_FIELDS = [
    "global_anchored_log2_fold_change",
    "query_gene_share",
    "control_gene_share",
    "gene_share_log2_fold_change",
    "gene_share_delta",
    "query_anchored_fragments",
    "control_anchored_fragments",
]
PAIRED_METRICS = [
    ("delta_query_gene_share", "DIS within-gene share"),
    ("delta_global_anchored_log2_fold_change", "Global log2(DIS/ILS)"),
    ("delta_gene_share_log2_fold_change", "Within-gene share log2(DIS/ILS)"),
]


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("--shared-introns", required=True)
    parser.add_argument("--site-counts", action="append", dest="site_counts", required=True)
    parser.add_argument("--control-condition", default="ILS")
    parser.add_argument("--query-condition", default="DIS")
    parser.add_argument("--query-branch-threshold", type=float, default=0.25)
    parser.add_argument("--output-summary", required=True)
    parser.add_argument("--output-intron-abundance", required=True)
    parser.add_argument("--output-group-comparison", required=True)
    parser.add_argument("--output-gene-matched-comparison", required=True)
    parser.add_argument("--output-plot-png", required=True)
    parser.add_argument("--output-plot-pdf", required=True)
    return parser.parse_args()


def open_text(path, mode="rt"):
    if path.endswith(".gz"):
        return gzip.open(path, mode)
    return open(path, mode)


def read_tsv_rows(path):
    with open_text(path) as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_rows(path, rows, fieldnames):
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    with open_text(path, "wt") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def count_value(row, field):
    raw_value = row.get(field, 0)
    if raw_value in {"", None}:
        return 0
    return int(float(raw_value))


def pooled_anchored_counts(paths):
    summary = defaultdict(
        lambda: {
            "total_anchored_fragments": 0,
            "genes": defaultdict(lambda: {"anchored_fragments": 0, "introns": set(), "gene_name": ""}),
            "introns": {},
        }
    )
    for path in paths:
        rows = read_tsv_rows(path)
        if not rows:
            continue
        condition = rows[0]["condition"]
        for row in rows:
            anchored_fragments = count_value(row, "anchored_fragments")
            gene = summary[condition]["genes"][row["gene_id"]]
            gene["anchored_fragments"] += anchored_fragments
            gene["introns"].add(row["intron_id"])
            gene["gene_name"] = row["gene_name"]
            intron = summary[condition]["introns"].setdefault(
                row["intron_id"],
                {
                    "gene_id": row["gene_id"],
                    "gene_name": row["gene_name"],
                    "anchored_fragments": 0,
                },
            )
            intron["anchored_fragments"] += anchored_fragments
            summary[condition]["total_anchored_fragments"] += anchored_fragments
    return summary


def classify_shared_rows(rows, query_branch_threshold):
    ordered = sorted(
        rows,
        key=lambda row: (
            -float(row["query_branch_fraction"]),
            -int(row["query_anchored_fragments"]),
            row["gene_name"],
            row["intron_id"],
        ),
    )
    classified = []
    for rank, row in enumerate(ordered, start=1):
        high = float(row["query_branch_fraction"]) > query_branch_threshold
        classified.append(
            {
                **row,
                "query_branch_rank": rank,
                "threshold_group": "threshold_high" if high else "other",
                "meets_query_branch_threshold": "1" if high else "0",
            }
        )
    return classified


def mean(values):
    return float(np.mean(values)) if len(values) else math.nan


def median(values):
    return float(np.median(values)) if len(values) else math.nan


def standardized_mean_difference(a_values, b_values):
    a_values = np.array(a_values, dtype=float)
    b_values = np.array(b_values, dtype=float)
    if not len(a_values) or not len(b_values):
        return math.nan
    a_var = float(np.var(a_values, ddof=1)) if len(a_values) > 1 else 0.0
    b_var = float(np.var(b_values, ddof=1)) if len(b_values) > 1 else 0.0
    pooled_sd = math.sqrt(max((a_var + b_var) / 2.0, 0.0))
    if pooled_sd == 0.0:
        return 0.0
    return (mean(a_values) - mean(b_values)) / pooled_sd


def exact_sign_test_two_sided(values):
    nonzero = [value for value in values if value != 0.0]
    positives = int(sum(value > 0.0 for value in nonzero))
    negatives = int(sum(value < 0.0 for value in nonzero))
    trials = int(positives + negatives)
    if trials == 0:
        return 1.0, 0.0
    smaller_tail = min(positives, negatives)
    tail_probability = sum(math.comb(trials, k) for k in range(smaller_tail + 1)) / (2.0**trials)
    return min(1.0, 2.0 * tail_probability), float(positives / trials)


def build_intron_rows(classified_rows, pooled_counts, control_condition, query_condition):
    if control_condition not in pooled_counts or query_condition not in pooled_counts:
        raise ValueError("Missing pooled site-count rows for one or both requested conditions")
    control_total = pooled_counts[control_condition]["total_anchored_fragments"]
    query_total = pooled_counts[query_condition]["total_anchored_fragments"]
    if control_total <= 0 or query_total <= 0:
        raise ValueError("Pooled anchored fragment totals must be positive in both conditions")

    output_rows = []
    for row in classified_rows:
        gene_id = row["gene_id"]
        if gene_id not in pooled_counts[control_condition]["genes"] or gene_id not in pooled_counts[query_condition]["genes"]:
            raise ValueError(f"Missing observed-gene totals for {gene_id}")

        control_anchored = int(row["control_anchored_fragments"])
        query_anchored = int(row["query_anchored_fragments"])
        control_gene_total = pooled_counts[control_condition]["genes"][gene_id]["anchored_fragments"]
        query_gene_total = pooled_counts[query_condition]["genes"][gene_id]["anchored_fragments"]
        control_gene_introns = len(pooled_counts[control_condition]["genes"][gene_id]["introns"])
        query_gene_introns = len(pooled_counts[query_condition]["genes"][gene_id]["introns"])

        control_global_fraction = control_anchored / control_total
        query_global_fraction = query_anchored / query_total
        control_gene_share = control_anchored / control_gene_total
        query_gene_share = query_anchored / query_gene_total

        output_rows.append(
            {
                "query_branch_rank": row["query_branch_rank"],
                "threshold_group": row["threshold_group"],
                "meets_query_branch_threshold": row["meets_query_branch_threshold"],
                "intron_id": row["intron_id"],
                "gene_id": gene_id,
                "gene_name": row["gene_name"],
                "transcript_id": row["transcript_id"],
                "intron_number": row["intron_number"],
                "control_anchored_fragments": control_anchored,
                "query_anchored_fragments": query_anchored,
                "control_branch_fraction": float(row["control_branch_fraction"]),
                "query_branch_fraction": float(row["query_branch_fraction"]),
                "control_global_anchored_fraction": control_global_fraction,
                "query_global_anchored_fraction": query_global_fraction,
                "global_anchored_log2_fold_change": math.log2(query_global_fraction / control_global_fraction),
                "control_gene_observed_anchored_fragments": control_gene_total,
                "query_gene_observed_anchored_fragments": query_gene_total,
                "control_gene_observed_introns": control_gene_introns,
                "query_gene_observed_introns": query_gene_introns,
                "control_gene_share": control_gene_share,
                "query_gene_share": query_gene_share,
                "gene_share_log2_fold_change": math.log2(query_gene_share / control_gene_share),
                "gene_share_delta": query_gene_share - control_gene_share,
            }
        )
    return output_rows


def build_group_comparison_rows(rows):
    high_rows = [row for row in rows if row["threshold_group"] == "threshold_high"]
    other_rows = [row for row in rows if row["threshold_group"] == "other"]
    comparison_rows = []
    for field in GROUP_COMPARISON_FIELDS:
        high_values = np.array([float(row[field]) for row in high_rows], dtype=float)
        other_values = np.array([float(row[field]) for row in other_rows], dtype=float)
        comparison_rows.append(
            {
                "field": field,
                "label": METRIC_LABELS.get(field, field),
                "threshold_high_n": len(high_values),
                "other_n": len(other_values),
                "threshold_high_mean": mean(high_values),
                "other_mean": mean(other_values),
                "delta_mean": mean(high_values) - mean(other_values),
                "threshold_high_median": median(high_values),
                "other_median": median(other_values),
                "delta_median": median(high_values) - median(other_values),
                "standardized_mean_difference": standardized_mean_difference(high_values, other_values),
            }
        )
    comparison_rows.sort(key=lambda row: abs(float(row["standardized_mean_difference"])), reverse=True)
    return comparison_rows


def build_gene_matched_rows(rows):
    grouped = defaultdict(list)
    for row in rows:
        grouped[row["gene_id"]].append(row)

    output_rows = []
    for gene_id, gene_rows in grouped.items():
        high_rows = [row for row in gene_rows if row["threshold_group"] == "threshold_high"]
        other_rows = [row for row in gene_rows if row["threshold_group"] == "other"]
        if not high_rows or not other_rows:
            continue
        output_rows.append(
            {
                "gene_id": gene_id,
                "gene_name": gene_rows[0]["gene_name"],
                "threshold_high_introns": len(high_rows),
                "other_introns": len(other_rows),
                "mean_threshold_high_query_gene_share": mean([float(row["query_gene_share"]) for row in high_rows]),
                "mean_other_query_gene_share": mean([float(row["query_gene_share"]) for row in other_rows]),
                "delta_query_gene_share": mean([float(row["query_gene_share"]) for row in high_rows])
                - mean([float(row["query_gene_share"]) for row in other_rows]),
                "mean_threshold_high_global_anchored_log2_fold_change": mean(
                    [float(row["global_anchored_log2_fold_change"]) for row in high_rows]
                ),
                "mean_other_global_anchored_log2_fold_change": mean(
                    [float(row["global_anchored_log2_fold_change"]) for row in other_rows]
                ),
                "delta_global_anchored_log2_fold_change": mean(
                    [float(row["global_anchored_log2_fold_change"]) for row in high_rows]
                )
                - mean([float(row["global_anchored_log2_fold_change"]) for row in other_rows]),
                "mean_threshold_high_gene_share_log2_fold_change": mean(
                    [float(row["gene_share_log2_fold_change"]) for row in high_rows]
                ),
                "mean_other_gene_share_log2_fold_change": mean(
                    [float(row["gene_share_log2_fold_change"]) for row in other_rows]
                ),
                "delta_gene_share_log2_fold_change": mean(
                    [float(row["gene_share_log2_fold_change"]) for row in high_rows]
                )
                - mean([float(row["gene_share_log2_fold_change"]) for row in other_rows]),
            }
        )
    output_rows.sort(key=lambda row: (float(row["delta_gene_share_log2_fold_change"]), row["gene_name"]))
    return output_rows


def paired_metric_summary(rows, field):
    values = np.array([float(row[field]) for row in rows], dtype=float)
    pvalue, fraction_positive = exact_sign_test_two_sided(values)
    return {
        "n": len(values),
        "mean": mean(values),
        "median": median(values),
        "fraction_positive": fraction_positive,
        "sign_test_pvalue": pvalue,
    }


def style_group_boxplot(axis, other_values, high_values, label, query_condition, query_branch_threshold, note):
    boxplot = axis.boxplot(
        [other_values, high_values],
        vert=False,
        patch_artist=True,
        showfliers=False,
        widths=0.6,
        tick_labels=["Other", f"{query_condition} > {query_branch_threshold:.2f}"],
    )
    for patch, color in zip(boxplot["boxes"], GROUP_COLORS, strict=True):
        patch.set_facecolor(color)
        patch.set_alpha(0.8)
    for median_line in boxplot["medians"]:
        median_line.set_color("black")
        median_line.set_linewidth(1.5)
    axis.set_title(label)
    axis.text(0.98, 0.06, note, transform=axis.transAxes, ha="right", va="bottom")
    axis.spines["top"].set_visible(False)
    axis.spines["right"].set_visible(False)


def style_single_distribution_boxplot(axis, values, label, note):
    boxplot = axis.boxplot(
        [values],
        vert=False,
        patch_artist=True,
        showfliers=False,
        widths=0.6,
        tick_labels=["Genes"],
    )
    boxplot["boxes"][0].set_facecolor("#7570b3")
    boxplot["boxes"][0].set_alpha(0.8)
    boxplot["medians"][0].set_color("black")
    boxplot["medians"][0].set_linewidth(1.5)
    axis.axvline(0.0, color="#666666", linestyle="--", linewidth=1.0)
    axis.set_title(label)
    axis.text(0.98, 0.06, note, transform=axis.transAxes, ha="right", va="bottom")
    axis.spines["top"].set_visible(False)
    axis.spines["right"].set_visible(False)


def plot_abundance_panels(
    intron_rows,
    gene_rows,
    query_condition,
    query_branch_threshold,
    output_png,
    output_pdf,
):
    high_rows = [row for row in intron_rows if row["threshold_group"] == "threshold_high"]
    other_rows = [row for row in intron_rows if row["threshold_group"] == "other"]

    figure, axes = plt.subplots(2, 3, figsize=(15, 8.5), constrained_layout=True)
    axes = axes.ravel()

    group_specs = [
        ("global_anchored_log2_fold_change", "Global anchored abundance log2(DIS/ILS)"),
        ("query_gene_share", f"{query_condition} within-gene anchored share"),
        ("gene_share_log2_fold_change", "Within-gene anchored share log2(DIS/ILS)"),
    ]
    for axis, (field, label) in zip(axes[:3], group_specs, strict=True):
        high_values = np.array([float(row[field]) for row in high_rows], dtype=float)
        other_values = np.array([float(row[field]) for row in other_rows], dtype=float)
        note = (
            f"delta median {median(high_values) - median(other_values):.3g}\n"
            f"SMD {standardized_mean_difference(high_values, other_values):.3f}"
        )
        style_group_boxplot(axis, other_values, high_values, label, query_condition, query_branch_threshold, note)
        if "share" in field and "log2" not in field:
            axis.set_xlim(0.0, 1.0)

    for axis, (field, label) in zip(axes[3:], PAIRED_METRICS, strict=True):
        values = np.array([float(row[field]) for row in gene_rows], dtype=float)
        pvalue, fraction_positive = exact_sign_test_two_sided(values)
        note = f"n genes {len(values)}\nfrac > 0 {fraction_positive:.3f}\nsign p {pvalue:.2g}"
        style_single_distribution_boxplot(axis, values, f"Gene-matched delta: {label}", note)
        if field == "delta_query_gene_share":
            axis.set_xlim(-1.0, 1.0)

    figure.suptitle(
        f"{query_condition} high-branch intron abundance\n"
        f"High group = pooled {query_condition} branch fraction > {query_branch_threshold:.2f}\n"
        "Within-gene denominators use introns observed in each condition",
        fontsize=13,
    )
    Path(output_png).parent.mkdir(parents=True, exist_ok=True)
    figure.savefig(output_png, dpi=300)
    figure.savefig(output_pdf)
    plt.close(figure)


def main():
    args = parse_args()
    if not (0.0 < args.query_branch_threshold < 1.0):
        raise ValueError("--query-branch-threshold must be between 0 and 1")

    shared_rows = read_tsv_rows(args.shared_introns)
    if not shared_rows:
        raise ValueError("No shared introns available")

    pooled_counts = pooled_anchored_counts(args.site_counts)
    classified_rows = classify_shared_rows(shared_rows, args.query_branch_threshold)
    intron_rows = build_intron_rows(classified_rows, pooled_counts, args.control_condition, args.query_condition)
    intron_rows.sort(key=lambda row: int(row["query_branch_rank"]))

    group_rows = build_group_comparison_rows(intron_rows)
    gene_rows = build_gene_matched_rows(intron_rows)
    paired_query_gene_share = paired_metric_summary(gene_rows, "delta_query_gene_share")
    paired_global_log2fc = paired_metric_summary(gene_rows, "delta_global_anchored_log2_fold_change")
    paired_gene_share_log2fc = paired_metric_summary(gene_rows, "delta_gene_share_log2_fold_change")

    threshold_high_rows = [row for row in intron_rows if row["threshold_group"] == "threshold_high"]
    other_rows = [row for row in intron_rows if row["threshold_group"] == "other"]
    summary_rows = [
        {
            "query_condition": args.query_condition,
            "control_condition": args.control_condition,
            "query_branch_threshold": args.query_branch_threshold,
            "shared_introns": len(intron_rows),
            "threshold_high_introns": len(threshold_high_rows),
            "other_introns": len(other_rows),
            "query_total_observed_anchored_fragments": pooled_counts[args.query_condition]["total_anchored_fragments"],
            "control_total_observed_anchored_fragments": pooled_counts[args.control_condition]["total_anchored_fragments"],
            "query_observed_introns": len(pooled_counts[args.query_condition]["introns"]),
            "control_observed_introns": len(pooled_counts[args.control_condition]["introns"]),
            "query_observed_genes": len(pooled_counts[args.query_condition]["genes"]),
            "control_observed_genes": len(pooled_counts[args.control_condition]["genes"]),
            "threshold_high_global_anchored_log2fc_mean": mean(
                [float(row["global_anchored_log2_fold_change"]) for row in threshold_high_rows]
            ),
            "threshold_high_global_anchored_log2fc_median": median(
                [float(row["global_anchored_log2_fold_change"]) for row in threshold_high_rows]
            ),
            "other_global_anchored_log2fc_mean": mean(
                [float(row["global_anchored_log2_fold_change"]) for row in other_rows]
            ),
            "other_global_anchored_log2fc_median": median(
                [float(row["global_anchored_log2_fold_change"]) for row in other_rows]
            ),
            "threshold_high_query_gene_share_mean": mean([float(row["query_gene_share"]) for row in threshold_high_rows]),
            "threshold_high_query_gene_share_median": median(
                [float(row["query_gene_share"]) for row in threshold_high_rows]
            ),
            "other_query_gene_share_mean": mean([float(row["query_gene_share"]) for row in other_rows]),
            "other_query_gene_share_median": median([float(row["query_gene_share"]) for row in other_rows]),
            "threshold_high_gene_share_log2fc_mean": mean(
                [float(row["gene_share_log2_fold_change"]) for row in threshold_high_rows]
            ),
            "threshold_high_gene_share_log2fc_median": median(
                [float(row["gene_share_log2_fold_change"]) for row in threshold_high_rows]
            ),
            "other_gene_share_log2fc_mean": mean(
                [float(row["gene_share_log2_fold_change"]) for row in other_rows]
            ),
            "other_gene_share_log2fc_median": median(
                [float(row["gene_share_log2_fold_change"]) for row in other_rows]
            ),
            "genes_with_both_groups": len(gene_rows),
            "paired_query_gene_share_mean_delta": paired_query_gene_share["mean"],
            "paired_query_gene_share_median_delta": paired_query_gene_share["median"],
            "paired_query_gene_share_fraction_positive": paired_query_gene_share["fraction_positive"],
            "paired_query_gene_share_sign_test_pvalue": paired_query_gene_share["sign_test_pvalue"],
            "paired_global_anchored_log2fc_mean_delta": paired_global_log2fc["mean"],
            "paired_global_anchored_log2fc_median_delta": paired_global_log2fc["median"],
            "paired_global_anchored_log2fc_fraction_positive": paired_global_log2fc["fraction_positive"],
            "paired_global_anchored_log2fc_sign_test_pvalue": paired_global_log2fc["sign_test_pvalue"],
            "paired_gene_share_log2fc_mean_delta": paired_gene_share_log2fc["mean"],
            "paired_gene_share_log2fc_median_delta": paired_gene_share_log2fc["median"],
            "paired_gene_share_log2fc_fraction_positive": paired_gene_share_log2fc["fraction_positive"],
            "paired_gene_share_log2fc_sign_test_pvalue": paired_gene_share_log2fc["sign_test_pvalue"],
        }
    ]

    write_rows(args.output_summary, summary_rows, list(summary_rows[0].keys()))
    write_rows(args.output_intron_abundance, intron_rows, list(intron_rows[0].keys()))
    write_rows(args.output_group_comparison, group_rows, list(group_rows[0].keys()))
    write_rows(args.output_gene_matched_comparison, gene_rows, list(gene_rows[0].keys()))
    plot_abundance_panels(
        intron_rows,
        gene_rows,
        args.query_condition,
        args.query_branch_threshold,
        args.output_plot_png,
        args.output_plot_pdf,
    )


if __name__ == "__main__":
    main()
