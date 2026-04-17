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


FEATURE_LABELS = {
    "anchored_enrichment_log2_fold_change": "Anchored abundance log2(DIS/ILS)",
    "control_anchored_fragments": "ILS anchored fragments",
    "query_anchored_fragments": "DIS anchored fragments",
    "control_gene_share": "ILS within-gene anchored share",
    "query_gene_share": "DIS within-gene anchored share",
    "gene_share_log2_fold_change": "Within-gene anchored share log2(DIS/ILS)",
    "gene_share_delta": "Within-gene anchored share delta (DIS - ILS)",
    "control_branch_fraction": "ILS branch fraction",
    "query_branch_fraction": "DIS branch fraction",
    "baseline_residual_query_branch_fraction": "DIS branch-fraction residual",
    "log_intron_length": "Intron length",
    "branchpoint_to_3ss_nt": "BP to 3'SS distance",
    "branchpoint_score": "BP score",
    "donor_maxent": "5'SS MaxEnt",
    "acceptor_maxent": "3'SS MaxEnt",
    "intron_gc": "Intron GC",
    "intron_gc_5p_window": "5' intron GC (100 nt)",
    "intron_gc_3p_window": "3' intron GC (100 nt)",
    "bp_to_3ss_pyrimidine_fraction": "BP to 3'SS pyrimidine fraction",
    "three_prime_window_pyrimidine_fraction": "3' window pyrimidine fraction",
}
COMPARISON_FIELDS = [
    "anchored_enrichment_log2_fold_change",
    "gene_share_log2_fold_change",
    "gene_share_delta",
    "control_gene_share",
    "query_gene_share",
    "control_branch_fraction",
    "query_branch_fraction",
    "baseline_residual_query_branch_fraction",
    "control_anchored_fragments",
    "query_anchored_fragments",
    "log_intron_length",
    "branchpoint_to_3ss_nt",
    "branchpoint_score",
    "donor_maxent",
    "acceptor_maxent",
    "intron_gc",
    "intron_gc_5p_window",
    "intron_gc_3p_window",
    "bp_to_3ss_pyrimidine_fraction",
    "three_prime_window_pyrimidine_fraction",
]
PINNED_PLOT_FIELDS = [
    "anchored_enrichment_log2_fold_change",
    "gene_share_log2_fold_change",
    "baseline_residual_query_branch_fraction",
    "control_branch_fraction",
    "query_branch_fraction",
]


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("--shared-introns", required=True)
    parser.add_argument("--site-counts", action="append", dest="site_counts", required=True)
    parser.add_argument("--control-condition", default="ILS")
    parser.add_argument("--query-condition", default="DIS")
    parser.add_argument("--control-enriched-quantile-fraction", type=float, default=0.1)
    parser.add_argument("--output-summary", required=True)
    parser.add_argument("--output-ranked-introns", required=True)
    parser.add_argument("--output-feature-group-comparison", required=True)
    parser.add_argument("--output-gene-summary", required=True)
    parser.add_argument("--output-feature-distribution-plot-png", required=True)
    parser.add_argument("--output-feature-distribution-plot-pdf", required=True)
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
            summary[condition]["total_anchored_fragments"] += anchored_fragments
    return summary


def classify_rows(rows, quantile_fraction):
    ordered = sorted(
        rows,
        key=lambda row: (
            float(row["anchored_enrichment_log2_fold_change"]),
            -int(row["control_anchored_fragments"]),
            row["gene_name"],
            row["intron_id"],
        ),
    )
    group_size = max(1, int(round(len(ordered) * quantile_fraction)))
    control_enriched_ids = {row["intron_id"] for row in ordered[:group_size]}
    cutoff = float(ordered[group_size - 1]["anchored_enrichment_log2_fold_change"])

    classified = []
    for rank, row in enumerate(ordered, start=1):
        control_enriched = row["intron_id"] in control_enriched_ids
        classified.append(
            {
                **row,
                "control_enrichment_rank": rank,
                "control_enriched_group": "control_enriched" if control_enriched else "other",
                "in_control_enriched_group": "1" if control_enriched else "0",
            }
        )
    return classified, cutoff


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


def build_output_rows(classified_rows, pooled_counts, control_condition, query_condition):
    output_rows = []
    for row in classified_rows:
        gene_id = row["gene_id"]
        if gene_id not in pooled_counts[control_condition]["genes"] or gene_id not in pooled_counts[query_condition]["genes"]:
            raise ValueError(f"Missing observed-gene totals for {gene_id}")
        control_gene_total = pooled_counts[control_condition]["genes"][gene_id]["anchored_fragments"]
        query_gene_total = pooled_counts[query_condition]["genes"][gene_id]["anchored_fragments"]
        output_rows.append(
            {
                **row,
                "control_gene_observed_anchored_fragments": control_gene_total,
                "query_gene_observed_anchored_fragments": query_gene_total,
                "control_gene_observed_introns": len(pooled_counts[control_condition]["genes"][gene_id]["introns"]),
                "query_gene_observed_introns": len(pooled_counts[query_condition]["genes"][gene_id]["introns"]),
                "control_gene_share": float(row["control_anchored_fragments"]) / control_gene_total,
                "query_gene_share": float(row["query_anchored_fragments"]) / query_gene_total,
                "gene_share_log2_fold_change": math.log2(
                    (float(row["query_anchored_fragments"]) / query_gene_total)
                    / (float(row["control_anchored_fragments"]) / control_gene_total)
                ),
                "gene_share_delta": (float(row["query_anchored_fragments"]) / query_gene_total)
                - (float(row["control_anchored_fragments"]) / control_gene_total),
            }
        )
    return output_rows


def build_group_comparison_rows(rows):
    control_enriched_rows = [row for row in rows if row["control_enriched_group"] == "control_enriched"]
    other_rows = [row for row in rows if row["control_enriched_group"] == "other"]
    comparison_rows = []
    for field in COMPARISON_FIELDS:
        control_values = np.array([float(row[field]) for row in control_enriched_rows], dtype=float)
        other_values = np.array([float(row[field]) for row in other_rows], dtype=float)
        comparison_rows.append(
            {
                "field": field,
                "label": FEATURE_LABELS.get(field, field),
                "control_enriched_introns": len(control_values),
                "other_introns": len(other_values),
                "control_enriched_mean": mean(control_values),
                "other_mean": mean(other_values),
                "delta_mean": mean(control_values) - mean(other_values),
                "control_enriched_median": median(control_values),
                "other_median": median(other_values),
                "delta_median": median(control_values) - median(other_values),
                "standardized_mean_difference": standardized_mean_difference(control_values, other_values),
            }
        )
    comparison_rows.sort(key=lambda row: abs(float(row["standardized_mean_difference"])), reverse=True)
    return comparison_rows


def select_plot_rows(comparison_rows, max_features=8):
    selected = []
    selected_fields = set()
    for field in PINNED_PLOT_FIELDS:
        matching = next((row for row in comparison_rows if row["field"] == field), None)
        if matching is None or field in selected_fields:
            continue
        selected.append(matching)
        selected_fields.add(field)
    for row in comparison_rows:
        if row["field"] in selected_fields:
            continue
        selected.append(row)
        selected_fields.add(row["field"])
        if len(selected) >= max_features:
            break
    return selected[:max_features]


def build_gene_summary_rows(rows):
    grouped = defaultdict(list)
    for row in rows:
        grouped[row["gene_id"]].append(row)

    gene_rows = []
    for gene_id, gene_introns in grouped.items():
        control_enriched_rows = [row for row in gene_introns if row["control_enriched_group"] == "control_enriched"]
        if not control_enriched_rows:
            continue
        gene_rows.append(
            {
                "gene_id": gene_id,
                "gene_name": gene_introns[0]["gene_name"],
                "shared_introns": len(gene_introns),
                "control_enriched_introns": len(control_enriched_rows),
                "other_introns": len(gene_introns) - len(control_enriched_rows),
                "control_enriched_fraction_of_gene": len(control_enriched_rows) / len(gene_introns),
                "mean_control_enriched_anchored_enrichment_log2_fold_change": mean(
                    [float(row["anchored_enrichment_log2_fold_change"]) for row in control_enriched_rows]
                ),
                "min_control_enriched_anchored_enrichment_log2_fold_change": min(
                    float(row["anchored_enrichment_log2_fold_change"]) for row in control_enriched_rows
                ),
                "mean_control_enriched_query_branch_fraction": mean(
                    [float(row["query_branch_fraction"]) for row in control_enriched_rows]
                ),
                "mean_control_enriched_control_branch_fraction": mean(
                    [float(row["control_branch_fraction"]) for row in control_enriched_rows]
                ),
            }
        )
    gene_rows.sort(
        key=lambda row: (
            -int(row["control_enriched_introns"]),
            float(row["min_control_enriched_anchored_enrichment_log2_fold_change"]),
            row["gene_name"],
        )
    )
    return gene_rows


def plot_feature_distributions(
    rows,
    comparison_rows,
    control_condition,
    query_condition,
    quantile_fraction,
    output_png,
    output_pdf,
):
    control_enriched_rows = [row for row in rows if row["control_enriched_group"] == "control_enriched"]
    other_rows = [row for row in rows if row["control_enriched_group"] == "other"]
    plot_rows = select_plot_rows(comparison_rows, max_features=8)
    n_features = len(plot_rows)
    ncols = 2
    nrows = int(math.ceil(n_features / ncols))
    figure, axes = plt.subplots(nrows, ncols, figsize=(12.5, 2.8 * nrows), constrained_layout=True)
    axes = np.atleast_1d(axes).ravel()

    for axis, row in zip(axes, plot_rows, strict=False):
        field = row["field"]
        control_values = np.array([float(entry[field]) for entry in control_enriched_rows], dtype=float)
        other_values = np.array([float(entry[field]) for entry in other_rows], dtype=float)
        boxplot = axis.boxplot(
            [other_values, control_values],
            vert=False,
            patch_artist=True,
            showfliers=False,
            widths=0.6,
            tick_labels=["Other", f"Bottom {quantile_fraction:.0%}"],
        )
        for patch, color in zip(boxplot["boxes"], ["#c7c7c7", "#1f77b4"], strict=True):
            patch.set_facecolor(color)
            patch.set_alpha(0.78)
        for median_line in boxplot["medians"]:
            median_line.set_color("black")
            median_line.set_linewidth(1.5)
        if field == "log_intron_length":
            axis.set_xlabel("log10 nt")
        elif "fraction" in field or field.startswith("intron_gc") or field.endswith("_share"):
            axis.set_xlim(0.0, 1.0)
            axis.set_xlabel("fraction")
        else:
            axis.set_xlabel("value")
        axis.set_title(row["label"])
        axis.text(
            0.98,
            0.06,
            f"SMD {float(row['standardized_mean_difference']):.3f}\n"
            f"delta median {float(row['delta_median']):.3g}",
            transform=axis.transAxes,
            ha="right",
            va="bottom",
        )
        axis.spines["top"].set_visible(False)
        axis.spines["right"].set_visible(False)

    for axis in axes[n_features:]:
        axis.axis("off")

    figure.suptitle(
        f"{control_condition}-enriched introns relative to {query_condition}\n"
        f"Bottom {quantile_fraction:.0%} of anchored abundance log2({query_condition}/{control_condition})",
        fontsize=13,
    )
    Path(output_png).parent.mkdir(parents=True, exist_ok=True)
    figure.savefig(output_png, dpi=300)
    figure.savefig(output_pdf)
    plt.close(figure)


def summarize_group_separation(comparison_rows):
    if not comparison_rows:
        return "", 0.0, 0.0
    best_row = max(comparison_rows, key=lambda row: abs(float(row["standardized_mean_difference"])))
    mean_abs_smd = float(np.mean([abs(float(row["standardized_mean_difference"])) for row in comparison_rows]))
    return best_row["field"], abs(float(best_row["standardized_mean_difference"])), mean_abs_smd


def main():
    args = parse_args()
    if not (0.0 < args.control_enriched_quantile_fraction <= 0.5):
        raise ValueError("--control-enriched-quantile-fraction must be in (0, 0.5]")

    shared_rows = read_tsv_rows(args.shared_introns)
    if not shared_rows:
        raise ValueError("No shared introns available")

    pooled_counts = pooled_anchored_counts(args.site_counts)
    classified_rows, cutoff = classify_rows(shared_rows, args.control_enriched_quantile_fraction)
    output_rows = build_output_rows(classified_rows, pooled_counts, args.control_condition, args.query_condition)
    comparison_rows = build_group_comparison_rows(output_rows)
    gene_rows = build_gene_summary_rows(output_rows)
    best_feature, best_abs_smd, mean_abs_smd = summarize_group_separation(comparison_rows)

    control_enriched_rows = [row for row in output_rows if row["control_enriched_group"] == "control_enriched"]
    other_rows = [row for row in output_rows if row["control_enriched_group"] == "other"]
    summary_rows = [
        {
            "control_condition": args.control_condition,
            "query_condition": args.query_condition,
            "control_enriched_quantile_fraction": args.control_enriched_quantile_fraction,
            "control_enriched_cutoff_anchored_enrichment_log2_fold_change": cutoff,
            "shared_introns": len(output_rows),
            "control_enriched_introns": len(control_enriched_rows),
            "other_introns": len(other_rows),
            "genes_with_control_enriched_introns": len(gene_rows),
            "genes_with_multiple_control_enriched_introns": sum(
                int(row["control_enriched_introns"]) > 1 for row in gene_rows
            ),
            "control_enriched_median_anchored_enrichment_log2_fold_change": median(
                [float(row["anchored_enrichment_log2_fold_change"]) for row in control_enriched_rows]
            ),
            "other_median_anchored_enrichment_log2_fold_change": median(
                [float(row["anchored_enrichment_log2_fold_change"]) for row in other_rows]
            ),
            "control_enriched_median_query_branch_fraction": median(
                [float(row["query_branch_fraction"]) for row in control_enriched_rows]
            ),
            "other_median_query_branch_fraction": median([float(row["query_branch_fraction"]) for row in other_rows]),
            "control_enriched_median_control_branch_fraction": median(
                [float(row["control_branch_fraction"]) for row in control_enriched_rows]
            ),
            "other_median_control_branch_fraction": median(
                [float(row["control_branch_fraction"]) for row in other_rows]
            ),
            "control_enriched_median_baseline_residual_query_branch_fraction": median(
                [float(row["baseline_residual_query_branch_fraction"]) for row in control_enriched_rows]
            ),
            "other_median_baseline_residual_query_branch_fraction": median(
                [float(row["baseline_residual_query_branch_fraction"]) for row in other_rows]
            ),
            "best_feature_by_abs_standardized_mean_difference": best_feature,
            "best_feature_abs_standardized_mean_difference": best_abs_smd,
            "mean_abs_standardized_mean_difference": mean_abs_smd,
        }
    ]

    output_rows.sort(key=lambda row: int(row["control_enrichment_rank"]))
    write_rows(args.output_summary, summary_rows, list(summary_rows[0].keys()))
    write_rows(args.output_ranked_introns, output_rows, list(output_rows[0].keys()))
    write_rows(args.output_feature_group_comparison, comparison_rows, list(comparison_rows[0].keys()))
    write_rows(args.output_gene_summary, gene_rows, list(gene_rows[0].keys()))
    plot_feature_distributions(
        output_rows,
        comparison_rows,
        args.control_condition,
        args.query_condition,
        args.control_enriched_quantile_fraction,
        args.output_feature_distribution_plot_png,
        args.output_feature_distribution_plot_pdf,
    )


if __name__ == "__main__":
    main()
