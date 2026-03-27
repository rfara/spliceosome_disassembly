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


CONDITION_COLORS = {
    "ILS": "#1f77b4",
    "DIS": "#d95f02",
}

PANEL_GRID = np.linspace(0.0, 1.0, 101)


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("--site-counts", action="append", dest="site_counts", required=True)
    parser.add_argument("--signal-field", default="zero_or_plus_one_branchpoint_fragments")
    parser.add_argument("--sample-min-reads", type=int, default=10)
    parser.add_argument("--condition-min-reads", type=int, default=20)
    parser.add_argument("--bootstraps", type=int, default=500)
    parser.add_argument("--seed", type=int, default=1)
    parser.add_argument("--query-condition", default="DIS")
    parser.add_argument("--control-condition", default="ILS")
    parser.add_argument("--output-summary-by-sample", required=True)
    parser.add_argument("--output-summary-by-condition", required=True)
    parser.add_argument("--output-condition-introns", required=True)
    parser.add_argument("--output-condition-comparison-summary", required=True)
    parser.add_argument("--output-condition-comparison-introns", required=True)
    parser.add_argument("--output-control-deciles", required=True)
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


def infer_sample_name(rows, path):
    samples = {row["sample"] for row in rows if row.get("sample")}
    if not samples:
        return path.rsplit("/", 1)[-1].split(".")[0]
    if len(samples) != 1:
        raise ValueError(f"Expected one sample in {path}, found {sorted(samples)}")
    return next(iter(samples))


def signal_label(signal_field):
    labels = {
        "exact_branchpoint_fragments": "exact branchpoint-terminated fragments",
        "plus_one_branchpoint_fragments": "+1 branchpoint-proximal fragments",
        "zero_or_plus_one_branchpoint_fragments": "0/+1 branchpoint-proximal fragments",
    }
    return labels.get(signal_field, signal_field.replace("_", " "))


def build_record(row, signal_field):
    anchored_fragments = count_value(row, "anchored_fragments")
    branched_fragments = count_value(row, signal_field)
    branch_fraction = 0.0 if anchored_fragments == 0 else branched_fragments / anchored_fragments
    return {
        "intron_id": row["intron_id"],
        "gene_id": row["gene_id"],
        "gene_name": row["gene_name"],
        "transcript_id": row["transcript_id"],
        "intron_number": row["intron_number"],
        "chrom": row["chrom"],
        "strand": row["strand"],
        "intron_start": row["intron_start"],
        "intron_end": row["intron_end"],
        "three_prime_ss": row["three_prime_ss"],
        "branchpoint_position": row["branchpoint_position"],
        "branchpoint_score": row["branchpoint_score"],
        "branchpoint_to_3ss_nt": row["branchpoint_to_3ss_nt"],
        "branchpoint_type": row.get("branchpoint_type", ""),
        "branchpoint_candidates": row["branchpoint_candidates"],
        "anchored_fragments": anchored_fragments,
        "branched_fragments": branched_fragments,
        "branch_fraction": branch_fraction,
    }


def load_datasets(paths, signal_field):
    sample_rows = {}
    condition_rows = defaultdict(dict)

    for path in paths:
        rows = read_tsv_rows(path)
        sample = infer_sample_name(rows, path)
        if not rows:
            continue
        condition = rows[0]["condition"]
        sample_records = {}
        for row in rows:
            record = build_record(row, signal_field)
            intron_id = record["intron_id"]
            sample_records[intron_id] = record

            pooled = condition_rows[condition].setdefault(
                intron_id,
                {
                    key: record[key]
                    for key in record
                    if key not in {"anchored_fragments", "branched_fragments", "branch_fraction"}
                }
                | {"anchored_fragments": 0, "branched_fragments": 0, "branch_fraction": 0.0},
            )
            pooled["anchored_fragments"] += record["anchored_fragments"]
            pooled["branched_fragments"] += record["branched_fragments"]
        sample_rows[sample] = {
            "condition": condition,
            "records": sample_records,
        }

    for condition in condition_rows:
        for record in condition_rows[condition].values():
            anchored = record["anchored_fragments"]
            record["branch_fraction"] = 0.0 if anchored == 0 else record["branched_fragments"] / anchored

    return sample_rows, condition_rows


def filtered_records(records_by_intron, min_reads):
    return [record for record in records_by_intron.values() if record["anchored_fragments"] >= min_reads]


def pearson_chisq(counts, depths, probability):
    if counts.size == 0 or not (0.0 < probability < 1.0):
        return 0.0
    expected = depths * probability
    denominator = expected * (1.0 - probability)
    return float(np.sum((counts - expected) ** 2 / denominator))


def bootstrap_homogeneous_null(depths, probability, observed_chisq, bootstraps, seed, quantile_grid=None):
    if depths.size == 0 or not (0.0 < probability < 1.0) or bootstraps <= 0:
        zero_count = 0 if probability >= 1.0 else int(depths.size)
        quantiles = np.zeros(len(quantile_grid)) if quantile_grid is not None else None
        return {
            "bootstrap_pvalue": 1.0,
            "null_zero_branch_introns_mean": float(zero_count),
            "null_zero_branch_introns_q025": float(zero_count),
            "null_zero_branch_introns_q975": float(zero_count),
            "null_quantiles_lower": quantiles,
            "null_quantiles_upper": quantiles,
        }

    rng = np.random.default_rng(seed)
    expected = depths * probability
    denominator = expected * (1.0 - probability)
    zero_counts = np.empty(bootstraps, dtype=np.int32)
    null_quantiles = (
        np.empty((bootstraps, len(quantile_grid)), dtype=np.float64) if quantile_grid is not None else None
    )

    exceedances = 0
    batch_size = 32
    completed = 0
    while completed < bootstraps:
        batch = min(batch_size, bootstraps - completed)
        simulated = rng.binomial(depths[None, :], probability, size=(batch, depths.size))
        diff = simulated - expected
        batch_chisq = np.sum((diff * diff) / denominator, axis=1)
        exceedances += int(np.sum(batch_chisq >= observed_chisq))
        zero_counts[completed : completed + batch] = np.sum(simulated == 0, axis=1)

        if quantile_grid is not None:
            simulated_fractions = simulated / depths
            null_quantiles[completed : completed + batch, :] = np.quantile(
                simulated_fractions,
                quantile_grid,
                axis=1,
            ).T

        completed += batch

    return {
        "bootstrap_pvalue": (exceedances + 1.0) / (bootstraps + 1.0),
        "null_zero_branch_introns_mean": float(np.mean(zero_counts)),
        "null_zero_branch_introns_q025": float(np.quantile(zero_counts, 0.025)),
        "null_zero_branch_introns_q975": float(np.quantile(zero_counts, 0.975)),
        "null_quantiles_lower": None
        if null_quantiles is None
        else np.quantile(null_quantiles, 0.025, axis=0),
        "null_quantiles_upper": None
        if null_quantiles is None
        else np.quantile(null_quantiles, 0.975, axis=0),
    }


def analyse_dataset(
    records,
    dataset_name,
    dataset_type,
    condition,
    min_reads,
    signal_field,
    bootstraps,
    seed,
    quantile_grid=None,
):
    filtered = filtered_records(records, min_reads)
    filtered.sort(key=lambda record: (-record["branch_fraction"], -record["anchored_fragments"], record["gene_name"]))

    depths = np.array([record["anchored_fragments"] for record in filtered], dtype=np.int64)
    counts = np.array([record["branched_fragments"] for record in filtered], dtype=np.int64)
    fractions = np.array([record["branch_fraction"] for record in filtered], dtype=np.float64)

    overall_probability = 0.0 if depths.sum() == 0 else float(counts.sum() / depths.sum())
    observed_chisq = pearson_chisq(counts, depths, overall_probability)
    df = max(len(filtered) - 1, 1)
    dispersion = 0.0 if df == 0 else observed_chisq / df
    null_stats = bootstrap_homogeneous_null(
        depths,
        overall_probability,
        observed_chisq,
        bootstraps,
        seed,
        quantile_grid=quantile_grid,
    )

    observed_quantiles = (
        None
        if quantile_grid is None
        else (np.zeros(len(quantile_grid)) if fractions.size == 0 else np.quantile(fractions, quantile_grid))
    )
    zero_branch_introns = int(np.sum(counts == 0))

    summary_row = {
        "dataset_type": dataset_type,
        "dataset": dataset_name,
        "condition": condition,
        "signal_field": signal_field,
        "min_anchored_reads": min_reads,
        "introns_tested": len(filtered),
        "anchored_fragments": int(depths.sum()),
        "branched_fragments": int(counts.sum()),
        "overall_branch_fraction": overall_probability,
        "overall_branch_percent": overall_probability * 100.0,
        "mean_intron_branch_fraction": 0.0 if fractions.size == 0 else float(np.mean(fractions)),
        "median_intron_branch_fraction": 0.0 if fractions.size == 0 else float(np.median(fractions)),
        "zero_branch_introns_observed": zero_branch_introns,
        "zero_branch_introns_fraction_observed": 0.0
        if len(filtered) == 0
        else zero_branch_introns / len(filtered),
        "null_zero_branch_introns_mean": null_stats["null_zero_branch_introns_mean"],
        "null_zero_branch_introns_q025": null_stats["null_zero_branch_introns_q025"],
        "null_zero_branch_introns_q975": null_stats["null_zero_branch_introns_q975"],
        "pearson_chisq": observed_chisq,
        "pearson_df": df,
        "dispersion_ratio": dispersion,
        "bootstrap_pvalue_homogeneous": null_stats["bootstrap_pvalue"],
    }

    return summary_row, filtered, {
        "quantile_grid": quantile_grid,
        "observed_quantiles": observed_quantiles,
        "null_quantiles_lower": null_stats["null_quantiles_lower"],
        "null_quantiles_upper": null_stats["null_quantiles_upper"],
        "summary": summary_row,
    }


def rankdata(values):
    order = np.argsort(values, kind="mergesort")
    ranks = np.empty(values.size, dtype=np.float64)
    sorted_values = values[order]

    start = 0
    while start < values.size:
        end = start + 1
        while end < values.size and sorted_values[end] == sorted_values[start]:
            end += 1
        average_rank = (start + end - 1) / 2.0 + 1.0
        ranks[order[start:end]] = average_rank
        start = end
    return ranks


def correlation(x, y):
    if x.size < 2:
        return 0.0
    x_centered = x - np.mean(x)
    y_centered = y - np.mean(y)
    denominator = math.sqrt(float(np.sum(x_centered * x_centered) * np.sum(y_centered * y_centered)))
    if denominator == 0.0:
        return 0.0
    return float(np.sum(x_centered * y_centered) / denominator)


def build_condition_introns_rows(condition_records, min_reads):
    rows = []
    for condition in sorted(condition_records):
        records = filtered_records(condition_records[condition], min_reads)
        for rank, record in enumerate(
            sorted(records, key=lambda item: (-item["branch_fraction"], -item["anchored_fragments"], item["gene_name"])),
            start=1,
        ):
            rows.append(
                {
                    "condition": condition,
                    "rank_within_condition": rank,
                    **record,
                }
            )
    return rows


def compare_conditions(condition_records, query_condition, control_condition, min_reads):
    query_records = {
        record["intron_id"]: record for record in filtered_records(condition_records.get(query_condition, {}), min_reads)
    }
    control_records = {
        record["intron_id"]: record for record in filtered_records(condition_records.get(control_condition, {}), min_reads)
    }
    shared_introns = sorted(set(query_records) & set(control_records))

    comparison_rows = []
    for intron_id in shared_introns:
        query = query_records[intron_id]
        control = control_records[intron_id]
        comparison_rows.append(
            {
                "intron_id": intron_id,
                "gene_id": query["gene_id"],
                "gene_name": query["gene_name"],
                "transcript_id": query["transcript_id"],
                "intron_number": query["intron_number"],
                "chrom": query["chrom"],
                "strand": query["strand"],
                "intron_start": query["intron_start"],
                "intron_end": query["intron_end"],
                "branchpoint_position": query["branchpoint_position"],
                "branchpoint_to_3ss_nt": query["branchpoint_to_3ss_nt"],
                f"{query_condition}_anchored_fragments": query["anchored_fragments"],
                f"{query_condition}_branched_fragments": query["branched_fragments"],
                f"{query_condition}_branch_fraction": query["branch_fraction"],
                f"{control_condition}_anchored_fragments": control["anchored_fragments"],
                f"{control_condition}_branched_fragments": control["branched_fragments"],
                f"{control_condition}_branch_fraction": control["branch_fraction"],
                f"{query_condition}_minus_{control_condition}_branch_fraction": query["branch_fraction"]
                - control["branch_fraction"],
            }
        )

    comparison_rows.sort(
        key=lambda row: (
            -row[f"{query_condition}_branch_fraction"],
            -row[f"{query_condition}_minus_{control_condition}_branch_fraction"],
            row["gene_name"],
            row["intron_id"],
        )
    )

    query_fractions = np.array([row[f"{query_condition}_branch_fraction"] for row in comparison_rows], dtype=np.float64)
    control_fractions = np.array([row[f"{control_condition}_branch_fraction"] for row in comparison_rows], dtype=np.float64)
    pearson_r = correlation(query_fractions, control_fractions)
    spearman_r = correlation(rankdata(query_fractions), rankdata(control_fractions))

    decile_rows = []
    if comparison_rows:
        order = np.argsort(control_fractions)
        decile_bins = np.array_split(order, 10)
        for decile, indices in enumerate(decile_bins, start=1):
            if indices.size == 0:
                continue
            query_anchored = int(sum(comparison_rows[index][f"{query_condition}_anchored_fragments"] for index in indices))
            query_branched = int(sum(comparison_rows[index][f"{query_condition}_branched_fragments"] for index in indices))
            control_anchored = int(sum(comparison_rows[index][f"{control_condition}_anchored_fragments"] for index in indices))
            control_branched = int(sum(comparison_rows[index][f"{control_condition}_branched_fragments"] for index in indices))
            decile_rows.append(
                {
                    "control_decile": decile,
                    "intron_count": int(indices.size),
                    f"{query_condition}_anchored_fragments": query_anchored,
                    f"{query_condition}_branched_fragments": query_branched,
                    f"{query_condition}_branch_fraction": 0.0 if query_anchored == 0 else query_branched / query_anchored,
                    f"{control_condition}_anchored_fragments": control_anchored,
                    f"{control_condition}_branched_fragments": control_branched,
                    f"{control_condition}_branch_fraction": 0.0
                    if control_anchored == 0
                    else control_branched / control_anchored,
                }
            )

    query_top = decile_rows[-1][f"{query_condition}_branch_fraction"] if decile_rows else 0.0
    query_bottom = decile_rows[0][f"{query_condition}_branch_fraction"] if decile_rows else 0.0
    summary_row = {
        "query_condition": query_condition,
        "control_condition": control_condition,
        "shared_introns_compared": len(comparison_rows),
        "pearson_r": pearson_r,
        "spearman_r": spearman_r,
        f"{query_condition}_overall_branch_fraction_shared_introns": 0.0
        if not comparison_rows
        else float(
            sum(row[f"{query_condition}_branched_fragments"] for row in comparison_rows)
            / sum(row[f"{query_condition}_anchored_fragments"] for row in comparison_rows)
        ),
        f"{control_condition}_overall_branch_fraction_shared_introns": 0.0
        if not comparison_rows
        else float(
            sum(row[f"{control_condition}_branched_fragments"] for row in comparison_rows)
            / sum(row[f"{control_condition}_anchored_fragments"] for row in comparison_rows)
        ),
        f"{query_condition}_branch_fraction_top_{control_condition}_decile": query_top,
        f"{query_condition}_branch_fraction_bottom_{control_condition}_decile": query_bottom,
        f"{query_condition}_top_over_bottom_{control_condition}_decile_fold": 0.0
        if query_bottom == 0.0
        else query_top / query_bottom,
    }

    return summary_row, comparison_rows, decile_rows


def plot_results(
    condition_plot_data,
    condition_order,
    comparison_summary,
    decile_rows,
    query_condition,
    control_condition,
    signal_field,
    output_png,
    output_pdf,
):
    figure, axes = plt.subplots(1, 3, figsize=(15, 4.6), constrained_layout=True)

    pooled_conditions = [condition for condition in condition_order if condition in condition_plot_data]
    for axis, condition in zip(axes[:2], pooled_conditions[:2]):
        plot_data = condition_plot_data[condition]
        color = CONDITION_COLORS.get(condition, "#4c4c4c")
        percentiles = plot_data["quantile_grid"] * 100.0

        axis.fill_between(
            percentiles,
            plot_data["null_quantiles_lower"] * 100.0,
            plot_data["null_quantiles_upper"] * 100.0,
            color=color,
            alpha=0.16,
            label="Homogeneous null (95%)",
        )
        axis.plot(
            percentiles,
            plot_data["observed_quantiles"] * 100.0,
            color=color,
            linewidth=2.4,
            label="Observed",
        )
        summary = plot_data["summary"]
        axis.set_title(
            f"{condition} pooled\n{summary['introns_tested']} introns, "
            f"dispersion {summary['dispersion_ratio']:.2f}, "
            f"p={summary['bootstrap_pvalue_homogeneous']:.3g}"
        )
        axis.set_xlabel("Intron percentile")
        axis.set_ylabel("Branched fraction (%)")
        axis.set_xlim(0, 100)
        axis.legend(frameon=False)

    while len(pooled_conditions) < 2:
        axes[len(pooled_conditions)].axis("off")
        pooled_conditions.append(None)

    decile_axis = axes[2]
    deciles = [row["control_decile"] for row in decile_rows]
    query_values = [row[f"{query_condition}_branch_fraction"] * 100.0 for row in decile_rows]
    control_values = [row[f"{control_condition}_branch_fraction"] * 100.0 for row in decile_rows]
    decile_axis.plot(deciles, query_values, color=CONDITION_COLORS.get(query_condition, "#4c4c4c"), marker="o", linewidth=2.2, label=query_condition)
    decile_axis.plot(deciles, control_values, color=CONDITION_COLORS.get(control_condition, "#4c4c4c"), marker="o", linewidth=2.2, label=control_condition)
    decile_axis.set_xlabel(f"{control_condition} branching decile")
    decile_axis.set_ylabel("Branched fraction (%)")
    decile_axis.set_title(
        f"{query_condition} vs {control_condition} shared introns\n"
        f"Spearman {comparison_summary['spearman_r']:.2f} across "
        f"{comparison_summary['shared_introns_compared']} introns"
    )
    decile_axis.set_xticks(deciles)
    decile_axis.legend(frameon=False)

    for axis in axes:
        axis.spines["top"].set_visible(False)
        axis.spines["right"].set_visible(False)

    figure.suptitle(f"Intron-level branching heterogeneity using {signal_label(signal_field)}", fontsize=14)
    Path(output_png).parent.mkdir(parents=True, exist_ok=True)
    figure.savefig(output_png, dpi=300)
    figure.savefig(output_pdf)
    plt.close(figure)


def main():
    args = parse_args()
    sample_datasets, condition_datasets = load_datasets(args.site_counts, args.signal_field)

    sample_summary_rows = []
    condition_summary_rows = []
    condition_plot_data = {}

    sample_order = sorted(sample_datasets, key=lambda sample: (sample_datasets[sample]["condition"], sample))
    condition_order = sorted(condition_datasets)

    for index, sample in enumerate(sample_order):
        summary_row, _, _ = analyse_dataset(
            sample_datasets[sample]["records"],
            sample,
            "sample",
            sample_datasets[sample]["condition"],
            args.sample_min_reads,
            args.signal_field,
            args.bootstraps,
            args.seed + index,
        )
        sample_summary_rows.append(summary_row)

    for index, condition in enumerate(condition_order):
        summary_row, _, plot_data = analyse_dataset(
            condition_datasets[condition],
            condition,
            "condition",
            condition,
            args.condition_min_reads,
            args.signal_field,
            args.bootstraps,
            args.seed + 1000 + index,
            quantile_grid=PANEL_GRID,
        )
        condition_summary_rows.append(summary_row)
        condition_plot_data[condition] = plot_data

    condition_introns_rows = build_condition_introns_rows(condition_datasets, args.condition_min_reads)
    comparison_summary, comparison_rows, decile_rows = compare_conditions(
        condition_datasets,
        args.query_condition,
        args.control_condition,
        args.condition_min_reads,
    )

    if not decile_rows:
        raise ValueError(
            f"No shared introns passed the condition min-read threshold for {args.query_condition} vs {args.control_condition}"
        )
    if not condition_introns_rows:
        raise ValueError("No pooled-condition introns passed the configured minimum-read threshold")

    sample_summary_rows.sort(key=lambda row: (row["condition"], row["dataset"]))
    condition_summary_rows.sort(key=lambda row: row["condition"])

    write_rows(args.output_summary_by_sample, sample_summary_rows, list(sample_summary_rows[0].keys()))
    write_rows(args.output_summary_by_condition, condition_summary_rows, list(condition_summary_rows[0].keys()))
    write_rows(args.output_condition_introns, condition_introns_rows, list(condition_introns_rows[0].keys()))
    write_rows(
        args.output_condition_comparison_summary,
        [comparison_summary],
        list(comparison_summary.keys()),
    )
    write_rows(
        args.output_condition_comparison_introns,
        comparison_rows,
        list(comparison_rows[0].keys()),
    )
    write_rows(args.output_control_deciles, decile_rows, list(decile_rows[0].keys()))

    plot_results(
        condition_plot_data,
        condition_order,
        comparison_summary,
        decile_rows,
        args.query_condition,
        args.control_condition,
        args.signal_field,
        args.output_plot_png,
        args.output_plot_pdf,
    )

    print(f"Samples analysed: {len(sample_summary_rows)}")
    print(f"Conditions analysed: {len(condition_summary_rows)}")
    print(f"Shared introns in {args.query_condition} vs {args.control_condition}: {comparison_summary['shared_introns_compared']}")


if __name__ == "__main__":
    main()
