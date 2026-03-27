#!/usr/bin/env python3

import argparse
import csv
import gzip
import math
from collections import Counter, defaultdict
from pathlib import Path

import matplotlib
import numpy as np
import pysam

matplotlib.use("Agg")
import matplotlib.pyplot as plt


DNA_BASES = ["A", "C", "G", "T"]
RNA_BASE_LABELS = ["A", "C", "G", "U"]
CONDITION_COLORS = {
    "ILS": "#1f77b4",
    "DIS": "#d95f02",
}


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("--site-counts", action="append", dest="site_counts", required=True)
    parser.add_argument("--genome-fasta", required=True)
    parser.add_argument("--signal-field", default="zero_or_plus_one_branchpoint_fragments")
    parser.add_argument("--control-condition", default="ILS")
    parser.add_argument("--query-condition", default="DIS")
    parser.add_argument("--control-min-reads", type=int, default=20)
    parser.add_argument("--query-min-reads", type=int, default=20)
    parser.add_argument("--bp-upstream", type=int, default=5)
    parser.add_argument("--bp-downstream", type=int, default=5)
    parser.add_argument("--five-prime-intron-bases", type=int, default=6)
    parser.add_argument("--decile-fraction", type=float, default=0.1)
    parser.add_argument("--require-bp-center-base")
    parser.add_argument("--output-summary", required=True)
    parser.add_argument("--output-intron-contexts", required=True)
    parser.add_argument("--output-context-effects", required=True)
    parser.add_argument("--output-top-bottom-enrichment", required=True)
    parser.add_argument("--output-model-summary", required=True)
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


def reverse_complement(sequence):
    return sequence.translate(str.maketrans("ACGTNacgtn", "TGCANtgcan"))[::-1]


def signal_label(signal_field):
    labels = {
        "exact_branchpoint_fragments": "exact branchpoint-terminated fragments",
        "plus_one_branchpoint_fragments": "+1 branchpoint-proximal fragments",
        "zero_or_plus_one_branchpoint_fragments": "0/+1 branchpoint-proximal fragments",
    }
    return labels.get(signal_field, signal_field.replace("_", " "))


def load_condition_rows(paths, signal_field):
    condition_rows = defaultdict(dict)
    for path in paths:
        rows = read_tsv_rows(path)
        if not rows:
            continue
        condition = rows[0]["condition"]
        for row in rows:
            intron_id = row["intron_id"]
            pooled = condition_rows[condition].setdefault(
                intron_id,
                {
                    "condition": condition,
                    "intron_id": intron_id,
                    "gene_id": row["gene_id"],
                    "gene_name": row["gene_name"],
                    "transcript_id": row["transcript_id"],
                    "intron_number": row["intron_number"],
                    "chrom": row["chrom"],
                    "strand": row["strand"],
                    "intron_start": int(row["intron_start"]),
                    "intron_end": int(row["intron_end"]),
                    "intron_length": int(row["intron_length"]),
                    "three_prime_ss": int(row["three_prime_ss"]),
                    "branchpoint_position": int(row["branchpoint_position"]),
                    "branchpoint_score": float(row["branchpoint_score"]),
                    "branchpoint_to_3ss_nt": int(row["branchpoint_to_3ss_nt"]),
                    "branchpoint_type": row.get("branchpoint_type", ""),
                    "branchpoint_candidates": int(row["branchpoint_candidates"]),
                    "anchored_fragments": 0,
                    "signal_fragments": 0,
                },
            )
            pooled["anchored_fragments"] += count_value(row, "anchored_fragments")
            pooled["signal_fragments"] += count_value(row, signal_field)

    return condition_rows


def fetch_oriented_sequence(fasta, chrom, start_1based, end_1based, strand):
    sequence = fasta.fetch(chrom, start_1based - 1, end_1based).upper()
    if strand == "-":
        return reverse_complement(sequence)
    return sequence


def extract_sequences(record, fasta, bp_upstream, bp_downstream, five_prime_intron_bases):
    branchpoint_context = fetch_oriented_sequence(
        fasta,
        record["chrom"],
        record["branchpoint_position"] - bp_upstream,
        record["branchpoint_position"] + bp_downstream,
        record["strand"],
    )

    if record["strand"] == "+":
        intron_5p_start = record["intron_start"]
        intron_5p_end = record["intron_start"] + five_prime_intron_bases - 1
    else:
        intron_5p_start = record["intron_end"] - five_prime_intron_bases + 1
        intron_5p_end = record["intron_end"]

    five_prime_intron_context = fetch_oriented_sequence(
        fasta,
        record["chrom"],
        intron_5p_start,
        intron_5p_end,
        record["strand"],
    )
    return branchpoint_context, five_prime_intron_context


def build_condition_context_rows(
    condition_records,
    fasta_path,
    bp_upstream,
    bp_downstream,
    five_prime_intron_bases,
):
    rows = []
    with pysam.FastaFile(fasta_path) as fasta:
        for condition in sorted(condition_records):
            for intron_id in sorted(condition_records[condition], key=lambda key: (condition_records[condition][key]["gene_name"], key)):
                record = dict(condition_records[condition][intron_id])
                branchpoint_context, five_prime_intron_context = extract_sequences(
                    record,
                    fasta,
                    bp_upstream,
                    bp_downstream,
                    five_prime_intron_bases,
                )
                branch_fraction = (
                    0.0
                    if record["anchored_fragments"] == 0
                    else record["signal_fragments"] / record["anchored_fragments"]
                )
                record["branch_fraction"] = branch_fraction
                record["estimated_readthrough_fraction_if_all_branched"] = 1.0 - branch_fraction
                record["branchpoint_context"] = branchpoint_context
                record["five_prime_intron_context"] = five_prime_intron_context
                rows.append(record)
    return rows


def filter_context_rows(rows, condition, min_reads):
    return [
        row
        for row in rows
        if row["condition"] == condition and row["anchored_fragments"] >= min_reads
    ]


def overall_fraction(rows, field):
    numerator = sum(row[field] for row in rows)
    denominator = sum(row["anchored_fragments"] for row in rows)
    return 0.0 if denominator == 0 else numerator / denominator


def filter_rows_by_center_base(rows, required_base):
    if not required_base:
        return rows
    required_base = required_base.upper()
    filtered = []
    for row in rows:
        context = row["branchpoint_context"]
        if context[len(context) // 2] == required_base:
            filtered.append(row)
    return filtered


def build_context_effect_rows(rows, dataset_label, region_name, sequence_field, positions, overall_branch_fraction):
    context_rows = []
    overall_readthrough = 1.0 - overall_branch_fraction
    for index, position in enumerate(positions):
        for base in DNA_BASES:
            matching = [row for row in rows if row[sequence_field][index] == base]
            anchored = sum(row["anchored_fragments"] for row in matching)
            signal = sum(row["signal_fragments"] for row in matching)
            if not matching or anchored == 0:
                continue
            branch_fraction = 0.0 if anchored == 0 else signal / anchored
            readthrough = 1.0 - branch_fraction
            context_rows.append(
                {
                    "dataset": dataset_label,
                    "region": region_name,
                    "position": position,
                    "base": base,
                    "intron_count": len(matching),
                    "anchored_fragments": anchored,
                    "signal_fragments": signal,
                    "branch_fraction": branch_fraction,
                    "estimated_readthrough_fraction_if_all_branched": readthrough,
                    "delta_readthrough_from_dataset_mean": readthrough - overall_readthrough,
                    "delta_branch_fraction_from_dataset_mean": branch_fraction - overall_branch_fraction,
                }
            )
    return context_rows


def sigmoid(values):
    positive = values >= 0
    result = np.empty_like(values, dtype=float)
    result[positive] = 1.0 / (1.0 + np.exp(-values[positive]))
    exp_values = np.exp(values[~positive])
    result[~positive] = exp_values / (1.0 + exp_values)
    return result


def logit(probability):
    probability = min(max(probability, 1e-9), 1.0 - 1e-9)
    return math.log(probability / (1.0 - probability))


def build_design_matrix(rows, region_specs):
    columns = [np.ones(len(rows), dtype=float)]
    metadata = [{"feature_group": "intercept", "position": "", "base": "", "baseline_base": ""}]

    for feature_group, sequence_field, positions in region_specs:
        for index, position in enumerate(positions):
            bases = [row[sequence_field][index] for row in rows]
            counts = Counter(bases)
            if len(counts) < 2:
                continue
            baseline_base = sorted(counts, key=lambda base: (-counts[base], base))[0]
            for base in sorted(counts):
                if base == baseline_base:
                    continue
                column = np.array([1.0 if value == base else 0.0 for value in bases], dtype=float)
                if column.sum() == 0:
                    continue
                columns.append(column)
                metadata.append(
                    {
                        "feature_group": feature_group,
                        "position": position,
                        "base": base,
                        "baseline_base": baseline_base,
                    }
                )

    return np.column_stack(columns), metadata


def binomial_deviance(successes, totals, probabilities):
    probabilities = np.clip(probabilities, 1e-9, 1.0 - 1e-9)
    fitted_successes = totals * probabilities
    fitted_failures = totals - fitted_successes
    failures = totals - successes

    deviance = np.zeros_like(probabilities, dtype=float)
    nonzero_successes = successes > 0
    nonzero_failures = failures > 0
    deviance[nonzero_successes] += successes[nonzero_successes] * np.log(
        successes[nonzero_successes] / fitted_successes[nonzero_successes]
    )
    deviance[nonzero_failures] += failures[nonzero_failures] * np.log(
        failures[nonzero_failures] / fitted_failures[nonzero_failures]
    )
    return float(2.0 * np.sum(deviance))


def fit_binomial_glm(design_matrix, successes, totals, ridge=1e-6, max_iter=100, tolerance=1e-8):
    beta = np.zeros(design_matrix.shape[1], dtype=float)
    initial_probability = (successes.sum() + 0.5) / (totals.sum() + 1.0)
    beta[0] = logit(initial_probability)

    penalty = np.zeros(design_matrix.shape[1], dtype=float)
    penalty[1:] = ridge

    for _ in range(max_iter):
        eta = design_matrix @ beta
        probabilities = np.clip(sigmoid(eta), 1e-9, 1.0 - 1e-9)
        weights = totals * probabilities * (1.0 - probabilities)
        weights = np.maximum(weights, 1e-9)
        working_response = eta + (successes - totals * probabilities) / weights

        weighted_design = design_matrix * np.sqrt(weights)[:, None]
        weighted_response = working_response * np.sqrt(weights)
        lhs = weighted_design.T @ weighted_design
        lhs.flat[:: lhs.shape[0] + 1] += penalty
        rhs = weighted_design.T @ weighted_response
        try:
            new_beta = np.linalg.solve(lhs, rhs)
        except np.linalg.LinAlgError:
            new_beta = np.linalg.lstsq(lhs, rhs, rcond=None)[0]

        if np.max(np.abs(new_beta - beta)) < tolerance:
            beta = new_beta
            break
        beta = new_beta

    probabilities = np.clip(sigmoid(design_matrix @ beta), 1e-9, 1.0 - 1e-9)
    return beta, probabilities


def weighted_r_squared(observed_fraction, predicted_fraction, weights):
    weighted_mean = float(np.average(observed_fraction, weights=weights))
    residual_sum = float(np.sum(weights * (observed_fraction - predicted_fraction) ** 2))
    total_sum = float(np.sum(weights * (observed_fraction - weighted_mean) ** 2))
    if total_sum == 0.0:
        return 0.0
    return 1.0 - (residual_sum / total_sum)


def fit_sequence_models(rows, bp_positions, intron5_positions):
    successes = np.array([row["signal_fragments"] for row in rows], dtype=float)
    totals = np.array([row["anchored_fragments"] for row in rows], dtype=float)
    observed_fraction = successes / totals
    null_probability = float(successes.sum() / totals.sum())
    null_probabilities = np.full(len(rows), null_probability, dtype=float)
    null_deviance = binomial_deviance(successes, totals, null_probabilities)

    model_specs = [
        ("bp_flank_only", [("branchpoint_flank", "branchpoint_context", bp_positions)]),
        ("five_prime_intron_only", [("five_prime_intron", "five_prime_intron_context", intron5_positions)]),
        (
            "combined_additive",
            [
                ("branchpoint_flank", "branchpoint_context", bp_positions),
                ("five_prime_intron", "five_prime_intron_context", intron5_positions),
            ],
        ),
    ]

    summary_rows = []
    coefficient_rows = []
    for model_name, region_specs in model_specs:
        design_matrix, metadata = build_design_matrix(rows, region_specs)
        beta, probabilities = fit_binomial_glm(design_matrix, successes, totals)
        residual_deviance = binomial_deviance(successes, totals, probabilities)
        deviance_explained = 0.0 if null_deviance == 0.0 else 1.0 - (residual_deviance / null_deviance)
        weighted_r2 = weighted_r_squared(observed_fraction, probabilities, totals)

        summary_rows.append(
            {
                "model": model_name,
                "introns": len(rows),
                "anchored_fragments": int(totals.sum()),
                "signal_fragments": int(successes.sum()),
                "feature_columns": design_matrix.shape[1] - 1,
                "null_deviance": null_deviance,
                "residual_deviance": residual_deviance,
                "deviance_explained": deviance_explained,
                "weighted_r_squared": weighted_r2,
            }
        )

        for coefficient, meta in zip(beta, metadata):
            if meta["feature_group"] == "intercept":
                coefficient_rows.append(
                    {
                        "model": model_name,
                        "feature_group": "intercept",
                        "position": "",
                        "base": "",
                        "baseline_base": "",
                        "coefficient_log_odds": coefficient,
                        "odds_ratio": math.exp(coefficient),
                    }
                )
            else:
                coefficient_rows.append(
                    {
                        "model": model_name,
                        "feature_group": meta["feature_group"],
                        "position": meta["position"],
                        "base": meta["base"],
                        "baseline_base": meta["baseline_base"],
                        "coefficient_log_odds": coefficient,
                        "odds_ratio": math.exp(coefficient),
                    }
                )

    return summary_rows, coefficient_rows


def build_top_bottom_enrichment_rows(rows, region_name, sequence_field, positions, decile_fraction):
    ordered = sorted(
        rows,
        key=lambda row: (
            -row["estimated_readthrough_fraction_if_all_branched"],
            -row["anchored_fragments"],
            row["gene_name"],
            row["intron_id"],
        ),
    )
    group_size = max(1, int(round(len(ordered) * decile_fraction)))
    top_rows = ordered[:group_size]
    bottom_rows = ordered[-group_size:]

    enrichment_rows = []
    for index, position in enumerate(positions):
        for base in DNA_BASES:
            top_count = sum(1 for row in top_rows if row[sequence_field][index] == base)
            bottom_count = sum(1 for row in bottom_rows if row[sequence_field][index] == base)
            top_fraction = top_count / len(top_rows)
            bottom_fraction = bottom_count / len(bottom_rows)

            top_other = len(top_rows) - top_count
            bottom_other = len(bottom_rows) - bottom_count
            odds_ratio = ((top_count + 0.5) / (top_other + 0.5)) / ((bottom_count + 0.5) / (bottom_other + 0.5))
            enrichment_rows.append(
                {
                    "region": region_name,
                    "position": position,
                    "base": base,
                    "top_group_introns": len(top_rows),
                    "bottom_group_introns": len(bottom_rows),
                    "top_group_base_count": top_count,
                    "bottom_group_base_count": bottom_count,
                    "top_group_base_fraction": top_fraction,
                    "bottom_group_base_fraction": bottom_fraction,
                    "log2_odds_ratio_top_vs_bottom": math.log2(odds_ratio),
                }
            )
    return enrichment_rows, len(top_rows), len(bottom_rows)


def rows_to_matrix(rows, value_field, positions):
    position_to_index = {position: index for index, position in enumerate(positions)}
    base_to_index = {base: index for index, base in enumerate(DNA_BASES)}
    matrix = np.full((len(DNA_BASES), len(positions)), np.nan, dtype=float)
    for row in rows:
        matrix[base_to_index[row["base"]], position_to_index[row["position"]]] = row[value_field]
    return matrix


def draw_heatmap(axis, matrix, positions, title, colorbar_label, cmap, value_limit):
    image = axis.imshow(
        matrix,
        aspect="auto",
        interpolation="nearest",
        cmap=cmap,
        vmin=-value_limit,
        vmax=value_limit,
    )
    axis.set_xticks(range(len(positions)))
    axis.set_xticklabels(positions)
    axis.set_yticks(range(len(RNA_BASE_LABELS)))
    axis.set_yticklabels(RNA_BASE_LABELS)
    axis.set_title(title)
    axis.set_xlabel("Position")
    return image, colorbar_label


def plot_results(
    control_condition,
    signal_field,
    control_rows,
    bp_effect_rows,
    intron5_effect_rows,
    bp_enrichment_rows,
    intron5_enrichment_rows,
    bp_positions,
    intron5_positions,
    top_count,
    bottom_count,
    center_base_filter,
    model_summary_rows,
    output_png,
    output_pdf,
):
    overall_branch = overall_fraction(control_rows, "signal_fragments")
    overall_readthrough = 1.0 - overall_branch
    figure, axes = plt.subplots(2, 3, figsize=(15.5, 8), constrained_layout=True)

    bp_effect_matrix = rows_to_matrix(bp_effect_rows, "delta_readthrough_from_dataset_mean", bp_positions) * 100.0
    intron5_effect_matrix = rows_to_matrix(intron5_effect_rows, "delta_readthrough_from_dataset_mean", intron5_positions) * 100.0
    bp_enrichment_matrix = rows_to_matrix(bp_enrichment_rows, "log2_odds_ratio_top_vs_bottom", bp_positions)
    intron5_enrichment_matrix = rows_to_matrix(intron5_enrichment_rows, "log2_odds_ratio_top_vs_bottom", intron5_positions)

    effect_limit = max(
        float(np.nanmax(np.abs(bp_effect_matrix))) if not np.isnan(bp_effect_matrix).all() else 0.0,
        float(np.nanmax(np.abs(intron5_effect_matrix))) if not np.isnan(intron5_effect_matrix).all() else 0.0,
        1.0,
    )
    enrichment_limit = max(
        float(np.nanmax(np.abs(bp_enrichment_matrix))) if not np.isnan(bp_enrichment_matrix).all() else 0.0,
        float(np.nanmax(np.abs(intron5_enrichment_matrix))) if not np.isnan(intron5_enrichment_matrix).all() else 0.0,
        0.5,
    )

    image, label = draw_heatmap(
        axes[0, 0],
        bp_effect_matrix,
        bp_positions,
        f"{control_condition}: BP-flanking context\nDelta in estimated readthrough (percentage points)",
        "Delta readthrough (percentage points)",
        "coolwarm",
        effect_limit,
    )
    figure.colorbar(image, ax=axes[0, 0], label=label, shrink=0.86)

    image, label = draw_heatmap(
        axes[0, 1],
        intron5_effect_matrix,
        intron5_positions,
        f"{control_condition}: intron 5' start context\nDelta in estimated readthrough (percentage points)",
        "Delta readthrough (percentage points)",
        "coolwarm",
        effect_limit,
    )
    figure.colorbar(image, ax=axes[0, 1], label=label, shrink=0.86)

    image, label = draw_heatmap(
        axes[1, 0],
        bp_enrichment_matrix,
        bp_positions,
        f"{control_condition}: top vs bottom readthrough introns\nBP-flanking context, log2 odds ratio",
        "Log2 odds ratio",
        "PiYG",
        enrichment_limit,
    )
    figure.colorbar(image, ax=axes[1, 0], label=label, shrink=0.86)

    image, label = draw_heatmap(
        axes[1, 1],
        intron5_enrichment_matrix,
        intron5_positions,
        f"{control_condition}: top vs bottom readthrough introns\nIntron 5' start context, log2 odds ratio",
        "Log2 odds ratio",
        "PiYG",
        enrichment_limit,
    )
    figure.colorbar(image, ax=axes[1, 1], label=label, shrink=0.86)

    model_axis = axes[0, 2]
    model_names = [row["model"] for row in model_summary_rows]
    deviance_values = [row["deviance_explained"] * 100.0 for row in model_summary_rows]
    weighted_r2_values = [row["weighted_r_squared"] * 100.0 for row in model_summary_rows]
    x = np.arange(len(model_names))
    width = 0.38
    model_axis.bar(x - width / 2, deviance_values, width=width, color="#4c4c4c", alpha=0.85, label="Deviance explained")
    model_axis.bar(x + width / 2, weighted_r2_values, width=width, color="#9a9a9a", alpha=0.85, label="Weighted R^2")
    model_axis.set_xticks(x)
    model_axis.set_xticklabels(["BP", "5' intron", "Combined"])
    model_axis.set_ylabel("Explained variation (%)")
    model_axis.set_title(f"{control_condition}: additive sequence model")
    model_axis.legend(frameon=False)

    text_axis = axes[1, 2]
    text_axis.axis("off")
    summary_lines = [
        f"Signal: {signal_label(signal_field)}",
        f"Control introns: {len(control_rows)}",
        f"Mean estimated readthrough: {overall_readthrough * 100.0:.1f}%",
        f"Top / bottom groups: {top_count} / {bottom_count}",
    ]
    if center_base_filter:
        summary_lines.append(f"Restricted to BP center = {center_base_filter}")
    best_model = max(model_summary_rows, key=lambda row: row["deviance_explained"])
    summary_lines.append(
        f"Best additive model: {best_model['model']}, "
        f"{best_model['deviance_explained'] * 100.0:.1f}% deviance explained"
    )
    text_axis.text(0.0, 1.0, "\n".join(summary_lines), va="top", ha="left")

    for axis in axes.flat:
        axis.spines["top"].set_visible(False)
        axis.spines["right"].set_visible(False)

    figure.suptitle(
        f"Sequence context around branchpoint-associated RT arrest using {signal_label(signal_field)}\n"
        f"{control_condition} pooled, {len(control_rows)} introns, estimated mean readthrough {overall_readthrough * 100.0:.1f}%, "
        f"top/bottom groups {top_count} / {bottom_count} introns",
        fontsize=13,
    )
    Path(output_png).parent.mkdir(parents=True, exist_ok=True)
    figure.savefig(output_png, dpi=300)
    figure.savefig(output_pdf)
    plt.close(figure)


def main():
    args = parse_args()
    if not (0.0 < args.decile_fraction <= 0.5):
        raise ValueError("--decile-fraction must be in (0, 0.5]")

    condition_records = load_condition_rows(args.site_counts, args.signal_field)
    context_rows = build_condition_context_rows(
        condition_records,
        args.genome_fasta,
        args.bp_upstream,
        args.bp_downstream,
        args.five_prime_intron_bases,
    )
    context_rows = filter_rows_by_center_base(context_rows, args.require_bp_center_base)

    control_rows = filter_context_rows(context_rows, args.control_condition, args.control_min_reads)
    query_rows = filter_context_rows(context_rows, args.query_condition, args.query_min_reads)
    if not control_rows:
        raise ValueError(f"No introns passed control filter for {args.control_condition}")

    control_by_intron = {row["intron_id"]: row for row in control_rows}
    query_by_intron = {row["intron_id"]: row for row in query_rows}
    shared_introns = sorted(set(control_by_intron) & set(query_by_intron))
    shared_control_rows = [control_by_intron[intron_id] for intron_id in shared_introns]
    shared_query_rows = [query_by_intron[intron_id] for intron_id in shared_introns]

    bp_positions = list(range(-args.bp_upstream, args.bp_downstream + 1))
    intron5_positions = list(range(1, args.five_prime_intron_bases + 1))

    control_branch_fraction = overall_fraction(control_rows, "signal_fragments")
    query_branch_fraction = overall_fraction(query_rows, "signal_fragments")
    shared_control_branch_fraction = overall_fraction(shared_control_rows, "signal_fragments")
    shared_query_branch_fraction = overall_fraction(shared_query_rows, "signal_fragments")

    effect_rows = []
    effect_rows.extend(
        build_context_effect_rows(
            control_rows,
            f"{args.control_condition}_pooled",
            "branchpoint_flank",
            "branchpoint_context",
            bp_positions,
            control_branch_fraction,
        )
    )
    effect_rows.extend(
        build_context_effect_rows(
            control_rows,
            f"{args.control_condition}_pooled",
            "five_prime_intron",
            "five_prime_intron_context",
            intron5_positions,
            control_branch_fraction,
        )
    )

    if query_rows:
        effect_rows.extend(
            build_context_effect_rows(
                query_rows,
                f"{args.query_condition}_pooled",
                "branchpoint_flank",
                "branchpoint_context",
                bp_positions,
                query_branch_fraction,
            )
        )
        effect_rows.extend(
            build_context_effect_rows(
                query_rows,
                f"{args.query_condition}_pooled",
                "five_prime_intron",
                "five_prime_intron_context",
                intron5_positions,
                query_branch_fraction,
            )
        )

    if shared_control_rows and shared_query_rows:
        effect_rows.extend(
            build_context_effect_rows(
                shared_control_rows,
                f"{args.control_condition}_shared_with_{args.query_condition}",
                "branchpoint_flank",
                "branchpoint_context",
                bp_positions,
                shared_control_branch_fraction,
            )
        )
        effect_rows.extend(
            build_context_effect_rows(
                shared_control_rows,
                f"{args.control_condition}_shared_with_{args.query_condition}",
                "five_prime_intron",
                "five_prime_intron_context",
                intron5_positions,
                shared_control_branch_fraction,
            )
        )
        effect_rows.extend(
            build_context_effect_rows(
                shared_query_rows,
                f"{args.query_condition}_shared_with_{args.control_condition}",
                "branchpoint_flank",
                "branchpoint_context",
                bp_positions,
                shared_query_branch_fraction,
            )
        )
        effect_rows.extend(
            build_context_effect_rows(
                shared_query_rows,
                f"{args.query_condition}_shared_with_{args.control_condition}",
                "five_prime_intron",
                "five_prime_intron_context",
                intron5_positions,
                shared_query_branch_fraction,
            )
        )

    bp_enrichment_rows, top_count, bottom_count = build_top_bottom_enrichment_rows(
        control_rows,
        "branchpoint_flank",
        "branchpoint_context",
        bp_positions,
        args.decile_fraction,
    )
    intron5_enrichment_rows, _, _ = build_top_bottom_enrichment_rows(
        control_rows,
        "five_prime_intron",
        "five_prime_intron_context",
        intron5_positions,
        args.decile_fraction,
    )
    enrichment_rows = bp_enrichment_rows + intron5_enrichment_rows
    model_summary_rows, model_coefficient_rows = fit_sequence_models(control_rows, bp_positions, intron5_positions)

    summary_rows = [
        {
            "dataset": f"{args.control_condition}_pooled",
            "condition": args.control_condition,
            "min_anchored_reads": args.control_min_reads,
            "require_bp_center_base": args.require_bp_center_base or "",
            "introns": len(control_rows),
            "anchored_fragments": sum(row["anchored_fragments"] for row in control_rows),
            "signal_fragments": sum(row["signal_fragments"] for row in control_rows),
            "branch_fraction": control_branch_fraction,
            "estimated_readthrough_fraction_if_all_branched": 1.0 - control_branch_fraction,
        }
    ]
    if query_rows:
        summary_rows.append(
            {
                "dataset": f"{args.query_condition}_pooled",
                "condition": args.query_condition,
                "min_anchored_reads": args.query_min_reads,
                "require_bp_center_base": args.require_bp_center_base or "",
                "introns": len(query_rows),
                "anchored_fragments": sum(row["anchored_fragments"] for row in query_rows),
                "signal_fragments": sum(row["signal_fragments"] for row in query_rows),
                "branch_fraction": query_branch_fraction,
                "estimated_readthrough_fraction_if_all_branched": 1.0 - query_branch_fraction,
            }
        )
    if shared_control_rows and shared_query_rows:
        summary_rows.extend(
            [
                {
                    "dataset": f"{args.control_condition}_shared_with_{args.query_condition}",
                    "condition": args.control_condition,
                    "min_anchored_reads": args.control_min_reads,
                    "require_bp_center_base": args.require_bp_center_base or "",
                    "introns": len(shared_control_rows),
                    "anchored_fragments": sum(row["anchored_fragments"] for row in shared_control_rows),
                    "signal_fragments": sum(row["signal_fragments"] for row in shared_control_rows),
                    "branch_fraction": shared_control_branch_fraction,
                    "estimated_readthrough_fraction_if_all_branched": 1.0 - shared_control_branch_fraction,
                },
                {
                    "dataset": f"{args.query_condition}_shared_with_{args.control_condition}",
                    "condition": args.query_condition,
                    "min_anchored_reads": args.query_min_reads,
                    "require_bp_center_base": args.require_bp_center_base or "",
                    "introns": len(shared_query_rows),
                    "anchored_fragments": sum(row["anchored_fragments"] for row in shared_query_rows),
                    "signal_fragments": sum(row["signal_fragments"] for row in shared_query_rows),
                    "branch_fraction": shared_query_branch_fraction,
                    "estimated_readthrough_fraction_if_all_branched": 1.0 - shared_query_branch_fraction,
                },
            ]
        )

    effect_rows.sort(key=lambda row: (row["dataset"], row["region"], row["position"], row["base"]))
    enrichment_rows.sort(key=lambda row: (row["region"], row["position"], row["base"]))
    context_rows.sort(key=lambda row: (row["condition"], row["gene_name"], row["intron_id"]))

    write_rows(args.output_summary, summary_rows, list(summary_rows[0].keys()))
    write_rows(args.output_intron_contexts, context_rows, list(context_rows[0].keys()))
    write_rows(args.output_context_effects, effect_rows, list(effect_rows[0].keys()))
    write_rows(args.output_top_bottom_enrichment, enrichment_rows, list(enrichment_rows[0].keys()))
    write_rows(args.output_model_summary, model_summary_rows, list(model_summary_rows[0].keys()))

    control_bp_effect_rows = [
        row for row in effect_rows if row["dataset"] == f"{args.control_condition}_pooled" and row["region"] == "branchpoint_flank"
    ]
    control_intron5_effect_rows = [
        row for row in effect_rows if row["dataset"] == f"{args.control_condition}_pooled" and row["region"] == "five_prime_intron"
    ]
    plot_results(
        args.control_condition,
        args.signal_field,
        control_rows,
        control_bp_effect_rows,
        control_intron5_effect_rows,
        bp_enrichment_rows,
        intron5_enrichment_rows,
        bp_positions,
        intron5_positions,
        top_count,
        bottom_count,
        args.require_bp_center_base or "",
        model_summary_rows,
        args.output_plot_png,
        args.output_plot_pdf,
    )

    print(f"Control introns analysed: {len(control_rows)}")
    print(f"Query introns analysed: {len(query_rows)}")
    print(f"Shared introns analysed: {len(shared_introns)}")


if __name__ == "__main__":
    main()
