#!/usr/bin/env python3

import argparse
import csv
import gzip
import math
from collections import Counter
from pathlib import Path

import matplotlib
import numpy as np
import pysam

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.font_manager import FontProperties
from matplotlib.patches import PathPatch
from matplotlib.textpath import TextPath
from matplotlib.transforms import Affine2D


DNA_BASES = ["A", "C", "G", "T"]
BASE_COLORS = {
    "A": "#2ca02c",
    "C": "#1f77b4",
    "G": "#ff7f0e",
    "T": "#d62728",
}
BASE_DISPLAY = {
    "A": "A",
    "C": "C",
    "G": "G",
    "T": "U",
}
FEATURE_LABELS = {
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
    "branchpoint_candidates": "BP candidate count",
}
COMPARISON_FEATURES = [
    "intron_gc_3p_window",
    "intron_gc",
    "intron_gc_5p_window",
    "branchpoint_score",
    "log_intron_length",
    "bp_to_3ss_pyrimidine_fraction",
    "acceptor_maxent",
    "donor_maxent",
    "branchpoint_to_3ss_nt",
    "three_prime_window_pyrimidine_fraction",
]
PINNED_PLOT_FEATURES = ["branchpoint_to_3ss_nt"]
GROUPS = [
    ("threshold_high", "High-branch DIS introns"),
    ("other", "Other shared introns"),
]
GLYPH_FONT = FontProperties(family="DejaVu Sans", weight="bold")
GLYPHS = {base: TextPath((0, 0), BASE_DISPLAY[base], size=1, prop=GLYPH_FONT) for base in DNA_BASES}


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("--shared-introns", required=True)
    parser.add_argument("--genome-fasta", required=True)
    parser.add_argument("--control-condition", default="ILS")
    parser.add_argument("--query-condition", default="DIS")
    parser.add_argument("--query-branch-threshold", type=float, default=0.25)
    parser.add_argument("--reference-residual-quantile-fraction", type=float, default=0.1)
    parser.add_argument("--bp-upstream", type=int, default=10)
    parser.add_argument("--bp-downstream", type=int, default=10)
    parser.add_argument("--five-prime-upstream", type=int, default=5)
    parser.add_argument("--five-prime-downstream", type=int, default=10)
    parser.add_argument("--three-prime-upstream", type=int, default=20)
    parser.add_argument("--three-prime-downstream", type=int, default=3)
    parser.add_argument("--output-summary", required=True)
    parser.add_argument("--output-intron-groups", required=True)
    parser.add_argument("--output-feature-group-comparison", required=True)
    parser.add_argument("--output-position-stats", required=True)
    parser.add_argument("--output-feature-distribution-plot-png", required=True)
    parser.add_argument("--output-feature-distribution-plot-pdf", required=True)
    parser.add_argument("--output-sequence-logo-plot-png", required=True)
    parser.add_argument("--output-sequence-logo-plot-pdf", required=True)
    parser.add_argument("--output-sequence-frequency-plot-png", required=True)
    parser.add_argument("--output-sequence-frequency-plot-pdf", required=True)
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


def reverse_complement(sequence):
    return sequence.translate(str.maketrans("ACGTNacgtn", "TGCANtgcan"))[::-1]


def fetch_raw_sequence(fasta, chrom, start_0based, end_0based):
    if start_0based < 0 or end_0based <= start_0based:
        return None
    try:
        reference_length = fasta.get_reference_length(chrom)
    except KeyError:
        return None
    if end_0based > reference_length:
        return None
    return fasta.fetch(chrom, start_0based, end_0based).upper()


def fetch_oriented_sequence(fasta, chrom, start_1based, end_1based, strand):
    sequence = fetch_raw_sequence(fasta, chrom, start_1based - 1, end_1based)
    if sequence is None:
        return None
    if strand == "-":
        return reverse_complement(sequence)
    return sequence


def valid_dna(sequence):
    return sequence is not None and all(base in DNA_BASES for base in sequence)


def classify_rows(rows, query_branch_threshold, residual_quantile_fraction):
    residual_sorted = sorted(
        rows,
        key=lambda row: (
            -float(row["baseline_residual_query_branch_fraction"]),
            -int(row["query_anchored_fragments"]),
            row["gene_name"],
            row["intron_id"],
        ),
    )
    residual_group_size = max(1, int(round(len(residual_sorted) * residual_quantile_fraction)))
    residual_high_ids = {row["intron_id"] for row in residual_sorted[:residual_group_size]}

    threshold_sorted = sorted(
        rows,
        key=lambda row: (
            -float(row["query_branch_fraction"]),
            -int(row["query_anchored_fragments"]),
            row["gene_name"],
            row["intron_id"],
        ),
    )
    classified = []
    for rank, row in enumerate(threshold_sorted, start=1):
        query_branch_fraction = float(row["query_branch_fraction"])
        threshold_high = query_branch_fraction > query_branch_threshold
        classified.append(
            {
                **row,
                "query_branch_rank": rank,
                "threshold_group": "threshold_high" if threshold_high else "other",
                "meets_query_branch_threshold": "1" if threshold_high else "0",
                "in_reference_residual_high": "1" if row["intron_id"] in residual_high_ids else "0",
            }
        )
    return classified, residual_high_ids


def extract_branchpoint_sequence(row, fasta, upstream, downstream):
    sequence = fetch_oriented_sequence(
        fasta,
        row["chrom"],
        int(row["branchpoint_position"]) - upstream,
        int(row["branchpoint_position"]) + downstream,
        row["strand"],
    )
    expected_length = upstream + downstream + 1
    if sequence is None or len(sequence) != expected_length or not valid_dna(sequence):
        return None
    return sequence


def extract_five_prime_sequence(row, fasta, upstream, downstream):
    intron_start = int(row["intron_start"])
    intron_end = int(row["intron_end"])
    if row["strand"] == "+":
        start_1based = intron_start - upstream
        end_1based = intron_start + downstream
    else:
        start_1based = intron_end - downstream
        end_1based = intron_end + upstream
    sequence = fetch_oriented_sequence(fasta, row["chrom"], start_1based, end_1based, row["strand"])
    expected_length = upstream + downstream + 1
    if sequence is None or len(sequence) != expected_length or not valid_dna(sequence):
        return None
    return sequence


def extract_three_prime_sequence(row, fasta, upstream, downstream):
    intron_start = int(row["intron_start"])
    intron_end = int(row["intron_end"])
    if row["strand"] == "+":
        start_1based = intron_end - upstream + 1
        end_1based = intron_end + downstream
    else:
        start_1based = intron_start - downstream
        end_1based = intron_start + upstream - 1
    sequence = fetch_oriented_sequence(fasta, row["chrom"], start_1based, end_1based, row["strand"])
    expected_length = upstream + downstream
    if sequence is None or len(sequence) != expected_length or not valid_dna(sequence):
        return None
    return sequence


def build_output_rows(
    classified_rows,
    fasta_path,
    bp_upstream,
    bp_downstream,
    five_prime_upstream,
    five_prime_downstream,
    three_prime_upstream,
    three_prime_downstream,
):
    output_rows = []
    with pysam.FastaFile(fasta_path) as fasta:
        for row in classified_rows:
            bp_sequence = extract_branchpoint_sequence(row, fasta, bp_upstream, bp_downstream)
            five_prime_sequence = extract_five_prime_sequence(row, fasta, five_prime_upstream, five_prime_downstream)
            three_prime_sequence = extract_three_prime_sequence(row, fasta, three_prime_upstream, three_prime_downstream)

            exclusion_reasons = []
            if bp_sequence is None:
                exclusion_reasons.append("invalid_branchpoint_window")
            if five_prime_sequence is None:
                exclusion_reasons.append("invalid_five_prime_window")
            if three_prime_sequence is None:
                exclusion_reasons.append("invalid_three_prime_window")

            output_rows.append(
                {
                    "query_branch_rank": row["query_branch_rank"],
                    "threshold_group": row["threshold_group"],
                    "meets_query_branch_threshold": row["meets_query_branch_threshold"],
                    "in_reference_residual_high": row["in_reference_residual_high"],
                    "intron_id": row["intron_id"],
                    "gene_id": row["gene_id"],
                    "gene_name": row["gene_name"],
                    "transcript_id": row["transcript_id"],
                    "intron_number": row["intron_number"],
                    "chrom": row["chrom"],
                    "strand": row["strand"],
                    "intron_start": row["intron_start"],
                    "intron_end": row["intron_end"],
                    "intron_length": row["intron_length"],
                    "branchpoint_position": row["branchpoint_position"],
                    "branchpoint_to_3ss_nt": row["branchpoint_to_3ss_nt"],
                    "branchpoint_score": row["branchpoint_score"],
                    "branchpoint_candidates": row["branchpoint_candidates"],
                    "control_anchored_fragments": row["control_anchored_fragments"],
                    "control_branched_fragments": row["control_branched_fragments"],
                    "control_branch_fraction": row["control_branch_fraction"],
                    "query_anchored_fragments": row["query_anchored_fragments"],
                    "query_branched_fragments": row["query_branched_fragments"],
                    "query_branch_fraction": row["query_branch_fraction"],
                    "baseline_residual_query_branch_fraction": row["baseline_residual_query_branch_fraction"],
                    "donor_maxent": row["donor_maxent"],
                    "acceptor_maxent": row["acceptor_maxent"],
                    "intron_gc": row["intron_gc"],
                    "intron_gc_5p_window": row["intron_gc_5p_window"],
                    "intron_gc_3p_window": row["intron_gc_3p_window"],
                    "bp_to_3ss_pyrimidine_fraction": row["bp_to_3ss_pyrimidine_fraction"],
                    "three_prime_window_pyrimidine_fraction": row["three_prime_window_pyrimidine_fraction"],
                    "branchpoint_logo_sequence": bp_sequence or "",
                    "five_prime_ss_logo_sequence": five_prime_sequence or "",
                    "three_prime_ss_logo_sequence": three_prime_sequence or "",
                    "included_in_logo": "1" if not exclusion_reasons else "0",
                    "exclusion_reason": ",".join(exclusion_reasons),
                }
            )
    return output_rows


def quantile(values, q):
    return float(np.quantile(np.array(values, dtype=float), q))


def feature_specs():
    specs = []
    for model_feature in COMPARISON_FEATURES:
        plot_feature = "intron_length" if model_feature == "log_intron_length" else model_feature
        specs.append(
            {
                "model_feature": model_feature,
                "plot_feature": plot_feature,
                "label": FEATURE_LABELS.get(model_feature, model_feature.replace("_", " ")),
                "log_scale": plot_feature == "intron_length",
            }
        )
    return specs


def build_group_comparison_rows(rows, specs):
    high_rows = [row for row in rows if row["threshold_group"] == "threshold_high"]
    other_rows = [row for row in rows if row["threshold_group"] == "other"]
    comparison_rows = []
    for spec in specs:
        plot_feature = spec["plot_feature"]
        high_values = np.array([float(row[plot_feature]) for row in high_rows], dtype=float)
        other_values = np.array([float(row[plot_feature]) for row in other_rows], dtype=float)
        high_mean = float(np.mean(high_values))
        other_mean = float(np.mean(other_values))
        high_var = float(np.var(high_values, ddof=1)) if len(high_values) > 1 else 0.0
        other_var = float(np.var(other_values, ddof=1)) if len(other_values) > 1 else 0.0
        pooled_sd = math.sqrt(max((high_var + other_var) / 2.0, 1e-12))
        comparison_rows.append(
            {
                "model_feature": spec["model_feature"],
                "plot_feature": plot_feature,
                "label": spec["label"],
                "threshold_high_introns": len(high_rows),
                "other_introns": len(other_rows),
                "threshold_high_mean": high_mean,
                "other_mean": other_mean,
                "delta_mean": high_mean - other_mean,
                "threshold_high_median": float(np.median(high_values)),
                "other_median": float(np.median(other_values)),
                "delta_median": float(np.median(high_values) - np.median(other_values)),
                "threshold_high_q25": quantile(high_values, 0.25),
                "threshold_high_q75": quantile(high_values, 0.75),
                "other_q25": quantile(other_values, 0.25),
                "other_q75": quantile(other_values, 0.75),
                "standardized_mean_difference": (high_mean - other_mean) / pooled_sd,
                "abs_standardized_mean_difference": abs((high_mean - other_mean) / pooled_sd),
            }
        )
    comparison_rows.sort(key=lambda row: (-row["abs_standardized_mean_difference"], row["model_feature"]))
    return comparison_rows


def select_plot_rows(comparison_rows, max_features=8):
    selected = comparison_rows[:max_features]
    selected_by_feature = {row["model_feature"]: row for row in selected}
    all_by_feature = {row["model_feature"]: row for row in comparison_rows}

    for feature in PINNED_PLOT_FEATURES:
        if feature not in all_by_feature or feature in selected_by_feature:
            continue
        if not selected:
            selected.append(all_by_feature[feature])
            selected_by_feature[feature] = all_by_feature[feature]
            continue
        selected[-1] = all_by_feature[feature]
        selected_by_feature = {row["model_feature"]: row for row in selected}

    deduped = []
    seen = set()
    for row in selected:
        if row["model_feature"] in seen:
            continue
        deduped.append(row)
        seen.add(row["model_feature"])
    return deduped


def summarize_group_separation(rows):
    if not rows:
        return "", 0.0, 0.0
    best = rows[0]
    mean_abs = float(np.mean([row["abs_standardized_mean_difference"] for row in rows]))
    return best["model_feature"], best["abs_standardized_mean_difference"], mean_abs


def build_reference_residual_comparison(rows, residual_high_ids, specs):
    residual_rows = []
    for row in rows:
        residual_rows.append(
            {
                **row,
                "threshold_group": "threshold_high" if row["intron_id"] in residual_high_ids else "other",
            }
        )
    return build_group_comparison_rows(residual_rows, specs)


def build_position_stats(rows, sequence_field, positions, region_name):
    position_rows = []
    for group_name, _ in GROUPS:
        sequences = [
            row[sequence_field]
            for row in rows
            if row["threshold_group"] == group_name and row["included_in_logo"] == "1"
        ]
        if not sequences:
            continue
        for index, position in enumerate(positions):
            counts = Counter(sequence[index] for sequence in sequences)
            total = len(sequences)
            fractions = {base: counts.get(base, 0) / total for base in DNA_BASES}
            entropy = -sum(fraction * math.log2(fraction) for fraction in fractions.values() if fraction > 0.0)
            information = max(0.0, math.log2(len(DNA_BASES)) - entropy)
            for base in DNA_BASES:
                position_rows.append(
                    {
                        "region": region_name,
                        "threshold_group": group_name,
                        "position": position,
                        "base": base,
                        "display_base": BASE_DISPLAY[base],
                        "sequence_count": total,
                        "base_count": counts.get(base, 0),
                        "base_fraction": fractions[base],
                        "information_bits": information,
                        "logo_height_bits": fractions[base] * information,
                    }
                )
    return position_rows


def position_tick_labels(positions):
    labels = []
    for position in positions:
        if position in {positions[0], positions[-1], 0} or position % 5 == 0:
            labels.append(str(position))
        else:
            labels.append("")
    return labels


def draw_letter(axis, base, x_left, y_bottom, height_bits, width=0.9):
    if height_bits <= 0.0:
        return
    glyph = GLYPHS[base]
    bbox = glyph.get_extents()
    transform = (
        Affine2D()
        .translate(-bbox.xmin, -bbox.ymin)
        .scale(width / bbox.width, height_bits / bbox.height)
        .translate(x_left, y_bottom)
    )
    axis.add_patch(
        PathPatch(glyph, transform=transform + axis.transData, facecolor=BASE_COLORS[base], edgecolor="none")
    )


def plot_logo(axis, stats_rows, positions, title, xlabel, reference_mode):
    return plot_logo_generic(
        axis,
        stats_rows,
        positions,
        title,
        xlabel,
        reference_mode,
        height_field="logo_height_bits",
        y_label="Information (bits)",
        y_limit=2.05,
    )


def plot_logo_generic(axis, stats_rows, positions, title, xlabel, reference_mode, height_field, y_label, y_limit):
    sequence_count = max((int(row["sequence_count"]) for row in stats_rows), default=0)
    heights_by_position = {position: {base: 0.0 for base in DNA_BASES} for position in positions}
    for row in stats_rows:
        heights_by_position[int(row["position"])][row["base"]] = float(row[height_field])

    for index, position in enumerate(positions):
        y_bottom = 0.0
        ordered_bases = sorted(DNA_BASES, key=lambda base: heights_by_position[position][base])
        for base in ordered_bases:
            height = heights_by_position[position][base]
            draw_letter(axis, base, index + 0.05, y_bottom, height)
            y_bottom += height

    axis.set_xlim(0, len(positions))
    axis.set_ylim(0, y_limit)
    axis.set_xticks([index + 0.5 for index in range(len(positions))])
    axis.set_xticklabels(position_tick_labels(positions), fontsize=8)
    axis.set_ylabel(y_label)
    axis.set_xlabel(xlabel)
    axis.set_title(f"{title}\n(n = {sequence_count})", fontsize=11)
    reference_x = positions.index(0) + 0.5 if reference_mode == "center" else positions.index(0)
    axis.axvline(reference_x, color="black", linewidth=1.0, linestyle="--", alpha=0.45)
    axis.spines["top"].set_visible(False)
    axis.spines["right"].set_visible(False)


def plot_sequence_logos(
    position_rows,
    query_condition,
    query_branch_threshold,
    bp_positions,
    five_prime_positions,
    three_prime_positions,
    output_png,
    output_pdf,
):
    figure, axes = plt.subplots(3, 2, figsize=(15, 11), constrained_layout=True)
    for column, (group_name, group_label) in enumerate(GROUPS):
        bp_rows = [row for row in position_rows if row["region"] == "branchpoint_flank" and row["threshold_group"] == group_name]
        five_prime_rows = [row for row in position_rows if row["region"] == "five_prime_ss" and row["threshold_group"] == group_name]
        three_prime_rows = [row for row in position_rows if row["region"] == "three_prime_ss" and row["threshold_group"] == group_name]
        plot_logo(
            axes[0, column],
            bp_rows,
            bp_positions,
            f"{group_label}\nBranchpoint flank",
            "Position relative to branchpoint A (nt)",
            "center",
        )
        plot_logo(
            axes[1, column],
            five_prime_rows,
            five_prime_positions,
            f"{group_label}\n5' splice site window",
            "Position relative to first intronic base at the 5'SS (nt)",
            "between",
        )
        plot_logo(
            axes[2, column],
            three_prime_rows,
            three_prime_positions,
            f"{group_label}\n3' splice site window",
            "Position relative to first exonic base at the 3'SS (nt)",
            "between",
        )

    figure.suptitle(
        f"{query_condition} high-branch introns defined by absolute branch fraction\n"
        f"High group = pooled {query_condition} branch fraction > {query_branch_threshold:.2f}",
        fontsize=13,
    )
    Path(output_png).parent.mkdir(parents=True, exist_ok=True)
    figure.savefig(output_png, dpi=300)
    figure.savefig(output_pdf)
    plt.close(figure)


def plot_sequence_frequencies(
    position_rows,
    query_condition,
    query_branch_threshold,
    bp_positions,
    five_prime_positions,
    three_prime_positions,
    output_png,
    output_pdf,
):
    figure, axes = plt.subplots(3, 2, figsize=(15, 11), constrained_layout=True)
    for column, (group_name, group_label) in enumerate(GROUPS):
        bp_rows = [row for row in position_rows if row["region"] == "branchpoint_flank" and row["threshold_group"] == group_name]
        five_prime_rows = [row for row in position_rows if row["region"] == "five_prime_ss" and row["threshold_group"] == group_name]
        three_prime_rows = [row for row in position_rows if row["region"] == "three_prime_ss" and row["threshold_group"] == group_name]
        plot_logo_generic(
            axes[0, column],
            bp_rows,
            bp_positions,
            f"{group_label}\nBranchpoint flank",
            "Position relative to branchpoint A (nt)",
            "center",
            height_field="base_fraction",
            y_label="Base frequency",
            y_limit=1.02,
        )
        plot_logo_generic(
            axes[1, column],
            five_prime_rows,
            five_prime_positions,
            f"{group_label}\n5' splice site window",
            "Position relative to first intronic base at the 5'SS (nt)",
            "between",
            height_field="base_fraction",
            y_label="Base frequency",
            y_limit=1.02,
        )
        plot_logo_generic(
            axes[2, column],
            three_prime_rows,
            three_prime_positions,
            f"{group_label}\n3' splice site window",
            "Position relative to first exonic base at the 3'SS (nt)",
            "between",
            height_field="base_fraction",
            y_label="Base frequency",
            y_limit=1.02,
        )

    figure.suptitle(
        f"{query_condition} high-branch introns defined by absolute branch fraction\n"
        f"High group = pooled {query_condition} branch fraction > {query_branch_threshold:.2f}\n"
        f"Raw sequence-composition logos",
        fontsize=13,
    )
    Path(output_png).parent.mkdir(parents=True, exist_ok=True)
    figure.savefig(output_png, dpi=300)
    figure.savefig(output_pdf)
    plt.close(figure)


def plot_feature_distributions(
    rows,
    comparison_rows,
    query_condition,
    query_branch_threshold,
    output_png,
    output_pdf,
):
    high_rows = [row for row in rows if row["threshold_group"] == "threshold_high"]
    other_rows = [row for row in rows if row["threshold_group"] == "other"]
    top_rows = select_plot_rows(comparison_rows, max_features=8)
    n_features = len(top_rows)
    ncols = 2
    nrows = int(math.ceil(n_features / ncols))
    figure, axes = plt.subplots(nrows, ncols, figsize=(12.5, 2.8 * nrows), constrained_layout=True)
    axes = np.atleast_1d(axes).ravel()

    for axis, row in zip(axes, top_rows, strict=False):
        plot_feature = row["plot_feature"]
        high_values = np.array([float(entry[plot_feature]) for entry in high_rows], dtype=float)
        other_values = np.array([float(entry[plot_feature]) for entry in other_rows], dtype=float)
        boxplot = axis.boxplot(
            [other_values, high_values],
            vert=False,
            patch_artist=True,
            showfliers=False,
            widths=0.6,
            tick_labels=["Other", f"{query_condition} > {query_branch_threshold:.2f}"],
        )
        for patch, color in zip(boxplot["boxes"], ["#c7c7c7", "#d95f02"], strict=True):
            patch.set_facecolor(color)
            patch.set_alpha(0.75)
        for median in boxplot["medians"]:
            median.set_color("black")
            median.set_linewidth(1.5)
        if plot_feature == "intron_length":
            axis.set_xscale("log")
            axis.set_xlabel("nt")
        elif "fraction" in plot_feature or plot_feature.startswith("intron_gc"):
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
        f"{query_condition} branch fraction threshold groups\n"
        f"High group = pooled {query_condition} branch fraction > {query_branch_threshold:.2f}",
        fontsize=13,
    )
    Path(output_png).parent.mkdir(parents=True, exist_ok=True)
    figure.savefig(output_png, dpi=300)
    figure.savefig(output_pdf)
    plt.close(figure)


def weighted_group_fraction(rows, count_field, total_field):
    numerator = sum(float(row[count_field]) for row in rows)
    denominator = sum(float(row[total_field]) for row in rows)
    return 0.0 if denominator == 0.0 else numerator / denominator


def main():
    args = parse_args()
    if not (0.0 < args.query_branch_threshold < 1.0):
        raise ValueError("--query-branch-threshold must be between 0 and 1")
    if not (0.0 < args.reference_residual_quantile_fraction <= 0.5):
        raise ValueError("--reference-residual-quantile-fraction must be in (0, 0.5]")

    shared_rows = read_tsv_rows(args.shared_introns)
    if not shared_rows:
        raise ValueError("No shared introns available")

    classified_rows, residual_high_ids = classify_rows(
        shared_rows,
        args.query_branch_threshold,
        args.reference_residual_quantile_fraction,
    )
    output_rows = build_output_rows(
        classified_rows,
        args.genome_fasta,
        args.bp_upstream,
        args.bp_downstream,
        args.five_prime_upstream,
        args.five_prime_downstream,
        args.three_prime_upstream,
        args.three_prime_downstream,
    )

    specs = feature_specs()
    feature_comparison_rows = build_group_comparison_rows(output_rows, specs)
    residual_reference_rows = build_reference_residual_comparison(output_rows, residual_high_ids, specs)
    best_feature, best_abs_smd, mean_abs_smd = summarize_group_separation(feature_comparison_rows)
    residual_best_feature, residual_best_abs_smd, residual_mean_abs_smd = summarize_group_separation(
        residual_reference_rows
    )

    threshold_high_rows = [row for row in output_rows if row["threshold_group"] == "threshold_high"]
    other_rows = [row for row in output_rows if row["threshold_group"] == "other"]
    if not threshold_high_rows or not other_rows:
        raise ValueError("The configured query branch threshold did not produce two non-empty intron groups")
    included_logo_rows = [row for row in output_rows if row["included_in_logo"] == "1"]
    overlap_with_residual = sum(1 for row in threshold_high_rows if row["in_reference_residual_high"] == "1")

    bp_positions = list(range(-args.bp_upstream, args.bp_downstream + 1))
    five_prime_positions = list(range(-args.five_prime_upstream, args.five_prime_downstream + 1))
    three_prime_positions = list(range(-args.three_prime_upstream, args.three_prime_downstream))
    position_rows = build_position_stats(output_rows, "branchpoint_logo_sequence", bp_positions, "branchpoint_flank")
    position_rows.extend(build_position_stats(output_rows, "five_prime_ss_logo_sequence", five_prime_positions, "five_prime_ss"))
    position_rows.extend(build_position_stats(output_rows, "three_prime_ss_logo_sequence", three_prime_positions, "three_prime_ss"))

    summary_rows = [
        {
            "query_condition": args.query_condition,
            "control_condition": args.control_condition,
            "query_branch_threshold": args.query_branch_threshold,
            "shared_introns": len(output_rows),
            "threshold_high_introns": len(threshold_high_rows),
            "other_introns": len(other_rows),
            "threshold_high_fraction_of_shared_introns": len(threshold_high_rows) / len(output_rows),
            "reference_residual_high_introns": len(residual_high_ids),
            "overlap_with_reference_residual_high": overlap_with_residual,
            "overlap_fraction_of_threshold_high": 0.0 if not threshold_high_rows else overlap_with_residual / len(threshold_high_rows),
            "threshold_high_query_branch_fraction_weighted": weighted_group_fraction(
                threshold_high_rows,
                "query_branched_fragments",
                "query_anchored_fragments",
            ),
            "other_query_branch_fraction_weighted": weighted_group_fraction(
                other_rows,
                "query_branched_fragments",
                "query_anchored_fragments",
            ),
            "threshold_high_control_branch_fraction_weighted": weighted_group_fraction(
                threshold_high_rows,
                "control_branched_fragments",
                "control_anchored_fragments",
            ),
            "other_control_branch_fraction_weighted": weighted_group_fraction(
                other_rows,
                "control_branched_fragments",
                "control_anchored_fragments",
            ),
            "best_feature_by_abs_standardized_mean_difference": best_feature,
            "best_feature_abs_standardized_mean_difference": best_abs_smd,
            "mean_abs_standardized_mean_difference": mean_abs_smd,
            "reference_residual_best_feature": residual_best_feature,
            "reference_residual_best_feature_abs_standardized_mean_difference": residual_best_abs_smd,
            "reference_residual_mean_abs_standardized_mean_difference": residual_mean_abs_smd,
            "logo_eligible_introns": len(included_logo_rows),
        }
    ]

    output_rows.sort(key=lambda row: int(row["query_branch_rank"]))
    position_rows.sort(key=lambda row: (row["region"], row["threshold_group"], int(row["position"]), row["base"]))

    write_rows(args.output_summary, summary_rows, list(summary_rows[0].keys()))
    write_rows(args.output_intron_groups, output_rows, list(output_rows[0].keys()))
    write_rows(args.output_feature_group_comparison, feature_comparison_rows, list(feature_comparison_rows[0].keys()))
    write_rows(args.output_position_stats, position_rows, list(position_rows[0].keys()))
    plot_feature_distributions(
        output_rows,
        feature_comparison_rows,
        args.query_condition,
        args.query_branch_threshold,
        args.output_feature_distribution_plot_png,
        args.output_feature_distribution_plot_pdf,
    )
    plot_sequence_logos(
        position_rows,
        args.query_condition,
        args.query_branch_threshold,
        bp_positions,
        five_prime_positions,
        three_prime_positions,
        args.output_sequence_logo_plot_png,
        args.output_sequence_logo_plot_pdf,
    )
    plot_sequence_frequencies(
        position_rows,
        args.query_condition,
        args.query_branch_threshold,
        bp_positions,
        five_prime_positions,
        three_prime_positions,
        args.output_sequence_frequency_plot_png,
        args.output_sequence_frequency_plot_pdf,
    )

    print(f"Shared introns analysed: {len(output_rows)}")
    print(f"Threshold-high introns: {len(threshold_high_rows)}")
    print(f"Reference residual-high overlap: {overlap_with_residual}")
    print(f"Best threshold-group feature by |SMD|: {best_feature} ({best_abs_smd:.3f})")
    print(f"Reference residual-high best |SMD|: {residual_best_feature} ({residual_best_abs_smd:.3f})")


if __name__ == "__main__":
    main()
