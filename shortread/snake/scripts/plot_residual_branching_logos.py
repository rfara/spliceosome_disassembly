#!/usr/bin/env python3

import argparse
import csv
import gzip
import math
from collections import Counter
from pathlib import Path

import matplotlib
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
GROUPS = [
    ("residual_high", "Residual-high introns"),
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
    parser.add_argument("--residual-quantile-fraction", type=float, default=0.1)
    parser.add_argument("--bp-upstream", type=int, default=10)
    parser.add_argument("--bp-downstream", type=int, default=10)
    parser.add_argument("--five-prime-upstream", type=int, default=5)
    parser.add_argument("--five-prime-downstream", type=int, default=10)
    parser.add_argument("--three-prime-upstream", type=int, default=20)
    parser.add_argument("--three-prime-downstream", type=int, default=3)
    parser.add_argument("--output-intron-groups", required=True)
    parser.add_argument("--output-position-stats", required=True)
    parser.add_argument("--output-plot-png", required=True)
    parser.add_argument("--output-plot-pdf", required=True)
    parser.add_argument("--output-frequency-plot-png", required=True)
    parser.add_argument("--output-frequency-plot-pdf", required=True)
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


def classify_rows(rows, quantile_fraction):
    ordered = sorted(
        rows,
        key=lambda row: (
            -float(row["baseline_residual_query_branch_fraction"]),
            -int(row["query_anchored_fragments"]),
            row["gene_name"],
            row["intron_id"],
        ),
    )
    group_size = max(1, int(round(len(ordered) * quantile_fraction)))
    classified = []
    for rank, row in enumerate(ordered, start=1):
        classified.append(
            {
                **row,
                "residual_rank": rank,
                "residual_group": "residual_high" if rank <= group_size else "other",
            }
        )
    return classified, group_size


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


def extract_three_prime_sequence(row, fasta, upstream, downstream):
    intron_start = int(row["intron_start"])
    intron_end = int(row["intron_end"])
    if row["strand"] == "+":
        start_1based = intron_end - upstream + 1
        end_1based = intron_end + downstream
    else:
        start_1based = intron_start - downstream
        end_1based = intron_start + upstream - 1

    sequence = fetch_oriented_sequence(
        fasta,
        row["chrom"],
        start_1based,
        end_1based,
        row["strand"],
    )
    expected_length = upstream + downstream
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

    sequence = fetch_oriented_sequence(
        fasta,
        row["chrom"],
        start_1based,
        end_1based,
        row["strand"],
    )
    expected_length = upstream + downstream + 1
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
                    "residual_rank": row["residual_rank"],
                    "residual_group": row["residual_group"],
                    "intron_id": row["intron_id"],
                    "gene_id": row["gene_id"],
                    "gene_name": row["gene_name"],
                    "transcript_id": row["transcript_id"],
                    "intron_number": row["intron_number"],
                    "chrom": row["chrom"],
                    "strand": row["strand"],
                    "intron_start": row["intron_start"],
                    "intron_end": row["intron_end"],
                    "branchpoint_position": row["branchpoint_position"],
                    "branchpoint_to_3ss_nt": row["branchpoint_to_3ss_nt"],
                    "control_branch_fraction": row["control_branch_fraction"],
                    "query_branch_fraction": row["query_branch_fraction"],
                    "baseline_expected_query_branch_fraction": row["baseline_expected_query_branch_fraction"],
                    "baseline_residual_query_branch_fraction": row["baseline_residual_query_branch_fraction"],
                    "branchpoint_logo_sequence": bp_sequence or "",
                    "five_prime_ss_logo_sequence": five_prime_sequence or "",
                    "three_prime_ss_logo_sequence": three_prime_sequence or "",
                    "included_in_logo": "1" if not exclusion_reasons else "0",
                    "exclusion_reason": ",".join(exclusion_reasons),
                }
            )
    return output_rows


def build_position_stats(rows, sequence_field, positions, region_name):
    grouped = []
    for group_name, _ in GROUPS:
        sequences = [row[sequence_field] for row in rows if row["residual_group"] == group_name and row["included_in_logo"] == "1"]
        if not sequences:
            continue
        for index, position in enumerate(positions):
            counts = Counter(sequence[index] for sequence in sequences)
            total = len(sequences)
            fractions = {base: counts.get(base, 0) / total for base in DNA_BASES}
            entropy = -sum(
                fraction * math.log2(fraction)
                for fraction in fractions.values()
                if fraction > 0.0
            )
            information = max(0.0, math.log2(len(DNA_BASES)) - entropy)
            for base in DNA_BASES:
                grouped.append(
                    {
                        "region": region_name,
                        "residual_group": group_name,
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
    return grouped


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
        PathPatch(
            glyph,
            transform=transform + axis.transData,
            facecolor=BASE_COLORS[base],
            edgecolor="none",
        )
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
    heights_by_position = {
        position: {base: 0.0 for base in DNA_BASES}
        for position in positions
    }
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
    if reference_mode == "center":
        reference_x = positions.index(0) + 0.5
    else:
        reference_x = positions.index(0)
    axis.axvline(reference_x, color="black", linewidth=1.0, linestyle="--", alpha=0.45)
    axis.spines["top"].set_visible(False)
    axis.spines["right"].set_visible(False)


def plot_results(
    position_rows,
    query_condition,
    control_condition,
    quantile_fraction,
    bp_positions,
    five_prime_positions,
    three_prime_positions,
    output_png,
    output_pdf,
):
    figure, axes = plt.subplots(3, 2, figsize=(15, 11), constrained_layout=True)

    for column, (group_name, group_label) in enumerate(GROUPS):
        bp_rows = [
            row for row in position_rows if row["region"] == "branchpoint_flank" and row["residual_group"] == group_name
        ]
        five_prime_rows = [
            row for row in position_rows if row["region"] == "five_prime_ss" and row["residual_group"] == group_name
        ]
        three_prime_rows = [
            row for row in position_rows if row["region"] == "three_prime_ss" and row["residual_group"] == group_name
        ]
        plot_logo(
            axes[0, column],
            bp_rows,
            bp_positions,
            f"{group_label}\nBranchpoint flank",
            "Position relative to branchpoint A (nt)",
            reference_mode="center",
        )
        plot_logo(
            axes[1, column],
            five_prime_rows,
            five_prime_positions,
            f"{group_label}\n5' splice site window",
            "Position relative to first intronic base at the 5'SS (nt)",
            reference_mode="between",
        )
        plot_logo(
            axes[2, column],
            three_prime_rows,
            three_prime_positions,
            f"{group_label}\n3' splice site window",
            "Position relative to first exonic base at the 3'SS (nt)",
            reference_mode="between",
        )

    figure.suptitle(
        f"{query_condition} residual-high introns relative to {control_condition}\n"
        f"Residual-high set defined as the top {quantile_fraction * 100:.1f}% of shared introns by baseline residual",
        fontsize=13,
    )
    Path(output_png).parent.mkdir(parents=True, exist_ok=True)
    figure.savefig(output_png, dpi=300)
    figure.savefig(output_pdf)
    plt.close(figure)


def plot_frequency_results(
    position_rows,
    query_condition,
    control_condition,
    quantile_fraction,
    bp_positions,
    five_prime_positions,
    three_prime_positions,
    output_png,
    output_pdf,
):
    figure, axes = plt.subplots(3, 2, figsize=(15, 11), constrained_layout=True)

    for column, (group_name, group_label) in enumerate(GROUPS):
        bp_rows = [
            row for row in position_rows if row["region"] == "branchpoint_flank" and row["residual_group"] == group_name
        ]
        five_prime_rows = [
            row for row in position_rows if row["region"] == "five_prime_ss" and row["residual_group"] == group_name
        ]
        three_prime_rows = [
            row for row in position_rows if row["region"] == "three_prime_ss" and row["residual_group"] == group_name
        ]
        plot_logo_generic(
            axes[0, column],
            bp_rows,
            bp_positions,
            f"{group_label}\nBranchpoint flank",
            "Position relative to branchpoint A (nt)",
            reference_mode="center",
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
            reference_mode="between",
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
            reference_mode="between",
            height_field="base_fraction",
            y_label="Base frequency",
            y_limit=1.02,
        )

    figure.suptitle(
        f"{query_condition} residual-high introns relative to {control_condition}\n"
        f"Residual-high set defined as the top {quantile_fraction * 100:.1f}% of shared introns by baseline residual\n"
        f"Raw sequence-composition logos",
        fontsize=13,
    )
    Path(output_png).parent.mkdir(parents=True, exist_ok=True)
    figure.savefig(output_png, dpi=300)
    figure.savefig(output_pdf)
    plt.close(figure)


def main():
    args = parse_args()
    if not (0.0 < args.residual_quantile_fraction <= 0.5):
        raise ValueError("--residual-quantile-fraction must be in (0, 0.5]")
    if args.bp_upstream < 0 or args.bp_downstream < 0:
        raise ValueError("Branchpoint windows must be non-negative")
    if args.five_prime_upstream <= 0 or args.five_prime_downstream <= 0:
        raise ValueError("5'SS windows must be positive")
    if args.three_prime_upstream <= 0 or args.three_prime_downstream <= 0:
        raise ValueError("3'SS windows must be positive")

    classified_rows, top_group_size = classify_rows(
        read_tsv_rows(args.shared_introns),
        args.residual_quantile_fraction,
    )
    if not classified_rows:
        raise ValueError("No shared introns were available for logo generation")

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
    included_rows = [row for row in output_rows if row["included_in_logo"] == "1"]
    if not included_rows:
        raise ValueError("No introns retained valid sequence windows for logo generation")

    bp_positions = list(range(-args.bp_upstream, args.bp_downstream + 1))
    five_prime_positions = list(range(-args.five_prime_upstream, args.five_prime_downstream + 1))
    three_prime_positions = list(range(-args.three_prime_upstream, args.three_prime_downstream))
    position_rows = build_position_stats(
        output_rows,
        "branchpoint_logo_sequence",
        bp_positions,
        "branchpoint_flank",
    )
    position_rows.extend(
        build_position_stats(
            output_rows,
            "five_prime_ss_logo_sequence",
            five_prime_positions,
            "five_prime_ss",
        )
    )
    position_rows.extend(
        build_position_stats(
            output_rows,
            "three_prime_ss_logo_sequence",
            three_prime_positions,
            "three_prime_ss",
        )
    )

    output_rows.sort(key=lambda row: int(row["residual_rank"]))
    position_rows.sort(key=lambda row: (row["region"], row["residual_group"], int(row["position"]), row["base"]))

    write_rows(args.output_intron_groups, output_rows, list(output_rows[0].keys()))
    write_rows(args.output_position_stats, position_rows, list(position_rows[0].keys()))
    plot_results(
        position_rows,
        args.query_condition,
        args.control_condition,
        args.residual_quantile_fraction,
        bp_positions,
        five_prime_positions,
        three_prime_positions,
        args.output_plot_png,
        args.output_plot_pdf,
    )
    plot_frequency_results(
        position_rows,
        args.query_condition,
        args.control_condition,
        args.residual_quantile_fraction,
        bp_positions,
        five_prime_positions,
        three_prime_positions,
        args.output_frequency_plot_png,
        args.output_frequency_plot_pdf,
    )

    group_counts = Counter(row["residual_group"] for row in output_rows)
    included_counts = Counter(row["residual_group"] for row in included_rows)
    print(f"Shared introns ranked: {len(output_rows)}")
    print(f"Residual-high introns: {group_counts['residual_high']} (target top group size {top_group_size})")
    print(f"Residual-high introns with valid windows: {included_counts['residual_high']}")
    print(f"Other introns with valid windows: {included_counts['other']}")


if __name__ == "__main__":
    main()
