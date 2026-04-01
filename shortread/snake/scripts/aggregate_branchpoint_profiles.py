#!/usr/bin/env python3

import argparse
import csv
import gzip
import math
import statistics
from collections import Counter, defaultdict
from pathlib import Path

import matplotlib
import pysam

matplotlib.use("Agg")
import matplotlib.pyplot as plt


CONDITION_COLORS = {
    "ILS": "#1f77b4",
    "DIS": "#d95f02",
}
BRANCHPOINT_MOTIF_CATEGORIES = ("canonical", "remaining")
BRANCHPOINT_MOTIF_LABELS = {
    "canonical": "Canonical branchpoints",
    "remaining": "Remaining branchpoints",
}
IUPAC_RNA_BASES = {
    "A": {"A"},
    "C": {"C"},
    "G": {"G"},
    "U": {"U"},
    "T": {"U"},
    "R": {"A", "G"},
    "Y": {"C", "U"},
    "S": {"G", "C"},
    "W": {"A", "U"},
    "K": {"G", "U"},
    "M": {"A", "C"},
    "B": {"C", "G", "U"},
    "D": {"A", "G", "U"},
    "H": {"A", "C", "U"},
    "V": {"A", "C", "G"},
    "N": {"A", "C", "G", "U"},
}

T_CRITICAL_95_BY_DF = {
    1: 12.706,
    2: 4.303,
    3: 3.182,
    4: 2.776,
    5: 2.571,
    6: 2.447,
    7: 2.365,
    8: 2.306,
    9: 2.262,
    10: 2.228,
    11: 2.201,
    12: 2.179,
    13: 2.160,
    14: 2.145,
    15: 2.131,
    16: 2.120,
    17: 2.110,
    18: 2.101,
    19: 2.093,
    20: 2.086,
    21: 2.080,
    22: 2.074,
    23: 2.069,
    24: 2.064,
    25: 2.060,
    26: 2.056,
    27: 2.052,
    28: 2.048,
    29: 2.045,
    30: 2.042,
}


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("--metaprofile", action="append", dest="metaprofiles", required=True)
    parser.add_argument("--summary", action="append", dest="summaries", required=True)
    parser.add_argument("--site-counts", action="append", dest="site_counts", required=True)
    parser.add_argument("--intron-offsets", action="append", dest="intron_offsets", required=True)
    parser.add_argument("--three-prime-coverage", action="append", dest="three_prime_coverages", required=True)
    parser.add_argument("--genome-fasta", required=True)
    parser.add_argument("--canonical-branchpoint-motif", default="YUNAY")
    parser.add_argument("--plot-upstream", type=int, default=50)
    parser.add_argument("--plot-downstream", type=int, default=10)
    parser.add_argument("--shared-min-reads", type=int, default=0)
    parser.add_argument("--output-metaprofile-by-sample", required=True)
    parser.add_argument("--output-metaprofile-by-condition", required=True)
    parser.add_argument("--output-summary-by-sample", required=True)
    parser.add_argument("--output-summary-by-condition", required=True)
    parser.add_argument("--output-shared-introns", required=True)
    parser.add_argument("--output-plot-png", required=True)
    parser.add_argument("--output-plot-pdf", required=True)
    parser.add_argument("--output-coverage-plot-png")
    parser.add_argument("--output-coverage-plot-pdf")
    parser.add_argument("--output-coverage-by-motif-by-sample", required=True)
    parser.add_argument("--output-coverage-by-motif-by-condition", required=True)
    parser.add_argument("--output-canonical-coverage-plot-png", required=True)
    parser.add_argument("--output-canonical-coverage-plot-pdf", required=True)
    parser.add_argument("--output-remaining-coverage-plot-png", required=True)
    parser.add_argument("--output-remaining-coverage-plot-pdf", required=True)
    parser.add_argument("--output-three-prime-coverage-by-sample", required=True)
    parser.add_argument("--output-three-prime-coverage-by-condition", required=True)
    parser.add_argument("--output-three-prime-coverage-plot-png", required=True)
    parser.add_argument("--output-three-prime-coverage-plot-pdf", required=True)
    return parser.parse_args()


def open_text(path, mode="rt"):
    if path.endswith(".gz"):
        return gzip.open(path, mode)
    return open(path, mode)


def iter_tsv_rows(path):
    with open_text(path) as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            yield row


def read_tsv_rows(path):
    return list(iter_tsv_rows(path))


def write_rows(path, rows, fieldnames):
    with open_text(path, "wt") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def float_mean(values):
    return 0.0 if not values else sum(values) / len(values)


def float_sd(values):
    if len(values) < 2:
        return 0.0
    return statistics.stdev(values)


def float_sem(values):
    if len(values) < 2:
        return 0.0
    return statistics.stdev(values) / math.sqrt(len(values))


def t_critical_95(sample_size):
    if sample_size < 2:
        return 0.0
    degrees_freedom = sample_size - 1
    return T_CRITICAL_95_BY_DF.get(degrees_freedom, 1.96)


def float_ci95_half_width(values):
    if len(values) < 2:
        return 0.0
    return float_sem(values) * t_critical_95(len(values))


def reverse_complement(sequence):
    return sequence.translate(str.maketrans("ACGTNacgtn", "TGCANtgcan"))[::-1]


def fetch_oriented_rna_sequence(fasta, chrom, start_1based, end_1based, strand):
    sequence = fasta.fetch(chrom, start_1based - 1, end_1based).upper()
    if strand == "-":
        sequence = reverse_complement(sequence)
    return sequence.replace("T", "U")


def motif_branchpoint_index(motif):
    motif = motif.upper().replace("T", "U")
    branchpoint_positions = [index for index, base in enumerate(motif) if base == "A"]
    if len(branchpoint_positions) != 1:
        raise ValueError(
            f"Canonical branchpoint motif must contain exactly one unambiguous A marking the branchpoint: {motif}"
        )
    return branchpoint_positions[0]


def matches_iupac_rna_motif(sequence, motif):
    sequence = sequence.upper().replace("T", "U")
    motif = motif.upper().replace("T", "U")
    if len(sequence) != len(motif):
        return False
    for observed_base, motif_base in zip(sequence, motif):
        if observed_base not in IUPAC_RNA_BASES.get(motif_base, set()):
            return False
    return True


def annotate_branchpoint_motif_categories(site_metadata, fasta_path, canonical_motif):
    canonical_motif = canonical_motif.upper().replace("T", "U")
    branchpoint_index = motif_branchpoint_index(canonical_motif)
    upstream = branchpoint_index
    downstream = len(canonical_motif) - branchpoint_index - 1

    with pysam.FastaFile(fasta_path) as fasta:
        for intron_id, metadata in site_metadata.items():
            branchpoint_position = count_value(metadata, "branchpoint_position")
            start_1based = branchpoint_position - upstream
            end_1based = branchpoint_position + downstream
            motif_sequence = ""
            if start_1based >= 1:
                try:
                    motif_sequence = fetch_oriented_rna_sequence(
                        fasta,
                        metadata["chrom"],
                        start_1based,
                        end_1based,
                        metadata["strand"],
                    )
                except ValueError:
                    motif_sequence = ""
            metadata["branchpoint_motif_sequence"] = motif_sequence
            metadata["branchpoint_motif_category"] = (
                "canonical" if matches_iupac_rna_motif(motif_sequence, canonical_motif) else "remaining"
            )


def infer_sample_name_from_path(path):
    return Path(path).name.split(".")[0]


def map_input_paths_by_sample(paths):
    sample_paths = {}
    for path in paths:
        sample = infer_sample_name_from_path(path)
        if sample in sample_paths:
            raise ValueError(f"Duplicate input detected for sample {sample}: {path}")
        sample_paths[sample] = path
    return sample_paths


def infer_sample_name(rows, path):
    samples = {row["sample"] for row in rows if "sample" in row and row["sample"]}
    if not samples:
        return Path(path).name.split(".")[0]
    if len(samples) != 1:
        raise ValueError(f"Expected one sample in {path}, found {sorted(samples)}")
    return next(iter(samples))


def count_value(row, field):
    raw_value = row.get(field, 0)
    if raw_value in {"", None}:
        return 0
    return int(float(raw_value))


def require_single_window(summary_rows, upstream_field, downstream_field, label):
    upstream_values = {count_value(row, upstream_field) for row in summary_rows if upstream_field in row}
    downstream_values = {count_value(row, downstream_field) for row in summary_rows if downstream_field in row}
    if len(upstream_values) != 1 or len(downstream_values) != 1:
        raise ValueError(f"Expected one shared {label} window across summaries")
    return next(iter(upstream_values)), next(iter(downstream_values))


def summarise_condition_rows(summary_rows, condition_order):
    grouped = defaultdict(list)
    for row in summary_rows:
        grouped[row["condition"]].append(row)

    numeric_fields = [field for field in summary_rows[0] if field not in {"sample", "condition"}]
    rows = []
    for condition in condition_order:
        entries = grouped[condition]
        condition_row = {
            "condition": condition,
            "replicate_count": len(entries),
        }
        for field in numeric_fields:
            values = [float(entry[field]) for entry in entries]
            condition_row[f"mean_{field}"] = float_mean(values)
            condition_row[f"sd_{field}"] = float_sd(values)
        rows.append(condition_row)
    return rows


def summarise_condition_profiles(metaprofile_rows, condition_order, extra_group_fields=()):
    grouped = defaultdict(lambda: defaultdict(lambda: defaultdict(list)))
    for row in metaprofile_rows:
        group_key = tuple(row[field] for field in extra_group_fields)
        grouped[group_key][row["condition"]][int(row["offset_nt"])].append(row)

    numeric_fields = [
        field
        for field in metaprofile_rows[0]
        if field not in {"sample", "condition", "offset_nt", *extra_group_fields}
    ]
    condition_rows = []
    for group_key in sorted(grouped):
        group_row = dict(zip(extra_group_fields, group_key))
        for condition in condition_order:
            offsets = sorted(grouped[group_key][condition])
            for offset in offsets:
                entries = grouped[group_key][condition][offset]
                condition_row = {
                    **group_row,
                    "condition": condition,
                    "offset_nt": offset,
                    "replicate_count": len(entries),
                }
                for field in numeric_fields:
                    values = [float(entry[field]) for entry in entries]
                    condition_row[f"mean_{field}"] = float_mean(values)
                    condition_row[f"sd_{field}"] = float_sd(values)
                    condition_row[f"sem_{field}"] = float_sem(values)
                    condition_row[f"ci95_{field}"] = float_ci95_half_width(values)
                condition_rows.append(condition_row)
    return condition_rows


def build_shared_intron_set(site_counts_by_sample, sample_order, shared_min_reads):
    qualifying_sets = []
    for sample in sample_order:
        sample_rows = site_counts_by_sample.get(sample, {})
        qualifying_sets.append(
            {
                intron_id
                for intron_id, row in sample_rows.items()
                if count_value(row, "anchored_fragments") >= shared_min_reads
            }
        )
    if not qualifying_sets:
        return set()
    return set.intersection(*qualifying_sets)


def build_shared_introns_rows(shared_introns, site_counts_by_sample, site_metadata, sample_order):
    fieldnames = [
        "intron_id",
        "gene_id",
        "gene_name",
        "transcript_id",
        "intron_number",
        "chrom",
        "strand",
        "intron_start",
        "intron_end",
        "three_prime_ss",
        "branchpoint_position",
        "branchpoint_score",
        "branchpoint_to_3ss_nt",
        "branchpoint_candidates",
        "branchpoint_motif_sequence",
        "branchpoint_motif_category",
        "min_anchored_fragments_all_samples",
    ] + [f"{sample}_anchored_fragments" for sample in sample_order]

    rows = []
    for intron_id in sorted(shared_introns, key=lambda key: (site_metadata[key]["gene_name"], key)):
        metadata = site_metadata[intron_id]
        per_sample_counts = [
            count_value(site_counts_by_sample.get(sample, {}).get(intron_id, {}), "anchored_fragments")
            for sample in sample_order
        ]
        row = {
            "intron_id": intron_id,
            "gene_id": metadata["gene_id"],
            "gene_name": metadata["gene_name"],
            "transcript_id": metadata["transcript_id"],
            "intron_number": metadata["intron_number"],
            "chrom": metadata["chrom"],
            "strand": metadata["strand"],
            "intron_start": metadata["intron_start"],
            "intron_end": metadata["intron_end"],
            "three_prime_ss": metadata["three_prime_ss"],
            "branchpoint_position": metadata["branchpoint_position"],
            "branchpoint_score": metadata["branchpoint_score"],
            "branchpoint_to_3ss_nt": metadata["branchpoint_to_3ss_nt"],
            "branchpoint_candidates": metadata["branchpoint_candidates"],
            "branchpoint_motif_sequence": metadata.get("branchpoint_motif_sequence", ""),
            "branchpoint_motif_category": metadata.get("branchpoint_motif_category", ""),
            "min_anchored_fragments_all_samples": min(per_sample_counts) if per_sample_counts else 0,
        }
        for sample, count in zip(sample_order, per_sample_counts):
            row[f"{sample}_anchored_fragments"] = count
        rows.append(row)
    return rows, fieldnames


def split_shared_introns_by_motif_category(shared_introns, site_metadata):
    shared_introns_by_category = {category: set() for category in BRANCHPOINT_MOTIF_CATEGORIES}
    for intron_id in shared_introns:
        category = site_metadata[intron_id].get("branchpoint_motif_category", "remaining")
        if category not in shared_introns_by_category:
            category = "remaining"
        shared_introns_by_category[category].add(intron_id)
    return shared_introns_by_category


def aggregate_offset_counts(intron_offset_counts, shared_introns, offset_range):
    offset_set = set(offset_range)
    total_counts = Counter()
    for intron_id in shared_introns:
        for offset, read_count in intron_offset_counts.get(intron_id, {}).items():
            if offset in offset_set:
                total_counts[offset] += read_count
    return total_counts


def aggregate_coverage_counts(intron_offset_counts, site_metadata, shared_introns, offset_range):
    coverage_counts = Counter()
    ordered_offsets = sorted(offset_range)
    if not ordered_offsets:
        return coverage_counts

    min_offset = ordered_offsets[0]

    for intron_id in shared_introns:
        intron_offsets = intron_offset_counts.get(intron_id, {})
        three_prime_offset = count_value(site_metadata[intron_id], "branchpoint_to_3ss_nt")

        # Offset tables retain all anchored 5' ends that fall within the intron,
        # including starts upstream of the plotted window.
        cumulative_fragments = sum(
            read_count
            for offset, read_count in intron_offsets.items()
            if offset < min_offset
        )
        for offset in ordered_offsets:
            cumulative_fragments += intron_offsets.get(offset, 0)
            if offset <= three_prime_offset:
                coverage_counts[offset] += cumulative_fragments

    return coverage_counts


def build_sample_metaprofile_rows(
    sample,
    condition,
    library_fragments,
    anchored_fragments,
    offset_range,
    total_offset_counts,
    total_coverage_counts,
):
    rows = []
    for offset in offset_range:
        read_count = total_offset_counts[offset]
        anchored_fraction = 0.0 if anchored_fragments == 0 else read_count / anchored_fragments
        coverage_count = total_coverage_counts[offset]
        coverage_anchored_fraction = 0.0 if anchored_fragments == 0 else coverage_count / anchored_fragments
        rows.append(
            {
                "sample": sample,
                "condition": condition,
                "offset_nt": offset,
                "read_count": read_count,
                "cpm": 0.0 if library_fragments == 0 else (read_count * 1_000_000.0 / library_fragments),
                "anchored_fraction": anchored_fraction,
                "anchored_percent": anchored_fraction * 100.0,
                "coverage_count": coverage_count,
                "coverage_cpm": 0.0
                if library_fragments == 0
                else (coverage_count * 1_000_000.0 / library_fragments),
                "coverage_anchored_fraction": coverage_anchored_fraction,
                "coverage_anchored_percent": coverage_anchored_fraction * 100.0,
            }
        )
    return rows


def build_sample_three_prime_coverage_rows(
    sample,
    condition,
    library_fragments,
    spanning_fragments,
    offset_range,
    total_coverage_counts,
):
    rows = []
    for offset in offset_range:
        coverage_count = total_coverage_counts[offset]
        coverage_spanning_fraction = 0.0 if spanning_fragments == 0 else coverage_count / spanning_fragments
        rows.append(
            {
                "sample": sample,
                "condition": condition,
                "offset_nt": offset,
                "coverage_count": coverage_count,
                "coverage_cpm": 0.0
                if library_fragments == 0
                else (coverage_count * 1_000_000.0 / library_fragments),
                "coverage_spanning_fraction": coverage_spanning_fraction,
                "coverage_spanning_percent": coverage_spanning_fraction * 100.0,
            }
        )
    return rows


def build_sample_motif_coverage_rows(
    sample,
    condition,
    motif_category,
    library_fragments,
    anchored_fragments,
    offset_range,
    total_coverage_counts,
):
    rows = []
    for offset in offset_range:
        coverage_count = total_coverage_counts[offset]
        coverage_anchored_fraction = 0.0 if anchored_fragments == 0 else coverage_count / anchored_fragments
        rows.append(
            {
                "sample": sample,
                "condition": condition,
                "branchpoint_motif_category": motif_category,
                "offset_nt": offset,
                "coverage_count": coverage_count,
                "coverage_cpm": 0.0
                if library_fragments == 0
                else (coverage_count * 1_000_000.0 / library_fragments),
                "coverage_anchored_fraction": coverage_anchored_fraction,
                "coverage_anchored_percent": coverage_anchored_fraction * 100.0,
            }
        )
    return rows


def build_sample_summary_row(
    raw_summary_row,
    site_counts,
    total_offset_counts,
    shared_introns,
    shared_min_reads,
):
    summary_row = dict(raw_summary_row)
    library_fragments = count_value(raw_summary_row, "library_fragments")
    raw_anchored_fragments = count_value(raw_summary_row, "anchored_fragments")
    raw_anchored_introns = count_value(raw_summary_row, "anchored_introns_with_reads")

    eligible_rows = [site_counts[intron_id] for intron_id in shared_introns if intron_id in site_counts]
    anchored_fragments = sum(count_value(row, "anchored_fragments") for row in eligible_rows)
    exact_branchpoint_fragments = sum(count_value(row, "exact_branchpoint_fragments") for row in eligible_rows)
    plus_one_branchpoint_fragments = sum(count_value(row, "plus_one_branchpoint_fragments") for row in eligible_rows)
    zero_or_plus_one_branchpoint_fragments = exact_branchpoint_fragments + plus_one_branchpoint_fragments
    profile_window_fragments = sum(total_offset_counts.values())

    summary_row["shared_min_reads_all_samples"] = shared_min_reads
    summary_row["shared_introns"] = len(shared_introns)
    summary_row["raw_anchored_fragments"] = raw_anchored_fragments
    summary_row["raw_anchored_introns_with_reads"] = raw_anchored_introns
    summary_row["anchored_fragments"] = anchored_fragments
    summary_row["anchored_fragments_cpm"] = (
        0.0 if library_fragments == 0 else anchored_fragments * 1_000_000.0 / library_fragments
    )
    summary_row["anchored_introns_with_reads"] = len(eligible_rows)
    summary_row["exact_branchpoint_fragments"] = exact_branchpoint_fragments
    summary_row["exact_branchpoint_cpm"] = (
        0.0 if library_fragments == 0 else exact_branchpoint_fragments * 1_000_000.0 / library_fragments
    )
    summary_row["exact_branchpoint_fraction_anchored"] = (
        0.0 if anchored_fragments == 0 else exact_branchpoint_fragments / anchored_fragments
    )
    summary_row["exact_branchpoint_percent_anchored"] = summary_row["exact_branchpoint_fraction_anchored"] * 100.0
    summary_row["plus_one_branchpoint_fragments"] = plus_one_branchpoint_fragments
    summary_row["plus_one_branchpoint_cpm"] = (
        0.0 if library_fragments == 0 else plus_one_branchpoint_fragments * 1_000_000.0 / library_fragments
    )
    summary_row["plus_one_branchpoint_fraction_anchored"] = (
        0.0 if anchored_fragments == 0 else plus_one_branchpoint_fragments / anchored_fragments
    )
    summary_row["plus_one_branchpoint_percent_anchored"] = (
        summary_row["plus_one_branchpoint_fraction_anchored"] * 100.0
    )
    summary_row["zero_or_plus_one_branchpoint_fragments"] = zero_or_plus_one_branchpoint_fragments
    summary_row["zero_or_plus_one_branchpoint_cpm"] = (
        0.0 if library_fragments == 0 else zero_or_plus_one_branchpoint_fragments * 1_000_000.0 / library_fragments
    )
    summary_row["zero_or_plus_one_branchpoint_fraction_anchored"] = (
        0.0 if anchored_fragments == 0 else zero_or_plus_one_branchpoint_fragments / anchored_fragments
    )
    summary_row["zero_or_plus_one_branchpoint_percent_anchored"] = (
        summary_row["zero_or_plus_one_branchpoint_fraction_anchored"] * 100.0
    )
    summary_row["profile_window_fragments"] = profile_window_fragments
    summary_row["profile_window_fraction_anchored"] = (
        0.0 if anchored_fragments == 0 else profile_window_fragments / anchored_fragments
    )
    summary_row["profile_window_percent_anchored"] = summary_row["profile_window_fraction_anchored"] * 100.0
    return summary_row


def plot_metaprofile_figure(
    sample_rows,
    condition_rows,
    summary_rows,
    condition_order,
    shared_min_reads,
    plot_upstream,
    plot_downstream,
    sample_value_field,
    condition_value_field,
    condition_ci95_field,
    xlabel,
    ylabel,
    title,
    output_png,
    output_pdf,
):
    profile_by_condition = defaultdict(list)
    for row in condition_rows:
        profile_by_condition[row["condition"]].append(row)

    figure, (ax_profile, ax_exact) = plt.subplots(
        1,
        2,
        figsize=(12, 4.5),
        gridspec_kw={"width_ratios": [3.4, 1.2]},
        constrained_layout=True,
    )

    for condition in condition_order:
        ordered_rows = sorted(profile_by_condition[condition], key=lambda row: int(row["offset_nt"]))
        color = CONDITION_COLORS.get(condition, "#4c4c4c")
        x_values = [int(row["offset_nt"]) for row in ordered_rows]
        y_values = [float(row[condition_value_field]) for row in ordered_rows]
        ci95_values = [float(row[condition_ci95_field]) for row in ordered_rows]
        lower = [max(y - ci95, 0.0) for y, ci95 in zip(y_values, ci95_values)]
        upper = [y + ci95 for y, ci95 in zip(y_values, ci95_values)]
        ax_profile.fill_between(x_values, lower, upper, color=color, alpha=0.16, linewidth=0)
        ax_profile.plot(
            x_values,
            y_values,
            color=color,
            linewidth=3.25,
            label=condition,
        )

    ax_profile.axvline(0, color="#4c4c4c", linestyle="--", linewidth=1)
    ax_profile.set_xlim(-plot_upstream, plot_downstream)
    ax_profile.set_xlabel(xlabel)
    ax_profile.set_ylabel(ylabel)
    if shared_min_reads > 0:
        title += f"\nShared introns with >= {shared_min_reads} anchored reads in every sample"
    ax_profile.set_title(title)
    ax_profile.legend(frameon=False)

    summary_by_condition = defaultdict(list)
    for row in summary_rows:
        summary_by_condition[row["condition"]].append(row)

    for idx, condition in enumerate(condition_order):
        entries = summary_by_condition[condition]
        values = [float(entry["zero_or_plus_one_branchpoint_percent_anchored"]) for entry in entries]
        color = CONDITION_COLORS.get(condition, "#4c4c4c")
        if values:
            if len(values) == 1:
                jitter = [idx]
            else:
                step = 0.24 / (len(values) - 1)
                jitter = [idx - 0.12 + step * i for i in range(len(values))]
            ax_exact.scatter(jitter, values, color=color, s=36, zorder=3)
            mean_value = float_mean(values)
            ci95_value = float_ci95_half_width(values)
            ax_exact.hlines(mean_value, idx - 0.18, idx + 0.18, color=color, linewidth=2.5)
            if ci95_value > 0:
                ax_exact.vlines(idx, mean_value - ci95_value, mean_value + ci95_value, color=color, linewidth=1.5)

    ax_exact.set_xticks(range(len(condition_order)))
    ax_exact.set_xticklabels(condition_order)
    ax_exact.set_ylabel("Branchpoint-proximal reads (%)")
    ax_exact.set_title("Offsets 0 / +1 fragments")

    for axis in (ax_profile, ax_exact):
        axis.spines["top"].set_visible(False)
        axis.spines["right"].set_visible(False)

    figure.savefig(output_png, dpi=300)
    figure.savefig(output_pdf)
    plt.close(figure)


def plot_single_panel_metaprofile_figure(
    sample_rows,
    condition_rows,
    condition_order,
    shared_min_reads,
    plot_upstream,
    plot_downstream,
    sample_value_field,
    condition_value_field,
    condition_ci95_field,
    xlabel,
    ylabel,
    title,
    output_png,
    output_pdf,
):
    profile_by_condition = defaultdict(list)
    for row in condition_rows:
        profile_by_condition[row["condition"]].append(row)

    ordered_offsets = sorted({int(row["offset_nt"]) for row in sample_rows})
    if not ordered_offsets:
        raise ValueError("No metaprofile rows available for plotting")

    figure, axis = plt.subplots(1, 1, figsize=(8.5, 4.5), constrained_layout=True)

    x_min = -plot_upstream
    x_max = plot_downstream
    visible_y_max = 0.0
    for condition in condition_order:
        ordered_rows = sorted(profile_by_condition[condition], key=lambda row: int(row["offset_nt"]))
        color = CONDITION_COLORS.get(condition, "#4c4c4c")
        x_values = [int(row["offset_nt"]) for row in ordered_rows]
        y_values = [float(row[condition_value_field]) for row in ordered_rows]
        ci95_values = [float(row[condition_ci95_field]) for row in ordered_rows]
        lower = [max(y - ci95, 0.0) for y, ci95 in zip(y_values, ci95_values)]
        upper = [y + ci95 for y, ci95 in zip(y_values, ci95_values)]
        visible_upper = [value for offset, value in zip(x_values, upper) if x_min <= offset <= x_max]
        visible_y_max = max(visible_y_max, max(visible_upper, default=0.0))
        axis.fill_between(x_values, lower, upper, color=color, alpha=0.16, linewidth=0)
        axis.plot(x_values, y_values, color=color, linewidth=3.25, label=condition)

    axis.axvline(0, color="#4c4c4c", linestyle="--", linewidth=1)
    axis.set_xlim(x_min, x_max)
    axis.set_ylim(0.0, 1.0 if visible_y_max == 0 else visible_y_max * 1.08)
    axis.set_xlabel(xlabel)
    axis.set_ylabel(ylabel)
    if shared_min_reads > 0:
        title += f"\nShared introns with >= {shared_min_reads} anchored reads in every sample"
    axis.set_title(title)
    axis.legend(frameon=False)
    axis.spines["top"].set_visible(False)
    axis.spines["right"].set_visible(False)

    figure.savefig(output_png, dpi=300)
    figure.savefig(output_pdf)
    plt.close(figure)


def plot_results(
    sample_rows,
    condition_rows,
    summary_rows,
    condition_order,
    shared_min_reads,
    plot_upstream,
    plot_downstream,
    output_png,
    output_pdf,
    coverage_output_png=None,
    coverage_output_pdf=None,
):
    plot_metaprofile_figure(
        sample_rows,
        condition_rows,
        summary_rows,
        condition_order,
        shared_min_reads,
        plot_upstream,
        plot_downstream,
        "anchored_percent",
        "mean_anchored_percent",
        "ci95_anchored_percent",
        "Read1 5' end offset from selected branchpoint (nt; + toward intron 3' end)",
        "Anchored shared-intron fragments (%)",
        "Branchpoint-centred 5' end metaprofile",
        output_png,
        output_pdf,
    )

    if coverage_output_png and coverage_output_pdf:
        plot_metaprofile_figure(
            sample_rows,
            condition_rows,
            summary_rows,
            condition_order,
            shared_min_reads,
            plot_upstream,
            plot_downstream,
            "coverage_anchored_percent",
            "mean_coverage_anchored_percent",
            "ci95_coverage_anchored_percent",
            "Offset from selected branchpoint (nt; + toward intron 3' end)",
            "Estimated anchored-fragment coverage (%)",
            "Branchpoint-centred fragment coverage\nAssuming anchored 3' ends align to the 3' splice site",
            coverage_output_png,
            coverage_output_pdf,
        )


def aggregate_three_prime_coverage_from_path(path, shared_introns):
    sample = infer_sample_name_from_path(path)
    totals = Counter()

    for row in iter_tsv_rows(path):
        row_sample = row.get("sample", "")
        if row_sample and row_sample != sample:
            raise ValueError(f"Expected one sample in {path}, found at least {sample!r} and {row_sample!r}")
        if row["intron_id"] not in shared_introns:
            continue
        totals[int(row["offset_nt"])] += int(row["coverage_count"])

    return sample, totals


def main():
    args = parse_args()
    if bool(args.output_coverage_plot_png) != bool(args.output_coverage_plot_pdf):
        raise ValueError("Coverage plot outputs must provide both PNG and PDF paths")

    raw_metaprofile_rows = []
    for path in args.metaprofiles:
        raw_metaprofile_rows.extend(read_tsv_rows(path))
    offset_range = sorted({int(row["offset_nt"]) for row in raw_metaprofile_rows})

    raw_summary_rows = []
    condition_order = []
    for path in args.summaries:
        rows = read_tsv_rows(path)
        if len(rows) != 1:
            raise ValueError(f"Expected exactly one summary row in {path}")
        row = rows[0]
        raw_summary_rows.append(row)
        if row["condition"] not in condition_order:
            condition_order.append(row["condition"])
    raw_summary_rows.sort(key=lambda row: (condition_order.index(row["condition"]), row["sample"]))
    sample_order = [row["sample"] for row in raw_summary_rows]
    raw_summary_by_sample = {row["sample"]: row for row in raw_summary_rows}
    three_prime_upstream, three_prime_downstream = require_single_window(
        raw_summary_rows,
        "three_prime_coverage_upstream_nt",
        "three_prime_coverage_downstream_nt",
        "3'SS coverage",
    )
    three_prime_offset_range = list(range(-three_prime_upstream, three_prime_downstream + 1))

    site_counts_by_sample = {}
    site_metadata = {}
    for path in args.site_counts:
        rows = read_tsv_rows(path)
        sample = infer_sample_name(rows, path)
        sample_rows = {}
        for row in rows:
            sample_rows[row["intron_id"]] = row
            site_metadata.setdefault(row["intron_id"], row)
        site_counts_by_sample[sample] = sample_rows
    annotate_branchpoint_motif_categories(
        site_metadata,
        args.genome_fasta,
        args.canonical_branchpoint_motif,
    )

    intron_offsets_by_sample = {}
    for path in args.intron_offsets:
        rows = read_tsv_rows(path)
        sample = infer_sample_name(rows, path)
        sample_offsets = defaultdict(Counter)
        for row in rows:
            sample_offsets[row["intron_id"]][int(row["offset_nt"])] = int(row["read_count"])
        intron_offsets_by_sample[sample] = sample_offsets

    three_prime_coverage_paths_by_sample = map_input_paths_by_sample(args.three_prime_coverages)
    missing_three_prime_samples = [sample for sample in sample_order if sample not in three_prime_coverage_paths_by_sample]
    if missing_three_prime_samples:
        raise ValueError(f"Missing 3'SS coverage inputs for samples: {', '.join(missing_three_prime_samples)}")

    shared_introns = build_shared_intron_set(site_counts_by_sample, sample_order, args.shared_min_reads)
    shared_intron_rows, shared_intron_fieldnames = build_shared_introns_rows(
        shared_introns,
        site_counts_by_sample,
        site_metadata,
        sample_order,
    )
    shared_introns_by_motif = split_shared_introns_by_motif_category(shared_introns, site_metadata)
    three_prime_coverage_by_sample = {}
    for sample in sample_order:
        _, three_prime_coverage_totals = aggregate_three_prime_coverage_from_path(
            three_prime_coverage_paths_by_sample[sample],
            shared_introns,
        )
        three_prime_coverage_by_sample[sample] = three_prime_coverage_totals

    sample_profile_rows = []
    sample_motif_coverage_rows = []
    sample_three_prime_rows = []
    sample_summary_rows = []
    for sample in sample_order:
        raw_summary_row = raw_summary_by_sample[sample]
        total_offset_counts = aggregate_offset_counts(
            intron_offsets_by_sample.get(sample, {}),
            shared_introns,
            offset_range,
        )
        total_coverage_counts = aggregate_coverage_counts(
            intron_offsets_by_sample.get(sample, {}),
            site_metadata,
            shared_introns,
            offset_range,
        )
        sample_profile_rows.extend(
            build_sample_metaprofile_rows(
                sample,
                raw_summary_row["condition"],
                count_value(raw_summary_row, "library_fragments"),
                sum(count_value(site_counts_by_sample.get(sample, {}).get(intron_id, {}), "anchored_fragments") for intron_id in shared_introns),
                offset_range,
                total_offset_counts,
                total_coverage_counts,
            )
        )
        sample_summary_rows.append(
            build_sample_summary_row(
                raw_summary_row,
                site_counts_by_sample.get(sample, {}),
                total_offset_counts,
                shared_introns,
                args.shared_min_reads,
            )
        )
        for motif_category in BRANCHPOINT_MOTIF_CATEGORIES:
            category_shared_introns = shared_introns_by_motif[motif_category]
            category_coverage_counts = aggregate_coverage_counts(
                intron_offsets_by_sample.get(sample, {}),
                site_metadata,
                category_shared_introns,
                offset_range,
            )
            sample_motif_coverage_rows.extend(
                build_sample_motif_coverage_rows(
                    sample,
                    raw_summary_row["condition"],
                    motif_category,
                    count_value(raw_summary_row, "library_fragments"),
                    sum(
                        count_value(site_counts_by_sample.get(sample, {}).get(intron_id, {}), "anchored_fragments")
                        for intron_id in category_shared_introns
                    ),
                    offset_range,
                    category_coverage_counts,
                )
            )
        sample_three_prime_rows.extend(
            build_sample_three_prime_coverage_rows(
                sample,
                raw_summary_row["condition"],
                count_value(raw_summary_row, "library_fragments"),
                three_prime_coverage_by_sample.get(sample, Counter())[0],
                three_prime_offset_range,
                three_prime_coverage_by_sample.get(sample, Counter()),
            )
        )

    sample_profile_rows.sort(
        key=lambda row: (condition_order.index(row["condition"]), row["sample"], int(row["offset_nt"]))
    )
    sample_motif_coverage_rows.sort(
        key=lambda row: (
            BRANCHPOINT_MOTIF_CATEGORIES.index(row["branchpoint_motif_category"]),
            condition_order.index(row["condition"]),
            row["sample"],
            int(row["offset_nt"]),
        )
    )
    sample_three_prime_rows.sort(
        key=lambda row: (condition_order.index(row["condition"]), row["sample"], int(row["offset_nt"]))
    )
    sample_summary_rows.sort(key=lambda row: (condition_order.index(row["condition"]), row["sample"]))

    condition_rows = summarise_condition_profiles(sample_profile_rows, condition_order)
    condition_motif_coverage_rows = summarise_condition_profiles(
        sample_motif_coverage_rows,
        condition_order,
        ("branchpoint_motif_category",),
    )
    condition_three_prime_rows = summarise_condition_profiles(sample_three_prime_rows, condition_order)
    condition_summary_rows = summarise_condition_rows(sample_summary_rows, condition_order)

    write_rows(args.output_metaprofile_by_sample, sample_profile_rows, list(sample_profile_rows[0].keys()))
    write_rows(args.output_metaprofile_by_condition, condition_rows, list(condition_rows[0].keys()))
    write_rows(
        args.output_coverage_by_motif_by_sample,
        sample_motif_coverage_rows,
        list(sample_motif_coverage_rows[0].keys()),
    )
    write_rows(
        args.output_coverage_by_motif_by_condition,
        condition_motif_coverage_rows,
        list(condition_motif_coverage_rows[0].keys()),
    )
    write_rows(
        args.output_three_prime_coverage_by_sample,
        sample_three_prime_rows,
        list(sample_three_prime_rows[0].keys()),
    )
    write_rows(
        args.output_three_prime_coverage_by_condition,
        condition_three_prime_rows,
        list(condition_three_prime_rows[0].keys()),
    )
    write_rows(args.output_summary_by_sample, sample_summary_rows, list(sample_summary_rows[0].keys()))
    write_rows(args.output_summary_by_condition, condition_summary_rows, list(condition_summary_rows[0].keys()))
    write_rows(args.output_shared_introns, shared_intron_rows, shared_intron_fieldnames)

    plot_results(
        sample_profile_rows,
        condition_rows,
        sample_summary_rows,
        condition_order,
        args.shared_min_reads,
        args.plot_upstream,
        args.plot_downstream,
        args.output_plot_png,
        args.output_plot_pdf,
        args.output_coverage_plot_png,
        args.output_coverage_plot_pdf,
    )
    for motif_category, output_png, output_pdf in (
        ("canonical", args.output_canonical_coverage_plot_png, args.output_canonical_coverage_plot_pdf),
        ("remaining", args.output_remaining_coverage_plot_png, args.output_remaining_coverage_plot_pdf),
    ):
        motif_sample_rows = [
            row for row in sample_motif_coverage_rows if row["branchpoint_motif_category"] == motif_category
        ]
        motif_condition_rows = [
            row for row in condition_motif_coverage_rows if row["branchpoint_motif_category"] == motif_category
        ]
        plot_single_panel_metaprofile_figure(
            motif_sample_rows,
            motif_condition_rows,
            condition_order,
            args.shared_min_reads,
            args.plot_upstream,
            args.plot_downstream,
            "coverage_anchored_percent",
            "mean_coverage_anchored_percent",
            "ci95_coverage_anchored_percent",
            "Offset from selected branchpoint (nt; + toward intron 3' end)",
            "Estimated anchored-fragment coverage (%)",
            "Branchpoint-centred fragment coverage\n"
            f"{BRANCHPOINT_MOTIF_LABELS[motif_category]} "
            f"({len(shared_introns_by_motif[motif_category])} shared introns)\n"
            f"Canonical motif: {args.canonical_branchpoint_motif.upper().replace('T', 'U')}",
            output_png,
            output_pdf,
        )
    plot_single_panel_metaprofile_figure(
        sample_three_prime_rows,
        condition_three_prime_rows,
        condition_order,
        args.shared_min_reads,
        args.plot_upstream,
        args.plot_downstream,
        "coverage_spanning_percent",
        "mean_coverage_spanning_percent",
        "ci95_coverage_spanning_percent",
        "Offset from intron 3' splice site (nt; + downstream of the 3'SS)",
        "3'SS-spanning fragment coverage (%)",
        "3' splice site-centred fragment coverage",
        args.output_three_prime_coverage_plot_png,
        args.output_three_prime_coverage_plot_pdf,
    )

    print(f"Shared introns retained: {len(shared_introns)}")
    print(
        "Shared introns by branchpoint motif: "
        + ", ".join(
            f"{category}={len(shared_introns_by_motif[category])}" for category in BRANCHPOINT_MOTIF_CATEGORIES
        )
    )
    print(f"Minimum anchored reads in all samples: {args.shared_min_reads}")
    print(f"Metaprofile rows aggregated: {len(sample_profile_rows)}")
    print(f"Branchpoint motif coverage rows aggregated: {len(sample_motif_coverage_rows)}")
    print(f"3'SS coverage rows aggregated: {len(sample_three_prime_rows)}")
    print(f"Summary rows aggregated: {len(sample_summary_rows)}")


if __name__ == "__main__":
    main()
