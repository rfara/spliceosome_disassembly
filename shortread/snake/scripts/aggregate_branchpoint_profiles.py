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
    parser.add_argument("--anchored-enrichment-control-condition")
    parser.add_argument("--anchored-enrichment-query-condition")
    parser.add_argument("--anchored-enrichment-denominator-scope", choices=("all", "shared"), default="all")
    parser.add_argument("--anchored-enrichment-min-log2-fold-change", type=float)
    parser.add_argument("--anchored-enrichment-max-log2-fold-change", type=float)
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
    parser.add_argument("--output-anchored-enrichment-histogram-png")
    parser.add_argument("--output-anchored-enrichment-histogram-pdf")
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

    numeric_fields = []
    metadata_fields = []
    for field in summary_rows[0]:
        if field in {"sample", "condition"}:
            continue
        non_empty_values = [entry[field] for entry in summary_rows if entry.get(field, "") not in {"", None}]
        if not non_empty_values:
            numeric_fields.append(field)
            continue
        try:
            for value in non_empty_values:
                float(value)
        except (TypeError, ValueError):
            metadata_fields.append(field)
        else:
            numeric_fields.append(field)

    rows = []
    for condition in condition_order:
        entries = grouped[condition]
        condition_row = {
            "condition": condition,
            "replicate_count": len(entries),
        }
        for field in numeric_fields:
            values = [float(entry[field]) for entry in entries if entry.get(field, "") not in {"", None}]
            if not values:
                condition_row[f"mean_{field}"] = ""
                condition_row[f"sd_{field}"] = ""
            else:
                condition_row[f"mean_{field}"] = float_mean(values)
                condition_row[f"sd_{field}"] = float_sd(values)
        for field in metadata_fields:
            values = [entry[field] for entry in entries if entry.get(field, "") not in {"", None}]
            condition_row[field] = values[0] if values else ""
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
    if shared_min_reads <= 0:
        return {
            intron_id
            for sample in sample_order
            for intron_id, row in site_counts_by_sample.get(sample, {}).items()
            if count_value(row, "anchored_fragments") > 0
        }

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


def build_condition_samples(summary_rows, condition_order):
    grouped = defaultdict(list)
    for row in summary_rows:
        grouped[row["condition"]].append(row["sample"])
    return {condition: grouped[condition] for condition in condition_order}


def pooled_anchored_fragments_for_intron(site_counts_by_sample, intron_id, sample_names):
    return sum(
        count_value(site_counts_by_sample.get(sample, {}).get(intron_id, {}), "anchored_fragments")
        for sample in sample_names
    )


def pooled_anchored_fragments_by_condition(site_counts_by_sample, condition_samples, intron_ids=None):
    totals = {}
    for condition, sample_names in condition_samples.items():
        if intron_ids is None:
            totals[condition] = sum(
                count_value(row, "anchored_fragments")
                for sample in sample_names
                for row in site_counts_by_sample.get(sample, {}).values()
            )
        else:
            totals[condition] = sum(
                count_value(site_counts_by_sample.get(sample, {}).get(intron_id, {}), "anchored_fragments")
                for sample in sample_names
                for intron_id in intron_ids
            )
    return totals


def compute_anchored_enrichment_log2_fold_change(
    control_anchored_fragments,
    query_anchored_fragments,
    control_total_anchored_fragments,
    query_total_anchored_fragments,
):
    if (
        control_anchored_fragments <= 0
        or query_anchored_fragments <= 0
        or control_total_anchored_fragments <= 0
        or query_total_anchored_fragments <= 0
    ):
        return None
    control_share = control_anchored_fragments / control_total_anchored_fragments
    query_share = query_anchored_fragments / query_total_anchored_fragments
    return math.log2(query_share / control_share)


def filter_shared_introns_by_anchored_enrichment(
    shared_introns,
    site_counts_by_sample,
    condition_samples,
    control_condition,
    query_condition,
    denominator_scope,
    min_log2_fold_change,
    max_log2_fold_change,
):
    if denominator_scope == "all":
        condition_totals = pooled_anchored_fragments_by_condition(site_counts_by_sample, condition_samples)
    elif denominator_scope == "shared":
        condition_totals = pooled_anchored_fragments_by_condition(
            site_counts_by_sample,
            condition_samples,
            shared_introns,
        )
    else:
        raise ValueError(f"Unsupported anchored-enrichment denominator scope: {denominator_scope}")
    filter_active = min_log2_fold_change is not None or max_log2_fold_change is not None

    if filter_active:
        if not control_condition or not query_condition:
            raise ValueError(
                "Anchored-enrichment filtering requires both control and query condition names"
            )
        if control_condition not in condition_samples:
            raise ValueError(f"Unknown anchored-enrichment control condition: {control_condition}")
        if query_condition not in condition_samples:
            raise ValueError(f"Unknown anchored-enrichment query condition: {query_condition}")

    intron_metrics = {}
    retained_introns = set()
    control_samples = condition_samples.get(control_condition, [])
    query_samples = condition_samples.get(query_condition, [])
    control_total = condition_totals.get(control_condition, 0)
    query_total = condition_totals.get(query_condition, 0)

    for intron_id in shared_introns:
        control_count = pooled_anchored_fragments_for_intron(site_counts_by_sample, intron_id, control_samples)
        query_count = pooled_anchored_fragments_for_intron(site_counts_by_sample, intron_id, query_samples)
        log2_fold_change = compute_anchored_enrichment_log2_fold_change(
            control_count,
            query_count,
            control_total,
            query_total,
        )
        passes_filter = True
        if min_log2_fold_change is not None and (
            log2_fold_change is None or log2_fold_change < min_log2_fold_change
        ):
            passes_filter = False
        if max_log2_fold_change is not None and (
            log2_fold_change is None or log2_fold_change > max_log2_fold_change
        ):
            passes_filter = False
        if passes_filter:
            retained_introns.add(intron_id)
        intron_metrics[intron_id] = {
            "pooled_control_anchored_fragments": control_count,
            "pooled_query_anchored_fragments": query_count,
            "anchored_enrichment_log2_fold_change": log2_fold_change,
            "passes_anchored_enrichment_filter": passes_filter,
        }

    if not filter_active:
        retained_introns = set(shared_introns)

    filter_metadata = {
        "filter_active": filter_active,
        "control_condition": control_condition or "",
        "query_condition": query_condition or "",
        "denominator_scope": denominator_scope,
        "min_log2_fold_change": min_log2_fold_change,
        "max_log2_fold_change": max_log2_fold_change,
        "pre_filter_shared_introns": len(shared_introns),
        "retained_shared_introns": len(retained_introns),
        "removed_shared_introns": len(shared_introns) - len(retained_introns),
    }
    return retained_introns, intron_metrics, filter_metadata


def build_shared_introns_rows(
    shared_introns,
    site_counts_by_sample,
    site_metadata,
    sample_order,
    anchored_enrichment_metrics,
):
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
        "pooled_control_anchored_fragments",
        "pooled_query_anchored_fragments",
        "anchored_enrichment_log2_fold_change",
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
            "pooled_control_anchored_fragments": anchored_enrichment_metrics[intron_id][
                "pooled_control_anchored_fragments"
            ],
            "pooled_query_anchored_fragments": anchored_enrichment_metrics[intron_id][
                "pooled_query_anchored_fragments"
            ],
            "anchored_enrichment_log2_fold_change": anchored_enrichment_metrics[intron_id][
                "anchored_enrichment_log2_fold_change"
            ],
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
    anchored_enrichment_filter,
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
    summary_row["shared_introns_before_anchored_enrichment_filter"] = anchored_enrichment_filter[
        "pre_filter_shared_introns"
    ]
    summary_row["shared_introns_removed_by_anchored_enrichment_filter"] = anchored_enrichment_filter[
        "removed_shared_introns"
    ]
    summary_row["anchored_enrichment_filter_applied"] = int(anchored_enrichment_filter["filter_active"])
    summary_row["anchored_enrichment_control_condition"] = anchored_enrichment_filter["control_condition"]
    summary_row["anchored_enrichment_query_condition"] = anchored_enrichment_filter["query_condition"]
    summary_row["anchored_enrichment_denominator_scope"] = anchored_enrichment_filter["denominator_scope"]
    summary_row["anchored_enrichment_min_log2_fold_change"] = (
        ""
        if anchored_enrichment_filter["min_log2_fold_change"] is None
        else anchored_enrichment_filter["min_log2_fold_change"]
    )
    summary_row["anchored_enrichment_max_log2_fold_change"] = (
        ""
        if anchored_enrichment_filter["max_log2_fold_change"] is None
        else anchored_enrichment_filter["max_log2_fold_change"]
    )
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


def format_anchored_enrichment_filter_label(anchored_enrichment_filter):
    if not anchored_enrichment_filter["filter_active"]:
        return ""

    metric = (
        f"log2({anchored_enrichment_filter['query_condition']} anchored share / "
        f"{anchored_enrichment_filter['control_condition']} anchored share)"
    )
    if anchored_enrichment_filter["denominator_scope"] == "shared":
        metric += " within shared introns"
    min_value = anchored_enrichment_filter["min_log2_fold_change"]
    max_value = anchored_enrichment_filter["max_log2_fold_change"]

    if min_value is not None and max_value is not None:
        return f"Anchored abundance filter: {min_value:g} <= {metric} <= {max_value:g}"
    if min_value is not None:
        return f"Anchored abundance filter: {metric} >= {min_value:g}"
    return f"Anchored abundance filter: {metric} <= {max_value:g}"


def plot_anchored_enrichment_histogram(
    anchored_enrichment_metrics,
    anchored_enrichment_filter,
    output_png,
    output_pdf,
):
    if not output_png and not output_pdf:
        return
    if bool(output_png) != bool(output_pdf):
        raise ValueError("Anchored-enrichment histogram outputs must provide both PNG and PDF paths")

    values = [
        metrics["anchored_enrichment_log2_fold_change"]
        for metrics in anchored_enrichment_metrics.values()
        if metrics["anchored_enrichment_log2_fold_change"] is not None
    ]
    if not values:
        raise ValueError("No finite anchored-enrichment values available for plotting")

    retained_values = [
        metrics["anchored_enrichment_log2_fold_change"]
        for metrics in anchored_enrichment_metrics.values()
        if metrics["anchored_enrichment_log2_fold_change"] is not None
        and metrics["passes_anchored_enrichment_filter"]
    ]
    filtered_values = [
        metrics["anchored_enrichment_log2_fold_change"]
        for metrics in anchored_enrichment_metrics.values()
        if metrics["anchored_enrichment_log2_fold_change"] is not None
        and not metrics["passes_anchored_enrichment_filter"]
    ]

    figure, axis = plt.subplots(1, 1, figsize=(7.5, 4.8), constrained_layout=True)
    bins = 60
    if anchored_enrichment_filter["filter_active"]:
        axis.hist(
            [filtered_values, retained_values],
            bins=bins,
            stacked=True,
            color=["#bdbdbd", "#2b8cbe"],
            label=[
                f"Filtered out (n={len(filtered_values)})",
                f"Retained (n={len(retained_values)})",
            ],
        )
    else:
        axis.hist(values, bins=bins, color="#2b8cbe", label=f"Shared introns (n={len(values)})")

    min_cutoff = anchored_enrichment_filter["min_log2_fold_change"]
    max_cutoff = anchored_enrichment_filter["max_log2_fold_change"]
    if min_cutoff is not None:
        axis.axvline(min_cutoff, color="#1f1f1f", linestyle="--", linewidth=1.5, label=f"min cutoff {min_cutoff:g}")
    if max_cutoff is not None:
        axis.axvline(max_cutoff, color="#1f1f1f", linestyle=":", linewidth=1.5, label=f"max cutoff {max_cutoff:g}")

    control_condition = anchored_enrichment_filter["control_condition"] or "control"
    query_condition = anchored_enrichment_filter["query_condition"] or "query"
    denominator_label = "shared introns" if anchored_enrichment_filter["denominator_scope"] == "shared" else "all introns"
    axis.set_xlabel(
        f"log2({query_condition} anchored share / {control_condition} anchored share; denominator: {denominator_label})"
    )
    axis.set_ylabel("Shared introns")
    axis.set_title(
        "Anchored abundance enrichment across pre-filter shared introns\n"
        f"n={len(values)} before filter; {len(retained_values)} retained"
    )
    axis.legend(frameon=False)
    axis.spines["top"].set_visible(False)
    axis.spines["right"].set_visible(False)

    figure.savefig(output_png, dpi=300)
    figure.savefig(output_pdf)
    plt.close(figure)


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
    filter_label,
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
    title_lines = [title]
    if shared_min_reads > 0:
        title_lines.append(f"Shared introns with >= {shared_min_reads} anchored reads in every sample")
    if filter_label:
        title_lines.append(filter_label)
    title = "\n".join(title_lines)
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
    filter_label,
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
    title_lines = [title]
    if shared_min_reads > 0:
        title_lines.append(f"Shared introns with >= {shared_min_reads} anchored reads in every sample")
    if filter_label:
        title_lines.append(filter_label)
    title = "\n".join(title_lines)
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
    filter_label,
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
        filter_label,
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
            filter_label,
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
    if bool(args.output_anchored_enrichment_histogram_png) != bool(args.output_anchored_enrichment_histogram_pdf):
        raise ValueError("Anchored-enrichment histogram outputs must provide both PNG and PDF paths")

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
    condition_samples = build_condition_samples(raw_summary_rows, condition_order)
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
    shared_introns, anchored_enrichment_metrics, anchored_enrichment_filter = (
        filter_shared_introns_by_anchored_enrichment(
            shared_introns,
            site_counts_by_sample,
            condition_samples,
            args.anchored_enrichment_control_condition,
            args.anchored_enrichment_query_condition,
            args.anchored_enrichment_denominator_scope,
            args.anchored_enrichment_min_log2_fold_change,
            args.anchored_enrichment_max_log2_fold_change,
        )
    )
    shared_intron_rows, shared_intron_fieldnames = build_shared_introns_rows(
        shared_introns,
        site_counts_by_sample,
        site_metadata,
        sample_order,
        anchored_enrichment_metrics,
    )
    shared_introns_by_motif = split_shared_introns_by_motif_category(shared_introns, site_metadata)
    filter_label = format_anchored_enrichment_filter_label(anchored_enrichment_filter)
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
                anchored_enrichment_filter,
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
        filter_label,
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
            filter_label,
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
        filter_label,
        args.output_three_prime_coverage_plot_png,
        args.output_three_prime_coverage_plot_pdf,
    )
    plot_anchored_enrichment_histogram(
        anchored_enrichment_metrics,
        anchored_enrichment_filter,
        args.output_anchored_enrichment_histogram_png,
        args.output_anchored_enrichment_histogram_pdf,
    )

    print(
        f"Shared introns retained: {len(shared_introns)} "
        f"(from {anchored_enrichment_filter['pre_filter_shared_introns']} before anchored-enrichment filtering)"
    )
    print(
        "Shared introns by branchpoint motif: "
        + ", ".join(
            f"{category}={len(shared_introns_by_motif[category])}" for category in BRANCHPOINT_MOTIF_CATEGORIES
        )
    )
    print(f"Minimum anchored reads in all samples: {args.shared_min_reads}")
    if anchored_enrichment_filter["filter_active"]:
        print(filter_label)
    print(f"Metaprofile rows aggregated: {len(sample_profile_rows)}")
    print(f"Branchpoint motif coverage rows aggregated: {len(sample_motif_coverage_rows)}")
    print(f"3'SS coverage rows aggregated: {len(sample_three_prime_rows)}")
    print(f"Summary rows aggregated: {len(sample_summary_rows)}")


if __name__ == "__main__":
    main()
