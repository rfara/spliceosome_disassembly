#!/usr/bin/env python3

import argparse
import csv
import gzip
import math
import statistics
from collections import Counter, defaultdict
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt


CONDITION_COLORS = {
    "ILS": "#1f77b4",
    "DIS": "#d95f02",
}

EVENT_SPECS = [
    ("Mismatch", "mismatch_percent_coverage", "mean_mismatch_percent_coverage", "sem_mismatch_percent_coverage"),
    ("Deletion", "deletion_percent_coverage", "mean_deletion_percent_coverage", "sem_deletion_percent_coverage"),
    ("Insertion", "insertion_percent_coverage", "mean_insertion_percent_coverage", "sem_insertion_percent_coverage"),
]


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("--summary", action="append", dest="summaries", required=True)
    parser.add_argument("--site-counts", action="append", dest="site_counts", required=True)
    parser.add_argument("--position-counts", action="append", dest="position_counts", required=True)
    parser.add_argument("--shared-min-reads", type=int, default=0)
    parser.add_argument("--blacklist-min-samples", type=int, default=0)
    parser.add_argument("--blacklist-min-traversing-fragments", type=int, default=0)
    parser.add_argument("--blacklist-deletion-percent-threshold", type=float, default=-1.0)
    parser.add_argument("--blacklist-insertion-percent-threshold", type=float, default=-1.0)
    parser.add_argument("--max-total-indel-percent-any-sample", type=float, default=-1.0)
    parser.add_argument("--blacklist-single-offset-min-samples", type=int, default=0)
    parser.add_argument("--blacklist-single-offset-min-coverage", type=int, default=0)
    parser.add_argument("--blacklist-single-offset-deletion-percent-threshold", type=float, default=-1.0)
    parser.add_argument("--blacklist-single-offset-insertion-percent-threshold", type=float, default=-1.0)
    parser.add_argument("--plot-upstream", type=int, default=None)
    parser.add_argument("--plot-downstream", type=int, default=None)
    parser.add_argument("--output-metaprofile-by-sample", required=True)
    parser.add_argument("--output-metaprofile-by-condition", required=True)
    parser.add_argument("--output-summary-by-sample", required=True)
    parser.add_argument("--output-summary-by-condition", required=True)
    parser.add_argument("--output-shared-introns", required=True)
    parser.add_argument("--output-plot-png", required=True)
    parser.add_argument("--output-plot-pdf", required=True)
    parser.add_argument("--output-blacklist")
    parser.add_argument("--output-filter-summary")
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


def read_single_tsv_row(path):
    row_iter = iter_tsv_rows(path)
    try:
        row = next(row_iter)
    except StopIteration as exc:
        raise ValueError(f"No rows found in {path}") from exc
    if next(row_iter, None) is not None:
        raise ValueError(f"Expected exactly one row in {path}")
    return row


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


def infer_sample_name_from_path(path):
    return Path(path).name.split(".")[0]


def load_site_counts(path, site_metadata):
    sample = None
    sample_rows = {}
    for row in iter_tsv_rows(path):
        row_sample = row.get("sample", "")
        if sample is None:
            sample = row_sample or infer_sample_name_from_path(path)
        elif row_sample and row_sample != sample:
            raise ValueError(f"Expected one sample in {path}, found at least {sample!r} and {row_sample!r}")
        sample_rows[row["intron_id"]] = row
        site_metadata.setdefault(row["intron_id"], row)
    if sample is None:
        sample = infer_sample_name_from_path(path)
    return sample, sample_rows


def map_input_paths_by_sample(paths):
    sample_paths = {}
    for path in paths:
        sample = infer_sample_name_from_path(path)
        if sample in sample_paths:
            raise ValueError(f"Duplicate input detected for sample {sample}: {path}")
        sample_paths[sample] = path
    return sample_paths


def count_value(row, field):
    raw_value = row.get(field, 0)
    if raw_value in {"", None}:
        return 0
    return int(float(raw_value))


def float_value(row, field):
    raw_value = row.get(field, 0.0)
    if raw_value in {"", None}:
        return 0.0
    return float(raw_value)


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


def summarise_condition_profiles(metaprofile_rows, condition_order):
    grouped = defaultdict(lambda: defaultdict(list))
    for row in metaprofile_rows:
        grouped[row["condition"]][int(row["offset_nt"])].append(row)

    numeric_fields = [field for field in metaprofile_rows[0] if field not in {"sample", "condition", "offset_nt"}]
    condition_rows = []
    for condition in condition_order:
        offsets = sorted(grouped[condition])
        for offset in offsets:
            entries = grouped[condition][offset]
            condition_row = {
                "condition": condition,
                "offset_nt": offset,
                "replicate_count": len(entries),
            }
            for field in numeric_fields:
                values = [float(entry[field]) for entry in entries]
                condition_row[f"mean_{field}"] = float_mean(values)
                condition_row[f"sd_{field}"] = float_sd(values)
                condition_row[f"sem_{field}"] = float_sem(values)
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
                if count_value(row, "traversing_fragments") >= shared_min_reads
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
        "min_traversing_fragments_all_samples",
    ]
    fieldnames.extend(f"{sample}_anchored_fragments" for sample in sample_order)
    fieldnames.extend(f"{sample}_traversing_fragments" for sample in sample_order)

    rows = []
    for intron_id in sorted(shared_introns, key=lambda key: (site_metadata[key]["gene_name"], key)):
        metadata = site_metadata[intron_id]
        traversing_counts = [
            count_value(site_counts_by_sample.get(sample, {}).get(intron_id, {}), "traversing_fragments")
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
            "min_traversing_fragments_all_samples": min(traversing_counts) if traversing_counts else 0,
        }
        for sample in sample_order:
            row[f"{sample}_anchored_fragments"] = count_value(
                site_counts_by_sample.get(sample, {}).get(intron_id, {}),
                "anchored_fragments",
            )
            row[f"{sample}_traversing_fragments"] = count_value(
                site_counts_by_sample.get(sample, {}).get(intron_id, {}),
                "traversing_fragments",
            )
        rows.append(row)
    return rows, fieldnames


def build_single_offset_indel_metrics_for_path(path, min_coverage):
    sample = None
    sample_metrics = {}
    for row in iter_tsv_rows(path):
        row_sample = row.get("sample", "")
        if sample is None:
            sample = row_sample or infer_sample_name_from_path(path)
        elif row_sample and row_sample != sample:
            raise ValueError(f"Expected one sample in {path}, found at least {sample!r} and {row_sample!r}")

        intron_metrics = sample_metrics.setdefault(
            row["intron_id"],
            {
                "max_deletion_percent": 0.0,
                "max_deletion_offset": "",
                "max_insertion_percent": 0.0,
                "max_insertion_offset": "",
            },
        )
        coverage_count = int(row["coverage_count"])
        if coverage_count == 0 or coverage_count < min_coverage:
            continue

        offset = int(row["offset_nt"])
        deletion_percent = int(row["deletion_count"]) * 100.0 / coverage_count
        insertion_percent = int(row["insertion_count"]) * 100.0 / coverage_count
        if deletion_percent > intron_metrics["max_deletion_percent"]:
            intron_metrics["max_deletion_percent"] = deletion_percent
            intron_metrics["max_deletion_offset"] = offset
        if insertion_percent > intron_metrics["max_insertion_percent"]:
            intron_metrics["max_insertion_percent"] = insertion_percent
            intron_metrics["max_insertion_offset"] = offset

    if sample is None:
        sample = infer_sample_name_from_path(path)
    return sample, sample_metrics


def build_readthrough_blacklist(
    site_counts_by_sample,
    single_offset_indel_metrics_by_sample,
    site_metadata,
    sample_order,
    min_samples,
    min_traversing_fragments,
    deletion_percent_threshold,
    insertion_percent_threshold,
    max_total_indel_percent_any_sample,
    single_offset_min_samples,
    single_offset_min_coverage,
    single_offset_deletion_percent_threshold,
    single_offset_insertion_percent_threshold,
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
        "samples_meeting_min_traversing",
        "recurrent_high_deletion_burden",
        "recurrent_high_insertion_burden",
        "extreme_total_indel_burden_any_sample",
        "recurrent_single_offset_high_deletion_spike",
        "recurrent_single_offset_high_insertion_spike",
        "deletion_flagged_sample_count",
        "insertion_flagged_sample_count",
        "single_offset_deletion_flagged_sample_count",
        "single_offset_insertion_flagged_sample_count",
        "deletion_flagged_samples",
        "insertion_flagged_samples",
        "single_offset_deletion_flagged_samples",
        "single_offset_insertion_flagged_samples",
        "max_deletion_percent_any_sample",
        "max_insertion_percent_any_sample",
        "max_total_indel_percent_any_sample",
        "max_single_offset_deletion_percent_any_sample",
        "max_single_offset_deletion_offset_any_sample",
        "max_single_offset_deletion_sample",
        "max_single_offset_insertion_percent_any_sample",
        "max_single_offset_insertion_offset_any_sample",
        "max_single_offset_insertion_sample",
        "total_deletion_events_all_samples",
        "total_insertion_events_all_samples",
        "blacklist_reasons",
    ]
    for sample in sample_order:
        fieldnames.extend(
            [
                f"{sample}_traversing_fragments",
                f"{sample}_deletion_percent",
                f"{sample}_insertion_percent",
                f"{sample}_total_indel_percent",
            ]
        )

    blacklist_rows = []
    blacklist_introns = set()
    reason_counter = Counter()

    for intron_id, metadata in site_metadata.items():
        deletion_flagged_samples = []
        insertion_flagged_samples = []
        single_offset_deletion_flagged_samples = []
        single_offset_insertion_flagged_samples = []
        max_deletion_percent = 0.0
        max_insertion_percent = 0.0
        max_total_indel_percent = 0.0
        max_single_offset_deletion_percent = 0.0
        max_single_offset_deletion_offset = ""
        max_single_offset_deletion_sample = ""
        max_single_offset_insertion_percent = 0.0
        max_single_offset_insertion_offset = ""
        max_single_offset_insertion_sample = ""
        samples_meeting_min_traversing = 0
        total_deletion_events = 0
        total_insertion_events = 0
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
        }

        for sample in sample_order:
            sample_row = site_counts_by_sample.get(sample, {}).get(intron_id, {})
            single_offset_metrics = single_offset_indel_metrics_by_sample.get(sample, {}).get(intron_id, {})
            traversing_fragments = count_value(sample_row, "traversing_fragments")
            deletion_percent = float_value(sample_row, "deletion_events_per_100_covered_positions")
            insertion_percent = float_value(sample_row, "insertion_events_per_100_covered_positions")
            total_indel_percent = deletion_percent + insertion_percent
            single_offset_deletion_percent = float(single_offset_metrics.get("max_deletion_percent", 0.0))
            single_offset_deletion_offset = single_offset_metrics.get("max_deletion_offset", "")
            single_offset_insertion_percent = float(single_offset_metrics.get("max_insertion_percent", 0.0))
            single_offset_insertion_offset = single_offset_metrics.get("max_insertion_offset", "")
            total_deletion_events += count_value(sample_row, "deletion_events")
            total_insertion_events += count_value(sample_row, "insertion_events")
            max_deletion_percent = max(max_deletion_percent, deletion_percent)
            max_insertion_percent = max(max_insertion_percent, insertion_percent)
            if single_offset_deletion_percent > max_single_offset_deletion_percent:
                max_single_offset_deletion_percent = single_offset_deletion_percent
                max_single_offset_deletion_offset = single_offset_deletion_offset
                max_single_offset_deletion_sample = sample
            if single_offset_insertion_percent > max_single_offset_insertion_percent:
                max_single_offset_insertion_percent = single_offset_insertion_percent
                max_single_offset_insertion_offset = single_offset_insertion_offset
                max_single_offset_insertion_sample = sample
            if traversing_fragments >= min_traversing_fragments:
                samples_meeting_min_traversing += 1
                max_total_indel_percent = max(max_total_indel_percent, total_indel_percent)
                if deletion_percent_threshold >= 0 and deletion_percent >= deletion_percent_threshold:
                    deletion_flagged_samples.append(sample)
                if insertion_percent_threshold >= 0 and insertion_percent >= insertion_percent_threshold:
                    insertion_flagged_samples.append(sample)
            if (
                single_offset_deletion_percent_threshold >= 0
                and single_offset_deletion_percent >= single_offset_deletion_percent_threshold
            ):
                single_offset_deletion_flagged_samples.append(sample)
            if (
                single_offset_insertion_percent_threshold >= 0
                and single_offset_insertion_percent >= single_offset_insertion_percent_threshold
            ):
                single_offset_insertion_flagged_samples.append(sample)

            row[f"{sample}_traversing_fragments"] = traversing_fragments
            row[f"{sample}_deletion_percent"] = deletion_percent
            row[f"{sample}_insertion_percent"] = insertion_percent
            row[f"{sample}_total_indel_percent"] = total_indel_percent

        reasons = []
        recurrent_high_deletion = (
            deletion_percent_threshold >= 0
            and min_samples > 0
            and len(deletion_flagged_samples) >= min_samples
        )
        recurrent_high_insertion = (
            insertion_percent_threshold >= 0
            and min_samples > 0
            and len(insertion_flagged_samples) >= min_samples
        )
        extreme_total_indel = (
            max_total_indel_percent_any_sample >= 0
            and max_total_indel_percent >= max_total_indel_percent_any_sample
        )
        recurrent_single_offset_high_deletion_spike = (
            single_offset_deletion_percent_threshold >= 0
            and single_offset_min_samples > 0
            and len(single_offset_deletion_flagged_samples) >= single_offset_min_samples
        )
        recurrent_single_offset_high_insertion_spike = (
            single_offset_insertion_percent_threshold >= 0
            and single_offset_min_samples > 0
            and len(single_offset_insertion_flagged_samples) >= single_offset_min_samples
        )

        if recurrent_high_deletion:
            reasons.append("recurrent_high_deletion_burden")
            reason_counter["recurrent_high_deletion_burden"] += 1
        if recurrent_high_insertion:
            reasons.append("recurrent_high_insertion_burden")
            reason_counter["recurrent_high_insertion_burden"] += 1
        if extreme_total_indel:
            reasons.append("extreme_total_indel_burden_any_sample")
            reason_counter["extreme_total_indel_burden_any_sample"] += 1
        if recurrent_single_offset_high_deletion_spike:
            reasons.append("recurrent_single_offset_high_deletion_spike")
            reason_counter["recurrent_single_offset_high_deletion_spike"] += 1
        if recurrent_single_offset_high_insertion_spike:
            reasons.append("recurrent_single_offset_high_insertion_spike")
            reason_counter["recurrent_single_offset_high_insertion_spike"] += 1

        if not reasons:
            continue

        blacklist_introns.add(intron_id)
        row.update(
            {
                "samples_meeting_min_traversing": samples_meeting_min_traversing,
                "recurrent_high_deletion_burden": int(recurrent_high_deletion),
                "recurrent_high_insertion_burden": int(recurrent_high_insertion),
                "extreme_total_indel_burden_any_sample": int(extreme_total_indel),
                "recurrent_single_offset_high_deletion_spike": int(recurrent_single_offset_high_deletion_spike),
                "recurrent_single_offset_high_insertion_spike": int(recurrent_single_offset_high_insertion_spike),
                "deletion_flagged_sample_count": len(deletion_flagged_samples),
                "insertion_flagged_sample_count": len(insertion_flagged_samples),
                "single_offset_deletion_flagged_sample_count": len(single_offset_deletion_flagged_samples),
                "single_offset_insertion_flagged_sample_count": len(single_offset_insertion_flagged_samples),
                "deletion_flagged_samples": ",".join(deletion_flagged_samples),
                "insertion_flagged_samples": ",".join(insertion_flagged_samples),
                "single_offset_deletion_flagged_samples": ",".join(single_offset_deletion_flagged_samples),
                "single_offset_insertion_flagged_samples": ",".join(single_offset_insertion_flagged_samples),
                "max_deletion_percent_any_sample": max_deletion_percent,
                "max_insertion_percent_any_sample": max_insertion_percent,
                "max_total_indel_percent_any_sample": max_total_indel_percent,
                "max_single_offset_deletion_percent_any_sample": max_single_offset_deletion_percent,
                "max_single_offset_deletion_offset_any_sample": max_single_offset_deletion_offset,
                "max_single_offset_deletion_sample": max_single_offset_deletion_sample,
                "max_single_offset_insertion_percent_any_sample": max_single_offset_insertion_percent,
                "max_single_offset_insertion_offset_any_sample": max_single_offset_insertion_offset,
                "max_single_offset_insertion_sample": max_single_offset_insertion_sample,
                "total_deletion_events_all_samples": total_deletion_events,
                "total_insertion_events_all_samples": total_insertion_events,
                "blacklist_reasons": ",".join(reasons),
            }
        )
        blacklist_rows.append(row)

    blacklist_rows.sort(
        key=lambda row: (
            -max(
                float(row["max_total_indel_percent_any_sample"]),
                float(row["max_single_offset_deletion_percent_any_sample"]),
                float(row["max_single_offset_insertion_percent_any_sample"]),
            ),
            -int(row["deletion_flagged_sample_count"]),
            -int(row["insertion_flagged_sample_count"]),
            row["gene_name"],
            row["intron_id"],
        )
    )
    return blacklist_introns, blacklist_rows, fieldnames, reason_counter


def build_filter_summary_row(
    sample_order,
    raw_shared_introns,
    filtered_shared_introns,
    blacklist_introns,
    blacklist_rows,
    reason_counter,
    args,
):
    blacklisted_shared_introns = raw_shared_introns.intersection(blacklist_introns)
    blacklist_only_nonshared = blacklist_introns.difference(raw_shared_introns)
    row = {
        "shared_min_reads_all_samples": args.shared_min_reads,
        "raw_shared_introns": len(raw_shared_introns),
        "filtered_shared_introns": len(filtered_shared_introns),
        "blacklist_introns_total": len(blacklist_introns),
        "blacklisted_shared_introns": len(blacklisted_shared_introns),
        "blacklisted_nonshared_introns": len(blacklist_only_nonshared),
        "blacklist_min_samples": args.blacklist_min_samples,
        "blacklist_min_traversing_fragments": args.blacklist_min_traversing_fragments,
        "blacklist_deletion_percent_threshold": args.blacklist_deletion_percent_threshold,
        "blacklist_insertion_percent_threshold": args.blacklist_insertion_percent_threshold,
        "max_total_indel_percent_any_sample": args.max_total_indel_percent_any_sample,
        "blacklist_single_offset_min_samples": args.blacklist_single_offset_min_samples,
        "blacklist_single_offset_min_coverage": args.blacklist_single_offset_min_coverage,
        "blacklist_single_offset_deletion_percent_threshold": (
            args.blacklist_single_offset_deletion_percent_threshold
        ),
        "blacklist_single_offset_insertion_percent_threshold": (
            args.blacklist_single_offset_insertion_percent_threshold
        ),
        "recurrent_high_deletion_burden_introns": reason_counter["recurrent_high_deletion_burden"],
        "recurrent_high_insertion_burden_introns": reason_counter["recurrent_high_insertion_burden"],
        "extreme_total_indel_burden_any_sample_introns": reason_counter["extreme_total_indel_burden_any_sample"],
        "recurrent_single_offset_high_deletion_spike_introns": (
            reason_counter["recurrent_single_offset_high_deletion_spike"]
        ),
        "recurrent_single_offset_high_insertion_spike_introns": (
            reason_counter["recurrent_single_offset_high_insertion_spike"]
        ),
    }
    for sample in sample_order:
        row[f"{sample}_shared_introns_after_filter"] = sum(1 for intron_id in filtered_shared_introns if intron_id)
    return row


def aggregate_position_counts_from_path(path, shared_introns):
    sample = None
    totals = defaultdict(Counter)
    observed_offsets = set()
    for row in iter_tsv_rows(path):
        row_sample = row.get("sample", "")
        if sample is None:
            sample = row_sample or infer_sample_name_from_path(path)
        elif row_sample and row_sample != sample:
            raise ValueError(f"Expected one sample in {path}, found at least {sample!r} and {row_sample!r}")

        offset = int(row["offset_nt"])
        observed_offsets.add(offset)
        if row["intron_id"] not in shared_introns:
            continue

        counts = totals[offset]
        counts["coverage_count"] += int(row["coverage_count"])
        counts["mismatch_count"] += int(row["mismatch_count"])
        counts["deletion_count"] += int(row["deletion_count"])
        counts["insertion_count"] += int(row["insertion_count"])

    if sample is None:
        sample = infer_sample_name_from_path(path)
    return sample, totals, observed_offsets


def build_sample_metaprofile_rows(sample, condition, traversing_fragments, offset_range, total_counts):
    rows = []
    for offset in offset_range:
        counts = total_counts[offset]
        coverage_count = counts["coverage_count"]
        mismatch_count = counts["mismatch_count"]
        deletion_count = counts["deletion_count"]
        insertion_count = counts["insertion_count"]
        rows.append(
            {
                "sample": sample,
                "condition": condition,
                "offset_nt": offset,
                "coverage_count": coverage_count,
                "coverage_fraction_traversing": 0.0
                if traversing_fragments == 0
                else coverage_count / traversing_fragments,
                "coverage_percent_traversing": 0.0
                if traversing_fragments == 0
                else coverage_count * 100.0 / traversing_fragments,
                "mismatch_count": mismatch_count,
                "mismatch_fraction_coverage": 0.0 if coverage_count == 0 else mismatch_count / coverage_count,
                "mismatch_percent_coverage": 0.0 if coverage_count == 0 else mismatch_count * 100.0 / coverage_count,
                "deletion_count": deletion_count,
                "deletion_fraction_coverage": 0.0 if coverage_count == 0 else deletion_count / coverage_count,
                "deletion_percent_coverage": 0.0 if coverage_count == 0 else deletion_count * 100.0 / coverage_count,
                "insertion_count": insertion_count,
                "insertion_fraction_coverage": 0.0 if coverage_count == 0 else insertion_count / coverage_count,
                "insertion_percent_coverage": 0.0 if coverage_count == 0 else insertion_count * 100.0 / coverage_count,
            }
        )
    return rows


def build_sample_summary_row(raw_summary_row, site_counts, total_counts, shared_introns, shared_min_reads):
    summary_row = dict(raw_summary_row)
    library_fragments = count_value(raw_summary_row, "library_fragments")
    raw_anchored_fragments = count_value(raw_summary_row, "anchored_fragments")
    raw_traversing_fragments = count_value(raw_summary_row, "traversing_fragments")
    raw_traversing_introns = count_value(raw_summary_row, "traversing_introns_with_reads")

    eligible_rows = [site_counts[intron_id] for intron_id in shared_introns if intron_id in site_counts]
    anchored_fragments = sum(count_value(row, "anchored_fragments") for row in eligible_rows)
    traversing_fragments = sum(count_value(row, "traversing_fragments") for row in eligible_rows)
    profiled_coverage_positions = sum(total_counts[offset]["coverage_count"] for offset in total_counts)
    mismatch_events = sum(total_counts[offset]["mismatch_count"] for offset in total_counts)
    deletion_events = sum(total_counts[offset]["deletion_count"] for offset in total_counts)
    insertion_events = sum(total_counts[offset]["insertion_count"] for offset in total_counts)

    summary_row["shared_min_reads_all_samples"] = shared_min_reads
    summary_row["shared_introns"] = len(shared_introns)
    summary_row["raw_anchored_fragments"] = raw_anchored_fragments
    summary_row["raw_traversing_fragments"] = raw_traversing_fragments
    summary_row["raw_traversing_introns_with_reads"] = raw_traversing_introns
    summary_row["anchored_fragments"] = anchored_fragments
    summary_row["anchored_fragments_cpm"] = 0.0 if library_fragments == 0 else anchored_fragments * 1_000_000.0 / library_fragments
    summary_row["traversing_fragments"] = traversing_fragments
    summary_row["traversing_fragments_cpm"] = (
        0.0 if library_fragments == 0 else traversing_fragments * 1_000_000.0 / library_fragments
    )
    summary_row["traversing_introns_with_reads"] = len(eligible_rows)
    summary_row["traversing_fraction_anchored"] = 0.0 if anchored_fragments == 0 else traversing_fragments / anchored_fragments
    summary_row["traversing_percent_anchored"] = summary_row["traversing_fraction_anchored"] * 100.0
    summary_row["profiled_coverage_positions"] = profiled_coverage_positions
    summary_row["mismatch_events"] = mismatch_events
    summary_row["mismatch_events_per_100_covered_positions"] = (
        0.0 if profiled_coverage_positions == 0 else mismatch_events * 100.0 / profiled_coverage_positions
    )
    summary_row["deletion_events"] = deletion_events
    summary_row["deletion_events_per_100_covered_positions"] = (
        0.0 if profiled_coverage_positions == 0 else deletion_events * 100.0 / profiled_coverage_positions
    )
    summary_row["insertion_events"] = insertion_events
    summary_row["insertion_events_per_100_covered_positions"] = (
        0.0 if profiled_coverage_positions == 0 else insertion_events * 100.0 / profiled_coverage_positions
    )
    return summary_row


def plot_results(
    sample_rows,
    condition_rows,
    condition_order,
    shared_min_reads,
    output_png,
    output_pdf,
    blacklisted_shared_introns=0,
    plot_upstream=None,
    plot_downstream=None,
):
    sample_profiles = defaultdict(list)
    for row in sample_rows:
        sample_profiles[row["sample"]].append(row)

    condition_profiles = defaultdict(list)
    for row in condition_rows:
        condition_profiles[row["condition"]].append(row)

    ordered_offsets = sorted({int(row["offset_nt"]) for row in sample_rows})
    if not ordered_offsets:
        raise ValueError("No metaprofile rows available for plotting")

    figure, axes = plt.subplots(1, len(EVENT_SPECS), figsize=(15, 4.8), constrained_layout=True, sharex=True)
    if len(EVENT_SPECS) == 1:
        axes = [axes]

    x_min = ordered_offsets[0]
    x_max = ordered_offsets[-1]
    if plot_upstream is not None:
        x_min = max(x_min, -plot_upstream)
    if plot_downstream is not None:
        x_max = min(x_max, plot_downstream)

    for axis, (title, sample_field, condition_field, sem_field) in zip(axes, EVENT_SPECS):
        visible_y_max = 0.0
        for sample in sorted(sample_profiles):
            rows = [
                row
                for row in sorted(sample_profiles[sample], key=lambda row: int(row["offset_nt"]))
                if x_min <= int(row["offset_nt"]) <= x_max
            ]
            if not rows:
                continue
            condition = rows[0]["condition"]
            color = CONDITION_COLORS.get(condition, "#4c4c4c")
            axis.plot(
                [int(row["offset_nt"]) for row in rows],
                [float(row[sample_field]) for row in rows],
                color=color,
                alpha=0.25,
                linewidth=1.0,
            )

        for condition in condition_order:
            rows = [
                row
                for row in sorted(condition_profiles[condition], key=lambda row: int(row["offset_nt"]))
                if x_min <= int(row["offset_nt"]) <= x_max
            ]
            if not rows:
                continue
            color = CONDITION_COLORS.get(condition, "#4c4c4c")
            x_values = [int(row["offset_nt"]) for row in rows]
            y_values = [float(row[condition_field]) for row in rows]
            sem_values = [float(row[sem_field]) for row in rows]
            lower = [max(y - sem, 0.0) for y, sem in zip(y_values, sem_values)]
            upper = [y + sem for y, sem in zip(y_values, sem_values)]
            visible_y_max = max(visible_y_max, max(upper, default=0.0))

            axis.fill_between(x_values, lower, upper, color=color, alpha=0.12, linewidth=0)
            axis.plot(x_values, y_values, color=color, linewidth=2.5, label=condition)

        axis.axvline(0, color="#4c4c4c", linestyle="--", linewidth=1)
        axis.set_xlim(x_min, x_max)
        axis.set_ylim(0.0, 1.0 if visible_y_max == 0 else visible_y_max * 1.08)
        axis.set_xlabel("Offset from selected branchpoint (nt; + toward intron 3' end)")
        axis.set_ylabel("Event frequency among traversing-read coverage (%)")
        axis.set_title(title)
        axis.spines["top"].set_visible(False)
        axis.spines["right"].set_visible(False)

    legend_handles, legend_labels = axes[0].get_legend_handles_labels()
    axes[0].legend(legend_handles, legend_labels, frameon=False, loc="upper left")

    title = "Branchpoint readthrough-associated sequence changes"
    if shared_min_reads > 0:
        title += f"\nShared introns with >= {shared_min_reads} traversing reads in every sample"
    if blacklisted_shared_introns > 0:
        title += f"\nExcluded {blacklisted_shared_introns} high-indel introns"
    figure.suptitle(title, y=1.03)
    figure.savefig(output_png, dpi=300, bbox_inches="tight")
    figure.savefig(output_pdf, bbox_inches="tight")
    plt.close(figure)


def main():
    args = parse_args()
    if bool(args.output_blacklist) != bool(args.output_filter_summary):
        raise ValueError("Blacklist outputs must provide both blacklist and filter-summary paths")

    raw_summary_rows = []
    condition_order = []
    for path in args.summaries:
        row = read_single_tsv_row(path)
        raw_summary_rows.append(row)
        if row["condition"] not in condition_order:
            condition_order.append(row["condition"])
    raw_summary_rows.sort(key=lambda row: (condition_order.index(row["condition"]), row["sample"]))
    sample_order = [row["sample"] for row in raw_summary_rows]
    raw_summary_by_sample = {row["sample"]: row for row in raw_summary_rows}

    site_counts_by_sample = {}
    site_metadata = {}
    for path in args.site_counts:
        sample, sample_rows = load_site_counts(path, site_metadata)
        site_counts_by_sample[sample] = sample_rows
    position_count_paths_by_sample = map_input_paths_by_sample(args.position_counts)
    missing_site_counts = [sample for sample in sample_order if sample not in site_counts_by_sample]
    if missing_site_counts:
        raise ValueError(f"Missing site-count inputs for samples: {', '.join(missing_site_counts)}")
    missing_position_counts = [sample for sample in sample_order if sample not in position_count_paths_by_sample]
    if missing_position_counts:
        raise ValueError(f"Missing position-count inputs for samples: {', '.join(missing_position_counts)}")

    raw_shared_introns = build_shared_intron_set(site_counts_by_sample, sample_order, args.shared_min_reads)
    blacklist_introns = set()
    blacklist_rows = []
    blacklist_fieldnames = []
    reason_counter = Counter()
    if args.output_blacklist and args.output_filter_summary:
        single_offset_indel_metrics_by_sample = {}
        for sample in sample_order:
            loaded_sample, sample_metrics = build_single_offset_indel_metrics_for_path(
                position_count_paths_by_sample[sample],
                args.blacklist_single_offset_min_coverage,
            )
            if loaded_sample != sample:
                raise ValueError(
                    f"Position-count input sample mismatch for {position_count_paths_by_sample[sample]}: "
                    f"expected {sample}, found {loaded_sample}"
                )
            single_offset_indel_metrics_by_sample[sample] = sample_metrics

        blacklist_introns, blacklist_rows, blacklist_fieldnames, reason_counter = build_readthrough_blacklist(
            site_counts_by_sample,
            single_offset_indel_metrics_by_sample,
            site_metadata,
            sample_order,
            args.blacklist_min_samples,
            args.blacklist_min_traversing_fragments,
            args.blacklist_deletion_percent_threshold,
            args.blacklist_insertion_percent_threshold,
            args.max_total_indel_percent_any_sample,
            args.blacklist_single_offset_min_samples,
            args.blacklist_single_offset_min_coverage,
            args.blacklist_single_offset_deletion_percent_threshold,
            args.blacklist_single_offset_insertion_percent_threshold,
        )
    shared_introns = raw_shared_introns.difference(blacklist_introns)
    if not shared_introns:
        raise ValueError("No shared introns retained for branchpoint readthrough event aggregation")

    shared_intron_rows, shared_intron_fieldnames = build_shared_introns_rows(
        shared_introns,
        site_counts_by_sample,
        site_metadata,
        sample_order,
    )

    sample_total_counts = {}
    offset_range = set()
    for sample in sample_order:
        loaded_sample, total_counts, observed_offsets = aggregate_position_counts_from_path(
            position_count_paths_by_sample[sample],
            shared_introns,
        )
        if loaded_sample != sample:
            raise ValueError(
                f"Position-count input sample mismatch for {position_count_paths_by_sample[sample]}: "
                f"expected {sample}, found {loaded_sample}"
            )
        sample_total_counts[sample] = total_counts
        offset_range.update(observed_offsets)
    offset_range = sorted(offset_range)

    sample_profile_rows = []
    sample_summary_rows = []
    for sample in sample_order:
        raw_summary_row = raw_summary_by_sample[sample]
        total_counts = sample_total_counts[sample]
        traversing_fragments = sum(
            count_value(site_counts_by_sample.get(sample, {}).get(intron_id, {}), "traversing_fragments")
            for intron_id in shared_introns
        )
        sample_profile_rows.extend(
            build_sample_metaprofile_rows(
                sample,
                raw_summary_row["condition"],
                traversing_fragments,
                offset_range,
                total_counts,
            )
        )
        sample_summary_rows.append(
            build_sample_summary_row(
                raw_summary_row,
                site_counts_by_sample.get(sample, {}),
                total_counts,
                shared_introns,
                args.shared_min_reads,
            )
        )

    sample_profile_rows.sort(
        key=lambda row: (condition_order.index(row["condition"]), row["sample"], int(row["offset_nt"]))
    )
    sample_summary_rows.sort(key=lambda row: (condition_order.index(row["condition"]), row["sample"]))

    condition_rows = summarise_condition_profiles(sample_profile_rows, condition_order)
    condition_summary_rows = summarise_condition_rows(sample_summary_rows, condition_order)

    write_rows(args.output_metaprofile_by_sample, sample_profile_rows, list(sample_profile_rows[0].keys()))
    write_rows(args.output_metaprofile_by_condition, condition_rows, list(condition_rows[0].keys()))
    write_rows(args.output_summary_by_sample, sample_summary_rows, list(sample_summary_rows[0].keys()))
    write_rows(args.output_summary_by_condition, condition_summary_rows, list(condition_summary_rows[0].keys()))
    write_rows(args.output_shared_introns, shared_intron_rows, shared_intron_fieldnames)
    if args.output_blacklist and args.output_filter_summary:
        filter_summary_row = build_filter_summary_row(
            sample_order,
            raw_shared_introns,
            shared_introns,
            blacklist_introns,
            blacklist_rows,
            reason_counter,
            args,
        )
        write_rows(args.output_blacklist, blacklist_rows, blacklist_fieldnames)
        write_rows(args.output_filter_summary, [filter_summary_row], list(filter_summary_row.keys()))

    plot_results(
        sample_profile_rows,
        condition_rows,
        condition_order,
        args.shared_min_reads,
        args.output_plot_png,
        args.output_plot_pdf,
        len(raw_shared_introns.intersection(blacklist_introns)),
        args.plot_upstream,
        args.plot_downstream,
    )

    print(f"Raw shared introns retained before blacklist: {len(raw_shared_introns)}")
    print(f"Blacklisted introns total: {len(blacklist_introns)}")
    print(f"Blacklisted shared introns: {len(raw_shared_introns.intersection(blacklist_introns))}")
    print(f"Shared introns retained after blacklist: {len(shared_introns)}")
    print(f"Minimum traversing reads in all samples: {args.shared_min_reads}")
    print(f"Metaprofile rows aggregated: {len(sample_profile_rows)}")
    print(f"Summary rows aggregated: {len(sample_summary_rows)}")


if __name__ == "__main__":
    main()
