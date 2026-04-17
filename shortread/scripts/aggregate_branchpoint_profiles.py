#!/usr/bin/env python3

import argparse
import csv
import gzip
import math
import statistics
from collections import Counter, defaultdict
from pathlib import Path

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
    parser.add_argument("--shared-min-reads", type=int, default=0)
    parser.add_argument("--output-metaprofile-by-sample", required=True)
    parser.add_argument("--output-metaprofile-by-condition", required=True)
    parser.add_argument("--output-summary-by-sample", required=True)
    parser.add_argument("--output-summary-by-condition", required=True)
    parser.add_argument("--output-shared-introns", required=True)
    parser.add_argument("--output-three-prime-coverage-by-sample", required=True)
    parser.add_argument("--output-three-prime-coverage-by-condition", required=True)
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

def build_shared_introns_rows(
    shared_introns,
    site_counts_by_sample,
    site_metadata,
    sample_order,
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
            "min_anchored_fragments_all_samples": min(per_sample_counts) if per_sample_counts else 0,
        }
        for sample, count in zip(sample_order, per_sample_counts):
            row[f"{sample}_anchored_fragments"] = count
        rows.append(row)
    return rows, fieldnames


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
    three_prime_coverage_by_sample = {}
    for sample in sample_order:
        _, three_prime_coverage_totals = aggregate_three_prime_coverage_from_path(
            three_prime_coverage_paths_by_sample[sample],
            shared_introns,
        )
        three_prime_coverage_by_sample[sample] = three_prime_coverage_totals

    sample_profile_rows = []
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
    sample_three_prime_rows.sort(
        key=lambda row: (condition_order.index(row["condition"]), row["sample"], int(row["offset_nt"]))
    )
    sample_summary_rows.sort(key=lambda row: (condition_order.index(row["condition"]), row["sample"]))

    condition_rows = summarise_condition_profiles(sample_profile_rows, condition_order)
    condition_three_prime_rows = summarise_condition_profiles(sample_three_prime_rows, condition_order)
    condition_summary_rows = summarise_condition_rows(sample_summary_rows, condition_order)

    write_rows(args.output_metaprofile_by_sample, sample_profile_rows, list(sample_profile_rows[0].keys()))
    write_rows(args.output_metaprofile_by_condition, condition_rows, list(condition_rows[0].keys()))
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

    print(
        f"Shared introns retained: {len(shared_introns)}"
    )
    print(f"Minimum anchored reads in all samples: {args.shared_min_reads}")
    print(f"Metaprofile rows aggregated: {len(sample_profile_rows)}")
    print(f"3'SS coverage rows aggregated: {len(sample_three_prime_rows)}")
    print(f"Summary rows aggregated: {len(sample_summary_rows)}")


if __name__ == "__main__":
    main()
