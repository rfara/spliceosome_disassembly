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
    parser.add_argument("--gene-counts", action="append", dest="gene_counts", required=True)
    parser.add_argument("--coverage", action="append", dest="coverages", required=True)
    parser.add_argument("--summary", action="append", dest="summaries", required=True)
    parser.add_argument("--shared-min-total-reads", type=int, required=True)
    parser.add_argument("--output-metaprofile-by-sample", required=True)
    parser.add_argument("--output-metaprofile-by-condition", required=True)
    parser.add_argument("--output-summary-by-sample", required=True)
    parser.add_argument("--output-summary-by-condition", required=True)
    parser.add_argument("--output-shared-genes", required=True)
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


def infer_sample_name_from_path(path):
    return Path(path).name.split(".")[0]


def infer_sample_name(rows, path):
    samples = {row["sample"] for row in rows if row.get("sample")}
    if not samples:
        return infer_sample_name_from_path(path)
    if len(samples) != 1:
        raise ValueError(f"Expected one sample in {path}, found {sorted(samples)}")
    return next(iter(samples))


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
    if raw_value in {"", None, "NA"}:
        return 0
    return int(float(raw_value))


def parse_optional_float(raw_value):
    if raw_value in {"", None, "NA", "nan", "NaN"}:
        return None
    value = float(raw_value)
    if not math.isfinite(value):
        return None
    return value


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
    return T_CRITICAL_95_BY_DF.get(sample_size - 1, 1.96)


def float_ci95_half_width(values):
    if len(values) < 2:
        return 0.0
    return float_sem(values) * t_critical_95(len(values))


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
            values = [parse_optional_float(entry[field]) for entry in entries]
            values = [value for value in values if value is not None]
            if not values:
                condition_row[f"mean_{field}"] = "NA"
                condition_row[f"sd_{field}"] = "NA"
                condition_row[f"sem_{field}"] = "NA"
                condition_row[f"ci95_{field}"] = "NA"
            else:
                condition_row[f"mean_{field}"] = float_mean(values)
                condition_row[f"sd_{field}"] = float_sd(values)
                condition_row[f"sem_{field}"] = float_sem(values)
                condition_row[f"ci95_{field}"] = float_ci95_half_width(values)
        rows.append(condition_row)
    return rows


def summarise_condition_profiles(metaprofile_rows, condition_order):
    grouped = defaultdict(lambda: defaultdict(list))
    for row in metaprofile_rows:
        grouped[row["condition"]][int(row["offset_nt"])].append(row)

    numeric_fields = [field for field in metaprofile_rows[0] if field not in {"sample", "condition", "offset_nt"}]
    condition_rows = []
    for condition in condition_order:
        for offset in sorted(grouped[condition]):
            entries = grouped[condition][offset]
            row = {
                "condition": condition,
                "offset_nt": offset,
                "replicate_count": len(entries),
            }
            for field in numeric_fields:
                values = [float(entry[field]) for entry in entries]
                row[f"mean_{field}"] = float_mean(values)
                row[f"sd_{field}"] = float_sd(values)
                row[f"sem_{field}"] = float_sem(values)
                row[f"ci95_{field}"] = float_ci95_half_width(values)
            condition_rows.append(row)
    return condition_rows


def require_single_window(summary_rows, upstream_field, downstream_field, label):
    upstream_values = {count_value(row, upstream_field) for row in summary_rows if upstream_field in row}
    downstream_values = {count_value(row, downstream_field) for row in summary_rows if downstream_field in row}
    if len(upstream_values) != 1 or len(downstream_values) != 1:
        raise ValueError(f"Expected one shared {label} window across summaries")
    return next(iter(upstream_values)), next(iter(downstream_values))


def build_shared_gene_set(gene_counts_by_sample, sample_order, shared_min_total_reads):
    qualifying_sets = []
    for sample in sample_order:
        sample_rows = gene_counts_by_sample.get(sample, {})
        qualifying_sets.append(
            {
                gene_id
                for gene_id, row in sample_rows.items()
                if count_value(row, "total_gene_fragments") >= shared_min_total_reads
            }
        )
    if not qualifying_sets:
        return set()
    return set.intersection(*qualifying_sets)


def build_shared_gene_rows(shared_genes, gene_counts_by_sample, gene_metadata, sample_order):
    fieldnames = [
        "gene_id",
        "gene_name",
        "transcript_id",
        "chrom",
        "strand",
        "transcript_start",
        "transcript_end",
        "three_prime_end",
        "terminal_exon_start",
        "terminal_exon_end",
        "min_total_gene_fragments_all_samples",
    ]
    for sample in sample_order:
        fieldnames.extend(
            [
                f"{sample}_total_gene_fragments",
                f"{sample}_mrna_fragments",
                f"{sample}_intronic_fragments",
            ]
        )

    rows = []
    for gene_id in sorted(shared_genes, key=lambda value: (gene_metadata[value]["gene_name"], value)):
        metadata = gene_metadata[gene_id]
        per_sample_totals = [
            count_value(gene_counts_by_sample.get(sample, {}).get(gene_id, {}), "total_gene_fragments")
            for sample in sample_order
        ]
        row = {
            "gene_id": gene_id,
            "gene_name": metadata["gene_name"],
            "transcript_id": metadata["transcript_id"],
            "chrom": metadata["chrom"],
            "strand": metadata["strand"],
            "transcript_start": metadata["transcript_start"],
            "transcript_end": metadata["transcript_end"],
            "three_prime_end": metadata["three_prime_end"],
            "terminal_exon_start": metadata["terminal_exon_start"],
            "terminal_exon_end": metadata["terminal_exon_end"],
            "min_total_gene_fragments_all_samples": min(per_sample_totals) if per_sample_totals else 0,
        }
        for sample in sample_order:
            sample_row = gene_counts_by_sample.get(sample, {}).get(gene_id, {})
            row[f"{sample}_total_gene_fragments"] = count_value(sample_row, "total_gene_fragments")
            row[f"{sample}_mrna_fragments"] = count_value(sample_row, "mrna_fragments")
            row[f"{sample}_intronic_fragments"] = count_value(sample_row, "intronic_fragments")
        rows.append(row)
    return rows, fieldnames


def aggregate_coverage_from_path(path, shared_genes):
    sample = infer_sample_name_from_path(path)
    totals = Counter()
    for row in iter_tsv_rows(path):
        row_sample = row.get("sample", "")
        if row_sample and row_sample != sample:
            raise ValueError(f"Expected one sample in {path}, found at least {sample!r} and {row_sample!r}")
        if row["gene_id"] not in shared_genes:
            continue
        totals[int(row["offset_nt"])] += int(row["coverage_count"])
    return sample, totals


def build_sample_metaprofile_rows(sample, condition, library_fragments, total_gene_fragments, offset_range, coverage_counts):
    rows = []
    for offset in offset_range:
        coverage_count = coverage_counts[offset]
        rows.append(
            {
                "sample": sample,
                "condition": condition,
                "offset_nt": offset,
                "coverage_count": coverage_count,
                "coverage_cpm": 0.0
                if library_fragments == 0
                else coverage_count * 1_000_000.0 / library_fragments,
                "coverage_fraction_gene_reads": 0.0
                if total_gene_fragments == 0
                else coverage_count / total_gene_fragments,
                "coverage_percent_gene_reads": 0.0
                if total_gene_fragments == 0
                else coverage_count * 100.0 / total_gene_fragments,
            }
        )
    return rows


def ratio_or_na(numerator, denominator):
    if denominator == 0:
        return "NA"
    return numerator / denominator


def build_sample_summary_row(raw_summary_row, gene_counts, shared_genes, shared_min_total_reads, coverage_counts):
    summary_row = dict(raw_summary_row)
    library_fragments = count_value(raw_summary_row, "library_fragments")
    raw_assigned_gene_fragments = count_value(raw_summary_row, "assigned_gene_fragments")
    raw_genes_with_reads = count_value(raw_summary_row, "genes_with_reads")
    raw_mrna_fragments = count_value(raw_summary_row, "mrna_fragments")
    raw_intronic_fragments = count_value(raw_summary_row, "intronic_fragments")

    eligible_rows = [gene_counts[gene_id] for gene_id in shared_genes if gene_id in gene_counts]
    total_gene_fragments = sum(count_value(row, "total_gene_fragments") for row in eligible_rows)
    mrna_fragments = sum(count_value(row, "mrna_fragments") for row in eligible_rows)
    intronic_fragments = sum(count_value(row, "intronic_fragments") for row in eligible_rows)

    summary_row["shared_min_total_reads_all_samples"] = shared_min_total_reads
    summary_row["shared_genes"] = len(shared_genes)
    summary_row["raw_assigned_gene_fragments"] = raw_assigned_gene_fragments
    summary_row["raw_genes_with_reads"] = raw_genes_with_reads
    summary_row["raw_mrna_fragments"] = raw_mrna_fragments
    summary_row["raw_intronic_fragments"] = raw_intronic_fragments
    summary_row["assigned_gene_fragments"] = total_gene_fragments
    summary_row["assigned_gene_fragments_cpm"] = 0.0
    if library_fragments > 0:
        summary_row["assigned_gene_fragments_cpm"] = total_gene_fragments * 1_000_000.0 / library_fragments
    summary_row["genes_with_reads"] = len(eligible_rows)
    summary_row["mrna_fragments"] = mrna_fragments
    summary_row["mrna_fragments_cpm"] = 0.0 if library_fragments == 0 else mrna_fragments * 1_000_000.0 / library_fragments
    summary_row["mrna_genes_with_reads"] = sum(count_value(row, "mrna_fragments") > 0 for row in eligible_rows)
    summary_row["mrna_fraction_gene_reads"] = 0.0 if total_gene_fragments == 0 else mrna_fragments / total_gene_fragments
    summary_row["mrna_percent_gene_reads"] = (
        0.0 if total_gene_fragments == 0 else mrna_fragments * 100.0 / total_gene_fragments
    )
    summary_row["intronic_fragments"] = intronic_fragments
    summary_row["intronic_fragments_cpm"] = (
        0.0 if library_fragments == 0 else intronic_fragments * 1_000_000.0 / library_fragments
    )
    summary_row["intronic_genes_with_reads"] = sum(count_value(row, "intronic_fragments") > 0 for row in eligible_rows)
    summary_row["intronic_fraction_gene_reads"] = (
        0.0 if total_gene_fragments == 0 else intronic_fragments / total_gene_fragments
    )
    summary_row["intronic_percent_gene_reads"] = (
        0.0 if total_gene_fragments == 0 else intronic_fragments * 100.0 / total_gene_fragments
    )
    summary_row["intronic_to_mrna_ratio"] = ratio_or_na(intronic_fragments, mrna_fragments)
    summary_row["profiled_coverage_positions"] = sum(coverage_counts.values())
    return summary_row


def main():
    args = parse_args()

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

    profile_upstream, profile_downstream = require_single_window(
        raw_summary_rows,
        "metaprofile_upstream_nt",
        "metaprofile_downstream_nt",
        "pre-mRNA metaprofile",
    )
    offset_range = list(range(-profile_upstream, profile_downstream + 1))

    gene_counts_by_sample = {}
    gene_metadata = {}
    for path in args.gene_counts:
        rows = read_tsv_rows(path)
        sample = infer_sample_name(rows, path)
        sample_rows = {}
        for row in rows:
            sample_rows[row["gene_id"]] = row
            gene_metadata.setdefault(row["gene_id"], row)
        gene_counts_by_sample[sample] = sample_rows

    coverage_paths_by_sample = map_input_paths_by_sample(args.coverages)
    missing_coverage_samples = [sample for sample in sample_order if sample not in coverage_paths_by_sample]
    if missing_coverage_samples:
        raise ValueError(f"Missing coverage inputs for samples: {', '.join(missing_coverage_samples)}")

    shared_genes = build_shared_gene_set(
        gene_counts_by_sample,
        sample_order,
        args.shared_min_total_reads,
    )
    if not shared_genes:
        raise ValueError("No shared genes retained for pre-mRNA aggregation")

    shared_gene_rows, shared_gene_fieldnames = build_shared_gene_rows(
        shared_genes,
        gene_counts_by_sample,
        gene_metadata,
        sample_order,
    )

    coverage_by_sample = {}
    for sample in sample_order:
        _, totals = aggregate_coverage_from_path(coverage_paths_by_sample[sample], shared_genes)
        coverage_by_sample[sample] = totals

    sample_metaprofile_rows = []
    sample_summary_rows = []
    for sample in sample_order:
        raw_summary_row = raw_summary_by_sample[sample]
        condition = raw_summary_row["condition"]
        coverage_counts = coverage_by_sample.get(sample, Counter())
        sample_summary_rows.append(
            build_sample_summary_row(
                raw_summary_row,
                gene_counts_by_sample.get(sample, {}),
                shared_genes,
                args.shared_min_total_reads,
                coverage_counts,
            )
        )
        sample_metaprofile_rows.extend(
            build_sample_metaprofile_rows(
                sample,
                condition,
                count_value(raw_summary_row, "library_fragments"),
                sample_summary_rows[-1]["assigned_gene_fragments"],
                offset_range,
                coverage_counts,
            )
        )

    sample_summary_rows.sort(key=lambda row: (condition_order.index(row["condition"]), row["sample"]))
    sample_metaprofile_rows.sort(key=lambda row: (condition_order.index(row["condition"]), row["sample"], row["offset_nt"]))

    condition_metaprofile_rows = summarise_condition_profiles(sample_metaprofile_rows, condition_order)
    condition_summary_rows = summarise_condition_rows(sample_summary_rows, condition_order)

    write_rows(
        args.output_metaprofile_by_sample,
        sample_metaprofile_rows,
        list(sample_metaprofile_rows[0].keys()),
    )
    write_rows(
        args.output_metaprofile_by_condition,
        condition_metaprofile_rows,
        list(condition_metaprofile_rows[0].keys()),
    )
    write_rows(
        args.output_summary_by_sample,
        sample_summary_rows,
        list(sample_summary_rows[0].keys()),
    )
    write_rows(
        args.output_summary_by_condition,
        condition_summary_rows,
        list(condition_summary_rows[0].keys()),
    )
    write_rows(args.output_shared_genes, shared_gene_rows, shared_gene_fieldnames)


if __name__ == "__main__":
    main()
