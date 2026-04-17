#!/usr/bin/env python3

import argparse
import csv
import gzip
import math
import statistics
from collections import Counter, defaultdict
from dataclasses import dataclass

import pysam


INTRON_INDEX_BIN_SIZE = 1000
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
}


@dataclass(frozen=True)
class SelectedIntron:
    intron_id: str
    transcript_id: str
    intron_number: int
    gene_id: str
    gene_name: str
    chrom: str
    strand: str
    intron_start: int
    intron_end: int
    intron_length: int
    three_prime_ss: int
    branchpoint_position: int
    branchpoint_score: float
    branchpoint_to_3ss_nt: int
    branchpoint_type: str
    branchpoint_candidates: int


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("--reference", required=True)
    parser.add_argument(
        "--sample-bam",
        nargs=3,
        action="append",
        metavar=("SAMPLE", "CONDITION", "BAM"),
        required=True,
    )
    parser.add_argument("--downstream-exon-min-offset", type=int, default=5)
    parser.add_argument("--profile-upstream", type=int, required=True)
    parser.add_argument("--profile-downstream", type=int, required=True)
    parser.add_argument("--output-metaprofile-by-sample", required=True)
    parser.add_argument("--output-metaprofile-by-condition", required=True)
    parser.add_argument("--output-summary-by-sample", required=True)
    parser.add_argument("--output-summary-by-condition", required=True)
    args = parser.parse_args()

    if args.downstream_exon_min_offset < 1:
        parser.error("--downstream-exon-min-offset must be positive")
    if args.profile_upstream < 0 or args.profile_downstream < 0:
        parser.error("Profile upstream/downstream values must be non-negative")

    return args


def open_text(path, mode="rt"):
    if path.endswith(".gz"):
        return gzip.open(path, mode)
    return open(path, mode)


def intron_index_bin(position, bin_size=INTRON_INDEX_BIN_SIZE):
    return (position - 1) // bin_size


def load_reference(path):
    introns = {}
    intron_interval_index = defaultdict(list)

    with open_text(path) as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            intron = SelectedIntron(
                intron_id=row["intron_id"],
                transcript_id=row["transcript_id"],
                intron_number=int(row["intron_number"]),
                gene_id=row["gene_id"],
                gene_name=row["gene_name"],
                chrom=row["chrom"],
                strand=row["strand"],
                intron_start=int(row["intron_start"]),
                intron_end=int(row["intron_end"]),
                intron_length=int(row["intron_length"]),
                three_prime_ss=int(row["three_prime_ss"]),
                branchpoint_position=int(row["branchpoint_position"]),
                branchpoint_score=float(row["branchpoint_score"]),
                branchpoint_to_3ss_nt=int(row["branchpoint_to_3ss_nt"]),
                branchpoint_type=row["branchpoint_type"],
                branchpoint_candidates=int(row["branchpoint_candidates"]),
            )
            introns[intron.intron_id] = intron
            start_bin = intron_index_bin(intron.intron_start)
            end_bin = intron_index_bin(intron.intron_end)
            for bin_number in range(start_bin, end_bin + 1):
                intron_interval_index[(intron.chrom, intron.strand, bin_number)].append(intron.intron_id)

    return introns, intron_interval_index


def fragment_strand(alignment):
    return "-" if alignment.is_reverse else "+"


def read1_five_prime_coordinate(alignment):
    if alignment.is_reverse:
        return alignment.reference_end
    return alignment.reference_start + 1


def fragment_three_prime_coordinate(alignment, strand):
    if strand == "+":
        if alignment.template_length == 0:
            return None
        return alignment.reference_start + abs(alignment.template_length)
    if alignment.next_reference_start < 0:
        return None
    return alignment.next_reference_start + 1


def oriented_offset(position, feature_position, strand):
    if strand == "+":
        return position - feature_position
    return feature_position - position


def position_within_intron(position, intron):
    return intron.intron_start <= position <= intron.intron_end


def introns_for_five_prime_position(introns, intron_interval_index, chrom, strand, position):
    intron_ids = intron_interval_index.get((chrom, strand, intron_index_bin(position)), [])
    return [introns[intron_id] for intron_id in intron_ids if position_within_intron(position, introns[intron_id])]


def downstream_exon_spanning_introns(
    introns,
    intron_interval_index,
    chrom,
    strand,
    read1_five_prime,
    fragment_three_prime,
    downstream_exon_min_offset,
):
    hits = []
    for intron in introns_for_five_prime_position(introns, intron_interval_index, chrom, strand, read1_five_prime):
        read1_offset_from_3ss = oriented_offset(read1_five_prime, intron.three_prime_ss, intron.strand)
        fragment_end_offset_from_3ss = oriented_offset(fragment_three_prime, intron.three_prime_ss, intron.strand)
        if read1_offset_from_3ss > 0:
            continue
        if fragment_end_offset_from_3ss < downstream_exon_min_offset:
            continue
        hits.append((intron, read1_offset_from_3ss, fragment_end_offset_from_3ss))
    return hits


def is_unique_primary_fragment(alignment):
    if alignment.is_secondary or alignment.is_supplementary:
        return False
    if alignment.is_unmapped or alignment.mate_is_unmapped:
        return False
    if not alignment.is_proper_pair:
        return False
    if alignment.reference_id != alignment.next_reference_id:
        return False
    if alignment.has_tag("NH"):
        return alignment.get_tag("NH") == 1
    return alignment.mapping_quality > 0


def update_coverage_counts(coverage_counts, start_offset, end_offset, upstream, downstream):
    coverage_start = max(start_offset, -upstream)
    coverage_end = min(end_offset, downstream)
    if coverage_start > coverage_end:
        return 0
    covered_positions = 0
    for offset in range(coverage_start, coverage_end + 1):
        coverage_counts[offset] += 1
        covered_positions += 1
    return covered_positions


def quantify_sample(
    sample,
    condition,
    bam_path,
    introns,
    intron_interval_index,
    downstream_exon_min_offset,
    profile_upstream,
    profile_downstream,
):
    counters = Counter()
    start_counts = Counter()
    coverage_counts = Counter()
    intron_counts = Counter()

    with pysam.AlignmentFile(bam_path, "rb") as bam:
        for alignment in bam.fetch(until_eof=True):
            if not alignment.is_read1:
                continue

            counters["read1_records_examined"] += 1

            if alignment.is_secondary or alignment.is_supplementary:
                counters["filtered_secondary_or_supplementary"] += 1
                continue
            if alignment.is_unmapped or alignment.mate_is_unmapped:
                counters["filtered_unmapped_or_mate_unmapped"] += 1
                continue
            if not alignment.is_proper_pair or alignment.reference_id != alignment.next_reference_id:
                counters["filtered_improper_pair"] += 1
                continue
            if alignment.has_tag("NH") and alignment.get_tag("NH") != 1:
                counters["filtered_nonunique"] += 1
                continue
            if not alignment.has_tag("NH") and alignment.mapping_quality <= 0:
                counters["filtered_nonunique"] += 1
                continue

            strand = fragment_strand(alignment)
            fragment_three_prime = fragment_three_prime_coordinate(alignment, strand)
            if fragment_three_prime is None:
                counters["filtered_missing_fragment_end"] += 1
                continue

            counters["library_fragments"] += 1
            read1_five_prime = read1_five_prime_coordinate(alignment)
            hits = downstream_exon_spanning_introns(
                introns,
                intron_interval_index,
                alignment.reference_name,
                strand,
                read1_five_prime,
                fragment_three_prime,
                downstream_exon_min_offset,
            )
            if not hits:
                counters["fragments_without_downstream_exon_assignment"] += 1
                continue
            if len(hits) > 1:
                counters["ambiguous_downstream_exon_fragments"] += 1
                continue

            intron, _, fragment_end_offset_from_3ss = hits[0]
            start_offset = oriented_offset(read1_five_prime, intron.branchpoint_position, intron.strand)
            end_offset = oriented_offset(fragment_three_prime, intron.branchpoint_position, intron.strand)
            if end_offset < start_offset:
                counters["filtered_inverted_fragment_offsets"] += 1
                continue

            counters["downstream_exon_spanning_fragments"] += 1
            counters["profiled_coverage_positions"] += update_coverage_counts(
                coverage_counts,
                start_offset,
                end_offset,
                profile_upstream,
                profile_downstream,
            )
            intron_counts[intron.intron_id] += 1
            counters["downstream_exon_offset_sum"] += fragment_end_offset_from_3ss

            if -profile_upstream <= start_offset <= profile_downstream:
                counters["profile_window_start_fragments"] += 1
                start_counts[start_offset] += 1
            if start_offset == 0:
                counters["exact_branchpoint_start_fragments"] += 1
            if start_offset == 1:
                counters["plus_one_branchpoint_start_fragments"] += 1
            if start_offset in {0, 1}:
                counters["zero_or_plus_one_branchpoint_start_fragments"] += 1

    selected_fragments = counters["downstream_exon_spanning_fragments"]
    library_fragments = counters["library_fragments"]
    summary_row = {
        "sample": sample,
        "condition": condition,
        "reference_introns": len(introns),
        "downstream_exon_min_offset_nt": downstream_exon_min_offset,
        "profile_upstream_nt": profile_upstream,
        "profile_downstream_nt": profile_downstream,
        "read1_records_examined": counters["read1_records_examined"],
        "library_fragments": library_fragments,
        "downstream_exon_spanning_fragments": selected_fragments,
        "downstream_exon_spanning_fragments_cpm": 0.0
        if library_fragments == 0
        else selected_fragments * 1_000_000.0 / library_fragments,
        "downstream_exon_spanning_introns": len(intron_counts),
        "mean_downstream_exon_end_offset_nt": 0.0
        if selected_fragments == 0
        else counters["downstream_exon_offset_sum"] / selected_fragments,
        "profile_window_start_fragments": counters["profile_window_start_fragments"],
        "profile_window_start_fraction": 0.0
        if selected_fragments == 0
        else counters["profile_window_start_fragments"] / selected_fragments,
        "profile_window_start_percent": 0.0
        if selected_fragments == 0
        else counters["profile_window_start_fragments"] * 100.0 / selected_fragments,
        "profiled_coverage_positions": counters["profiled_coverage_positions"],
        "exact_branchpoint_start_fragments": counters["exact_branchpoint_start_fragments"],
        "exact_branchpoint_start_fraction": 0.0
        if selected_fragments == 0
        else counters["exact_branchpoint_start_fragments"] / selected_fragments,
        "exact_branchpoint_start_percent": 0.0
        if selected_fragments == 0
        else counters["exact_branchpoint_start_fragments"] * 100.0 / selected_fragments,
        "plus_one_branchpoint_start_fragments": counters["plus_one_branchpoint_start_fragments"],
        "plus_one_branchpoint_start_fraction": 0.0
        if selected_fragments == 0
        else counters["plus_one_branchpoint_start_fragments"] / selected_fragments,
        "plus_one_branchpoint_start_percent": 0.0
        if selected_fragments == 0
        else counters["plus_one_branchpoint_start_fragments"] * 100.0 / selected_fragments,
        "zero_or_plus_one_branchpoint_start_fragments": counters["zero_or_plus_one_branchpoint_start_fragments"],
        "zero_or_plus_one_branchpoint_start_fraction": 0.0
        if selected_fragments == 0
        else counters["zero_or_plus_one_branchpoint_start_fragments"] / selected_fragments,
        "zero_or_plus_one_branchpoint_start_percent": 0.0
        if selected_fragments == 0
        else counters["zero_or_plus_one_branchpoint_start_fragments"] * 100.0 / selected_fragments,
        "filtered_secondary_or_supplementary": counters["filtered_secondary_or_supplementary"],
        "filtered_unmapped_or_mate_unmapped": counters["filtered_unmapped_or_mate_unmapped"],
        "filtered_improper_pair": counters["filtered_improper_pair"],
        "filtered_nonunique": counters["filtered_nonunique"],
        "filtered_missing_fragment_end": counters["filtered_missing_fragment_end"],
        "fragments_without_downstream_exon_assignment": counters["fragments_without_downstream_exon_assignment"],
        "ambiguous_downstream_exon_fragments": counters["ambiguous_downstream_exon_fragments"],
        "filtered_inverted_fragment_offsets": counters["filtered_inverted_fragment_offsets"],
    }

    return summary_row, start_counts, coverage_counts


def build_sample_metaprofile_rows(
    sample,
    condition,
    library_fragments,
    selected_fragments,
    start_counts,
    coverage_counts,
    upstream,
    downstream,
):
    rows = []
    for offset in range(-upstream, downstream + 1):
        start_count = start_counts[offset]
        coverage_count = coverage_counts[offset]
        rows.append(
            {
                "sample": sample,
                "condition": condition,
                "offset_nt": offset,
                "read_start_count": start_count,
                "read_start_cpm": 0.0 if library_fragments == 0 else start_count * 1_000_000.0 / library_fragments,
                "read_start_fraction": 0.0 if selected_fragments == 0 else start_count / selected_fragments,
                "read_start_percent": 0.0 if selected_fragments == 0 else start_count * 100.0 / selected_fragments,
                "coverage_count": coverage_count,
                "coverage_cpm": 0.0 if library_fragments == 0 else coverage_count * 1_000_000.0 / library_fragments,
                "coverage_fraction": 0.0 if selected_fragments == 0 else coverage_count / selected_fragments,
                "coverage_percent": 0.0 if selected_fragments == 0 else coverage_count * 100.0 / selected_fragments,
            }
        )
    return rows


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
        for offset in sorted(grouped[condition]):
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
                condition_row[f"ci95_{field}"] = float_ci95_half_width(values)
            condition_rows.append(condition_row)
    return condition_rows


def write_rows(path, rows, fieldnames):
    with open_text(path, "wt") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def condition_order_from_samples(sample_bams):
    order = []
    for _, condition, _ in sample_bams:
        if condition not in order:
            order.append(condition)
    return order


def main():
    args = parse_args()
    introns, intron_interval_index = load_reference(args.reference)

    sample_summary_rows = []
    sample_profile_rows = []
    sample_bams = [(sample, condition, bam) for sample, condition, bam in args.sample_bam]
    condition_order = condition_order_from_samples(sample_bams)

    for sample, condition, bam_path in sample_bams:
        summary_row, start_counts, coverage_counts = quantify_sample(
            sample,
            condition,
            bam_path,
            introns,
            intron_interval_index,
            args.downstream_exon_min_offset,
            args.profile_upstream,
            args.profile_downstream,
        )
        sample_summary_rows.append(summary_row)
        sample_profile_rows.extend(
            build_sample_metaprofile_rows(
                sample,
                condition,
                int(summary_row["library_fragments"]),
                int(summary_row["downstream_exon_spanning_fragments"]),
                start_counts,
                coverage_counts,
                args.profile_upstream,
                args.profile_downstream,
            )
        )
        print(
            f"{sample}: downstream-exon-spanning fragments = "
            f"{summary_row['downstream_exon_spanning_fragments']}"
        )

    sample_summary_rows.sort(key=lambda row: (condition_order.index(row["condition"]), row["sample"]))
    sample_profile_rows.sort(key=lambda row: (condition_order.index(row["condition"]), row["sample"], int(row["offset_nt"])))
    condition_summary_rows = summarise_condition_rows(sample_summary_rows, condition_order)
    condition_profile_rows = summarise_condition_profiles(sample_profile_rows, condition_order)

    write_rows(args.output_summary_by_sample, sample_summary_rows, list(sample_summary_rows[0].keys()))
    write_rows(args.output_summary_by_condition, condition_summary_rows, list(condition_summary_rows[0].keys()))
    write_rows(args.output_metaprofile_by_sample, sample_profile_rows, list(sample_profile_rows[0].keys()))
    write_rows(args.output_metaprofile_by_condition, condition_profile_rows, list(condition_profile_rows[0].keys()))

    print(f"Samples analysed: {len(sample_summary_rows)}")
    print(f"Metaprofile rows written: {len(sample_profile_rows)}")


if __name__ == "__main__":
    main()
