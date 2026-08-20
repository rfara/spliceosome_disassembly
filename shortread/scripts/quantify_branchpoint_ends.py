#!/usr/bin/env python3

import argparse
import csv
import gzip
import subprocess
from collections import Counter, defaultdict
from contextlib import contextmanager
from dataclasses import dataclass

import pysam


INTRON_INDEX_BIN_SIZE = 1000
CIGAR_REF_SKIP = 3


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
    parser.add_argument("--bam", required=True)
    parser.add_argument("--reference", required=True)
    parser.add_argument("--sample", required=True)
    parser.add_argument("--condition", required=True)
    parser.add_argument("--anchor-window", type=int)
    parser.add_argument("--anchor-upstream", type=int)
    parser.add_argument("--anchor-downstream", type=int)
    parser.add_argument("--profile-upstream", type=int, required=True)
    parser.add_argument("--profile-downstream", type=int, required=True)
    parser.add_argument("--three-prime-coverage-upstream", type=int, required=True)
    parser.add_argument("--three-prime-coverage-downstream", type=int, required=True)
    parser.add_argument("--output-site-counts", required=True)
    parser.add_argument("--output-offset-counts", required=True)
    parser.add_argument("--output-three-prime-coverage", required=True)
    parser.add_argument("--output-metaprofile", required=True)
    parser.add_argument("--output-summary", required=True)
    args = parser.parse_args()

    if args.anchor_window is None and (args.anchor_upstream is None or args.anchor_downstream is None):
        parser.error("Provide --anchor-window or both --anchor-upstream and --anchor-downstream")

    args.anchor_upstream = args.anchor_window if args.anchor_upstream is None else args.anchor_upstream
    args.anchor_downstream = args.anchor_window if args.anchor_downstream is None else args.anchor_downstream

    if args.anchor_upstream < 0 or args.anchor_downstream < 0:
        parser.error("Anchor upstream/downstream limits must be non-negative")

    return args


def open_text(path, mode="rt"):
    if path.endswith(".gz"):
        return gzip.open(path, mode)
    return open(path, mode)


def intron_index_bin(position, bin_size=INTRON_INDEX_BIN_SIZE):
    return (position - 1) // bin_size


def oriented_position_from_offset(feature_position, offset, strand):
    if strand == "+":
        return feature_position + offset
    return feature_position - offset


def load_reference(path, anchor_upstream, anchor_downstream):
    introns = {}
    anchor_index = defaultdict(list)
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
            for offset in range(-anchor_upstream, anchor_downstream + 1):
                position = oriented_position_from_offset(intron.three_prime_ss, offset, intron.strand)
                anchor_index[(intron.chrom, intron.strand, position)].append(intron.intron_id)
            start_bin = intron_index_bin(intron.intron_start)
            end_bin = intron_index_bin(intron.intron_end)
            for bin_number in range(start_bin, end_bin + 1):
                intron_interval_index[(intron.chrom, intron.strand, bin_number)].append(intron.intron_id)

    return introns, anchor_index, intron_interval_index


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


def is_spliced_alignment(alignment):
    if alignment.cigartuples is None:
        return False
    return any(operation == CIGAR_REF_SKIP and length > 0 for operation, length in alignment.cigartuples)


@contextmanager
def open_queryname_collated_bam(path):
    process = subprocess.Popen(
        ["samtools", "collate", "-u", "-O", path],
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
    )
    if process.stdout is None:
        raise RuntimeError(f"Failed to open samtools collate stdout for {path}")

    bam = pysam.AlignmentFile(process.stdout, "rb")
    try:
        yield bam
    finally:
        bam.close()
        stderr = b""
        if process.stderr is not None:
            stderr = process.stderr.read()
            process.stderr.close()
        return_code = process.wait()
        if return_code != 0:
            raise RuntimeError(
                f"samtools collate failed for {path}: {stderr.decode(errors='replace').strip()}"
            )


def iter_queryname_groups(bam):
    current_name = None
    current_group = []
    for alignment in bam.fetch(until_eof=True):
        if current_name is None or alignment.query_name == current_name:
            current_group.append(alignment)
            current_name = alignment.query_name
            continue
        yield current_group
        current_group = [alignment]
        current_name = alignment.query_name
    if current_group:
        yield current_group


def primary_alignment(group, read1):
    primary = None
    for alignment in group:
        if alignment.is_read1 != read1:
            continue
        if alignment.is_secondary or alignment.is_supplementary:
            continue
        if primary is not None:
            return None
        primary = alignment
    return primary


def position_within_intron(position, intron):
    return intron.intron_start <= position <= intron.intron_end


def introns_for_five_prime_position(introns, intron_interval_index, chrom, strand, position):
    intron_ids = intron_interval_index.get((chrom, strand, intron_index_bin(position)), [])
    return [introns[intron_id] for intron_id in intron_ids if position_within_intron(position, introns[intron_id])]


def three_prime_spanning_introns(
    introns,
    intron_interval_index,
    chrom,
    strand,
    read1_five_prime,
    fragment_three_prime,
):
    hits = []
    for intron in introns_for_five_prime_position(introns, intron_interval_index, chrom, strand, read1_five_prime):
        read1_five_prime_offset = oriented_offset(read1_five_prime, intron.three_prime_ss, intron.strand)
        fragment_three_prime_offset = oriented_offset(fragment_three_prime, intron.three_prime_ss, intron.strand)
        if read1_five_prime_offset > 0:
            continue
        if fragment_three_prime_offset < 0:
            continue
        hits.append((intron, read1_five_prime_offset, fragment_three_prime_offset))
    return hits


def update_three_prime_coverage_counts(
    intron_coverage_counts,
    intron_id,
    start_offset,
    end_offset,
    upstream,
    downstream,
):
    coverage_start = max(start_offset, -upstream)
    coverage_end = min(end_offset, downstream)
    if coverage_start > coverage_end:
        return
    for offset in range(coverage_start, coverage_end + 1):
        intron_coverage_counts[intron_id][offset] += 1


def write_site_counts(path, sample, condition, introns, site_counts, library_size):
    fieldnames = [
        "sample",
        "condition",
        "intron_id",
        "gene_id",
        "gene_name",
        "transcript_id",
        "intron_number",
        "chrom",
        "strand",
        "intron_start",
        "intron_end",
        "intron_length",
        "three_prime_ss",
        "branchpoint_position",
        "branchpoint_score",
        "branchpoint_to_3ss_nt",
        "branchpoint_type",
        "branchpoint_candidates",
        "anchored_fragments",
        "exact_branchpoint_fragments",
        "exact_branchpoint_cpm",
        "exact_branchpoint_fraction_anchored",
        "exact_branchpoint_percent_anchored",
        "plus_one_branchpoint_fragments",
        "plus_one_branchpoint_cpm",
        "plus_one_branchpoint_fraction_anchored",
        "plus_one_branchpoint_percent_anchored",
        "zero_or_plus_one_branchpoint_fragments",
        "zero_or_plus_one_branchpoint_cpm",
        "zero_or_plus_one_branchpoint_fraction_anchored",
        "zero_or_plus_one_branchpoint_percent_anchored",
    ]

    ordered_rows = sorted(
        site_counts.items(),
        key=lambda item: (
            -item[1]["exact_branchpoint_fragments"],
            -item[1]["anchored_fragments"],
            introns[item[0]].gene_name,
            item[0],
        ),
    )

    with open_text(path, "wt") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        for intron_id, counts in ordered_rows:
            intron = introns[intron_id]
            anchored = counts["anchored_fragments"]
            exact = counts["exact_branchpoint_fragments"]
            plus_one = counts["plus_one_branchpoint_fragments"]
            zero_or_plus_one = exact + plus_one
            exact_fraction = 0 if anchored == 0 else (exact / anchored)
            plus_one_fraction = 0 if anchored == 0 else (plus_one / anchored)
            zero_or_plus_one_fraction = 0 if anchored == 0 else (zero_or_plus_one / anchored)
            writer.writerow(
                {
                    "sample": sample,
                    "condition": condition,
                    "intron_id": intron.intron_id,
                    "gene_id": intron.gene_id,
                    "gene_name": intron.gene_name,
                    "transcript_id": intron.transcript_id,
                    "intron_number": intron.intron_number,
                    "chrom": intron.chrom,
                    "strand": intron.strand,
                    "intron_start": intron.intron_start,
                    "intron_end": intron.intron_end,
                    "intron_length": intron.intron_length,
                    "three_prime_ss": intron.three_prime_ss,
                    "branchpoint_position": intron.branchpoint_position,
                    "branchpoint_score": intron.branchpoint_score,
                    "branchpoint_to_3ss_nt": intron.branchpoint_to_3ss_nt,
                    "branchpoint_type": intron.branchpoint_type,
                    "branchpoint_candidates": intron.branchpoint_candidates,
                    "anchored_fragments": anchored,
                    "exact_branchpoint_fragments": exact,
                    "exact_branchpoint_cpm": 0 if library_size == 0 else (exact * 1_000_000.0 / library_size),
                    "exact_branchpoint_fraction_anchored": exact_fraction,
                    "exact_branchpoint_percent_anchored": exact_fraction * 100.0,
                    "plus_one_branchpoint_fragments": plus_one,
                    "plus_one_branchpoint_cpm": 0 if library_size == 0 else (plus_one * 1_000_000.0 / library_size),
                    "plus_one_branchpoint_fraction_anchored": plus_one_fraction,
                    "plus_one_branchpoint_percent_anchored": plus_one_fraction * 100.0,
                    "zero_or_plus_one_branchpoint_fragments": zero_or_plus_one,
                    "zero_or_plus_one_branchpoint_cpm": 0
                    if library_size == 0
                    else (zero_or_plus_one * 1_000_000.0 / library_size),
                    "zero_or_plus_one_branchpoint_fraction_anchored": zero_or_plus_one_fraction,
                    "zero_or_plus_one_branchpoint_percent_anchored": zero_or_plus_one_fraction * 100.0,
                }
            )


def write_metaprofile(path, sample, condition, library_size, anchored_fragments, profile_counts, upstream, downstream):
    fieldnames = [
        "sample",
        "condition",
        "offset_nt",
        "read_count",
        "cpm",
        "anchored_fraction",
        "anchored_percent",
    ]
    with open(path, "w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        for offset in range(-upstream, downstream + 1):
            read_count = profile_counts[offset]
            writer.writerow(
                {
                    "sample": sample,
                    "condition": condition,
                    "offset_nt": offset,
                    "read_count": read_count,
                    "cpm": 0 if library_size == 0 else (read_count * 1_000_000.0 / library_size),
                    "anchored_fraction": 0 if anchored_fragments == 0 else (read_count / anchored_fragments),
                    "anchored_percent": 0 if anchored_fragments == 0 else (read_count / anchored_fragments) * 100.0,
                }
            )


def write_offset_counts(path, sample, condition, offset_counts):
    fieldnames = ["sample", "condition", "intron_id", "offset_nt", "read_count"]
    with open_text(path, "wt") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        for intron_id in sorted(offset_counts):
            for offset in sorted(offset_counts[intron_id]):
                writer.writerow(
                    {
                        "sample": sample,
                        "condition": condition,
                        "intron_id": intron_id,
                        "offset_nt": offset,
                        "read_count": offset_counts[intron_id][offset],
                    }
                )


def write_three_prime_coverage(path, sample, condition, coverage_counts, spanning_counts):
    fieldnames = [
        "sample",
        "condition",
        "intron_id",
        "three_prime_spanning_fragments",
        "offset_nt",
        "coverage_count",
    ]
    with open_text(path, "wt") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        for intron_id in sorted(coverage_counts):
            for offset in sorted(coverage_counts[intron_id]):
                writer.writerow(
                    {
                        "sample": sample,
                        "condition": condition,
                        "intron_id": intron_id,
                        "three_prime_spanning_fragments": spanning_counts[intron_id],
                        "offset_nt": offset,
                        "coverage_count": coverage_counts[intron_id][offset],
                    }
                )


def write_summary(path, summary_row):
    fieldnames = list(summary_row.keys())
    with open(path, "w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerow(summary_row)


def main():
    args = parse_args()
    introns, anchor_index, intron_interval_index = load_reference(
        args.reference,
        args.anchor_upstream,
        args.anchor_downstream,
    )

    counters = Counter()
    profile_counts = Counter()
    offset_counts = defaultdict(Counter)
    intron_three_prime_coverage_counts = defaultdict(Counter)
    intron_three_prime_spanning_counts = Counter()
    site_counts = defaultdict(
        lambda: {
            "anchored_fragments": 0,
            "exact_branchpoint_fragments": 0,
            "plus_one_branchpoint_fragments": 0,
        }
    )

    with open_queryname_collated_bam(args.bam) as bam:
        for group in iter_queryname_groups(bam):
            for record in group:
                if not record.is_read1:
                    continue
                counters["read1_records_examined"] += 1
                if record.is_secondary or record.is_supplementary:
                    counters["filtered_secondary_or_supplementary"] += 1

            alignment = primary_alignment(group, True)
            if alignment is None:
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
            mate_alignment = primary_alignment(group, False)
            if mate_alignment is None:
                counters["filtered_missing_mate_alignment"] += 1
                continue
            read1_five_prime = read1_five_prime_coordinate(alignment)

            # These metaprofiles use continuous fragment spans, so drop any fragment whose
            # read1 or mate alignment crosses a splice junction.
            if is_spliced_alignment(alignment) or is_spliced_alignment(mate_alignment):
                counters["filtered_spliced_fragments"] += 1
                continue

            three_prime_hits = three_prime_spanning_introns(
                introns,
                intron_interval_index,
                alignment.reference_name,
                strand,
                read1_five_prime,
                fragment_three_prime,
            )
            if not three_prime_hits:
                counters["fragments_without_three_prime_coverage_assignment"] += 1
            elif len(three_prime_hits) > 1:
                counters["ambiguous_three_prime_coverage_fragments"] += 1
            else:
                coverage_intron, coverage_start_offset, coverage_end_offset = three_prime_hits[0]
                counters["three_prime_spanning_fragments"] += 1
                intron_three_prime_spanning_counts[coverage_intron.intron_id] += 1
                update_three_prime_coverage_counts(
                    intron_three_prime_coverage_counts,
                    coverage_intron.intron_id,
                    coverage_start_offset,
                    coverage_end_offset,
                    args.three_prime_coverage_upstream,
                    args.three_prime_coverage_downstream,
                )

            anchor_hits = anchor_index.get((alignment.reference_name, strand, fragment_three_prime), [])
            if not anchor_hits:
                counters["fragments_without_branchpoint_anchor"] += 1
                continue
            if len(anchor_hits) > 1:
                counters["ambiguous_anchor_fragments"] += 1
                continue

            intron = introns[anchor_hits[0]]
            if not position_within_intron(read1_five_prime, intron):
                counters["filtered_five_prime_outside_intron"] += 1
                continue
            offset = oriented_offset(read1_five_prime, intron.branchpoint_position, intron.strand)

            counters["anchored_fragments"] += 1
            site_counts[intron.intron_id]["anchored_fragments"] += 1

            if offset == 0:
                counters["exact_branchpoint_fragments"] += 1
                site_counts[intron.intron_id]["exact_branchpoint_fragments"] += 1
            if offset == 1:
                counters["plus_one_branchpoint_fragments"] += 1
                site_counts[intron.intron_id]["plus_one_branchpoint_fragments"] += 1
            if offset in {0, 1}:
                counters["zero_or_plus_one_branchpoint_fragments"] += 1

            offset_counts[intron.intron_id][offset] += 1
            if -args.profile_upstream <= offset <= args.profile_downstream:
                counters["profile_window_fragments"] += 1
                profile_counts[offset] += 1

    summary_row = {
        "sample": args.sample,
        "condition": args.condition,
        "reference_introns": len(introns),
        "anchor_window_nt": max(args.anchor_upstream, args.anchor_downstream),
        "anchor_upstream_nt": args.anchor_upstream,
        "anchor_downstream_nt": args.anchor_downstream,
        "profile_upstream_nt": args.profile_upstream,
        "profile_downstream_nt": args.profile_downstream,
        "three_prime_coverage_upstream_nt": args.three_prime_coverage_upstream,
        "three_prime_coverage_downstream_nt": args.three_prime_coverage_downstream,
        "read1_records_examined": counters["read1_records_examined"],
        "library_fragments": counters["library_fragments"],
        "three_prime_spanning_fragments": counters["three_prime_spanning_fragments"],
        "three_prime_spanning_fragments_cpm": 0
        if counters["library_fragments"] == 0
        else counters["three_prime_spanning_fragments"] * 1_000_000.0 / counters["library_fragments"],
        "three_prime_introns_with_coverage": len(intron_three_prime_coverage_counts),
        "anchored_fragments": counters["anchored_fragments"],
        "anchored_fragments_cpm": 0
        if counters["library_fragments"] == 0
        else counters["anchored_fragments"] * 1_000_000.0 / counters["library_fragments"],
        "anchored_introns_with_reads": len(site_counts),
        "exact_branchpoint_fragments": counters["exact_branchpoint_fragments"],
        "exact_branchpoint_cpm": 0
        if counters["library_fragments"] == 0
        else counters["exact_branchpoint_fragments"] * 1_000_000.0 / counters["library_fragments"],
        "exact_branchpoint_fraction_anchored": 0
        if counters["anchored_fragments"] == 0
        else counters["exact_branchpoint_fragments"] / counters["anchored_fragments"],
        "exact_branchpoint_percent_anchored": 0
        if counters["anchored_fragments"] == 0
        else (counters["exact_branchpoint_fragments"] / counters["anchored_fragments"]) * 100.0,
        "plus_one_branchpoint_fragments": counters["plus_one_branchpoint_fragments"],
        "plus_one_branchpoint_cpm": 0
        if counters["library_fragments"] == 0
        else counters["plus_one_branchpoint_fragments"] * 1_000_000.0 / counters["library_fragments"],
        "plus_one_branchpoint_fraction_anchored": 0
        if counters["anchored_fragments"] == 0
        else counters["plus_one_branchpoint_fragments"] / counters["anchored_fragments"],
        "plus_one_branchpoint_percent_anchored": 0
        if counters["anchored_fragments"] == 0
        else (counters["plus_one_branchpoint_fragments"] / counters["anchored_fragments"]) * 100.0,
        "zero_or_plus_one_branchpoint_fragments": counters["zero_or_plus_one_branchpoint_fragments"],
        "zero_or_plus_one_branchpoint_cpm": 0
        if counters["library_fragments"] == 0
        else counters["zero_or_plus_one_branchpoint_fragments"] * 1_000_000.0 / counters["library_fragments"],
        "zero_or_plus_one_branchpoint_fraction_anchored": 0
        if counters["anchored_fragments"] == 0
        else counters["zero_or_plus_one_branchpoint_fragments"] / counters["anchored_fragments"],
        "zero_or_plus_one_branchpoint_percent_anchored": 0
        if counters["anchored_fragments"] == 0
        else (counters["zero_or_plus_one_branchpoint_fragments"] / counters["anchored_fragments"]) * 100.0,
        "profile_window_fragments": counters["profile_window_fragments"],
        "profile_window_fraction_anchored": 0
        if counters["anchored_fragments"] == 0
        else counters["profile_window_fragments"] / counters["anchored_fragments"],
        "profile_window_percent_anchored": 0
        if counters["anchored_fragments"] == 0
        else (counters["profile_window_fragments"] / counters["anchored_fragments"]) * 100.0,
        "filtered_secondary_or_supplementary": counters["filtered_secondary_or_supplementary"],
        "filtered_unmapped_or_mate_unmapped": counters["filtered_unmapped_or_mate_unmapped"],
        "filtered_improper_pair": counters["filtered_improper_pair"],
        "filtered_nonunique": counters["filtered_nonunique"],
        "filtered_missing_fragment_end": counters["filtered_missing_fragment_end"],
        "filtered_missing_mate_alignment": counters["filtered_missing_mate_alignment"],
        "filtered_spliced_fragments": counters["filtered_spliced_fragments"],
        "fragments_without_three_prime_coverage_assignment": counters[
            "fragments_without_three_prime_coverage_assignment"
        ],
        "ambiguous_three_prime_coverage_fragments": counters["ambiguous_three_prime_coverage_fragments"],
        "fragments_without_branchpoint_anchor": counters["fragments_without_branchpoint_anchor"],
        "ambiguous_anchor_fragments": counters["ambiguous_anchor_fragments"],
        "filtered_five_prime_outside_intron": counters["filtered_five_prime_outside_intron"],
    }

    write_site_counts(
        args.output_site_counts,
        args.sample,
        args.condition,
        introns,
        site_counts,
        counters["library_fragments"],
    )
    write_offset_counts(args.output_offset_counts, args.sample, args.condition, offset_counts)
    write_three_prime_coverage(
        args.output_three_prime_coverage,
        args.sample,
        args.condition,
        intron_three_prime_coverage_counts,
        intron_three_prime_spanning_counts,
    )
    write_metaprofile(
        args.output_metaprofile,
        args.sample,
        args.condition,
        counters["library_fragments"],
        counters["anchored_fragments"],
        profile_counts,
        args.profile_upstream,
        args.profile_downstream,
    )
    write_summary(args.output_summary, summary_row)

    print(f"Sample: {args.sample}")
    print(f"Condition: {args.condition}")
    print(f"Anchor window: -{args.anchor_upstream}..+{args.anchor_downstream} nt from 3'SS")
    print(f"Library fragments: {counters['library_fragments']}")
    print(f"Filtered spliced fragments: {counters['filtered_spliced_fragments']}")
    print(f"3'SS-spanning fragments: {counters['three_prime_spanning_fragments']}")
    print(f"Anchored fragments: {counters['anchored_fragments']}")
    print(f"Exact branchpoint fragments: {counters['exact_branchpoint_fragments']}")
    print(f"Plus-one branchpoint fragments: {counters['plus_one_branchpoint_fragments']}")
    print(f"Zero-or-plus-one branchpoint fragments: {counters['zero_or_plus_one_branchpoint_fragments']}")


if __name__ == "__main__":
    main()
