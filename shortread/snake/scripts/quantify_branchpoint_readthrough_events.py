#!/usr/bin/env python3

import argparse
import csv
import gzip
from collections import Counter, defaultdict
from dataclasses import dataclass

import pysam


MATCH_CIGAR_OPS = {0, 7, 8}
DELETION_CIGAR_OP = 2
SKIP_CIGAR_OP = 3
INSERTION_CIGAR_OP = 1


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
    parser.add_argument("--anchor-window", type=int, required=True)
    parser.add_argument("--profile-upstream", type=int, required=True)
    parser.add_argument("--profile-downstream", type=int, required=True)
    parser.add_argument("--output-site-counts", required=True)
    parser.add_argument("--output-position-counts", required=True)
    parser.add_argument("--output-metaprofile", required=True)
    parser.add_argument("--output-summary", required=True)
    return parser.parse_args()


def open_text(path, mode="rt"):
    if path.endswith(".gz"):
        return gzip.open(path, mode)
    return open(path, mode)


def load_reference(path, anchor_window):
    introns = {}
    anchor_index = defaultdict(list)

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
            for position in range(intron.three_prime_ss - anchor_window, intron.three_prime_ss + anchor_window + 1):
                anchor_index[(intron.chrom, intron.strand, position)].append(intron.intron_id)

    return introns, anchor_index


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


def parse_structural_events(alignment, intron, upstream, downstream):
    window_start = intron.branchpoint_position - upstream
    window_end = intron.branchpoint_position + downstream
    branchpoint_zero_based = intron.branchpoint_position - 1

    coverage_offsets = []
    deletion_offsets = []
    insertion_offsets = []
    branchpoint_covered = False

    ref_pos = alignment.reference_start
    alignment_start = alignment.reference_start
    alignment_end = alignment.reference_end

    for op, length in alignment.cigartuples or []:
        if op in MATCH_CIGAR_OPS or op == DELETION_CIGAR_OP:
            segment_start = ref_pos
            segment_end = ref_pos + length
            if segment_start <= branchpoint_zero_based < segment_end:
                branchpoint_covered = True

            overlap_start = max(segment_start, window_start - 1)
            overlap_end = min(segment_end, window_end)
            if overlap_start < overlap_end:
                for zero_based_position in range(overlap_start, overlap_end):
                    reference_position = zero_based_position + 1
                    offset = oriented_offset(reference_position, intron.branchpoint_position, intron.strand)
                    coverage_offsets.append(offset)
                    if op == DELETION_CIGAR_OP:
                        deletion_offsets.append(offset)
            ref_pos = segment_end
        elif op == SKIP_CIGAR_OP:
            ref_pos += length
        elif op == INSERTION_CIGAR_OP:
            anchor_position = None
            if intron.strand == "+":
                if ref_pos > alignment_start:
                    anchor_position = ref_pos
            else:
                if ref_pos < alignment_end:
                    anchor_position = ref_pos + 1
            if anchor_position is not None and window_start <= anchor_position <= window_end:
                insertion_offsets.append(
                    oriented_offset(anchor_position, intron.branchpoint_position, intron.strand)
                )

    return coverage_offsets, deletion_offsets, insertion_offsets, branchpoint_covered


def parse_mismatch_offsets(alignment, intron, upstream, downstream):
    if not alignment.has_tag("MD"):
        return None

    window_start = intron.branchpoint_position - upstream
    window_end = intron.branchpoint_position + downstream
    mismatch_offsets = []

    for query_position, reference_position, reference_base in alignment.get_aligned_pairs(with_seq=True):
        if query_position is None or reference_position is None or reference_base is None:
            continue

        reference_position += 1
        if reference_position < window_start or reference_position > window_end:
            continue
        if reference_base.islower():
            mismatch_offsets.append(
                oriented_offset(reference_position, intron.branchpoint_position, intron.strand)
            )

    return mismatch_offsets


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
        "traversing_fragments",
        "traversing_cpm",
        "traversing_fraction_anchored",
        "traversing_percent_anchored",
        "profiled_coverage_positions",
        "mismatch_events",
        "mismatch_events_per_100_covered_positions",
        "deletion_events",
        "deletion_events_per_100_covered_positions",
        "insertion_events",
        "insertion_events_per_100_covered_positions",
    ]

    ordered_rows = sorted(
        site_counts.items(),
        key=lambda item: (
            -item[1]["traversing_fragments"],
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
            anchored_fragments = counts["anchored_fragments"]
            traversing_fragments = counts["traversing_fragments"]
            profiled_coverage_positions = counts["profiled_coverage_positions"]
            mismatch_events = counts["mismatch_events"]
            deletion_events = counts["deletion_events"]
            insertion_events = counts["insertion_events"]

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
                    "anchored_fragments": anchored_fragments,
                    "traversing_fragments": traversing_fragments,
                    "traversing_cpm": 0
                    if library_size == 0
                    else traversing_fragments * 1_000_000.0 / library_size,
                    "traversing_fraction_anchored": 0
                    if anchored_fragments == 0
                    else traversing_fragments / anchored_fragments,
                    "traversing_percent_anchored": 0
                    if anchored_fragments == 0
                    else traversing_fragments * 100.0 / anchored_fragments,
                    "profiled_coverage_positions": profiled_coverage_positions,
                    "mismatch_events": mismatch_events,
                    "mismatch_events_per_100_covered_positions": 0
                    if profiled_coverage_positions == 0
                    else mismatch_events * 100.0 / profiled_coverage_positions,
                    "deletion_events": deletion_events,
                    "deletion_events_per_100_covered_positions": 0
                    if profiled_coverage_positions == 0
                    else deletion_events * 100.0 / profiled_coverage_positions,
                    "insertion_events": insertion_events,
                    "insertion_events_per_100_covered_positions": 0
                    if profiled_coverage_positions == 0
                    else insertion_events * 100.0 / profiled_coverage_positions,
                }
            )


def write_position_counts(path, sample, condition, position_counts):
    fieldnames = [
        "sample",
        "condition",
        "intron_id",
        "offset_nt",
        "coverage_count",
        "mismatch_count",
        "deletion_count",
        "insertion_count",
    ]
    with open_text(path, "wt") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        for intron_id in sorted(position_counts):
            for offset in sorted(position_counts[intron_id]):
                counts = position_counts[intron_id][offset]
                writer.writerow(
                    {
                        "sample": sample,
                        "condition": condition,
                        "intron_id": intron_id,
                        "offset_nt": offset,
                        "coverage_count": counts["coverage_count"],
                        "mismatch_count": counts["mismatch_count"],
                        "deletion_count": counts["deletion_count"],
                        "insertion_count": counts["insertion_count"],
                    }
                )


def write_metaprofile(
    path,
    sample,
    condition,
    traversing_fragments,
    coverage_counts,
    mismatch_counts,
    deletion_counts,
    insertion_counts,
    upstream,
    downstream,
):
    fieldnames = [
        "sample",
        "condition",
        "offset_nt",
        "coverage_count",
        "coverage_fraction_traversing",
        "coverage_percent_traversing",
        "mismatch_count",
        "mismatch_fraction_coverage",
        "mismatch_percent_coverage",
        "deletion_count",
        "deletion_fraction_coverage",
        "deletion_percent_coverage",
        "insertion_count",
        "insertion_fraction_coverage",
        "insertion_percent_coverage",
    ]

    with open(path, "w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        for offset in range(-upstream, downstream + 1):
            coverage_count = coverage_counts[offset]
            mismatch_count = mismatch_counts[offset]
            deletion_count = deletion_counts[offset]
            insertion_count = insertion_counts[offset]
            writer.writerow(
                {
                    "sample": sample,
                    "condition": condition,
                    "offset_nt": offset,
                    "coverage_count": coverage_count,
                    "coverage_fraction_traversing": 0
                    if traversing_fragments == 0
                    else coverage_count / traversing_fragments,
                    "coverage_percent_traversing": 0
                    if traversing_fragments == 0
                    else coverage_count * 100.0 / traversing_fragments,
                    "mismatch_count": mismatch_count,
                    "mismatch_fraction_coverage": 0
                    if coverage_count == 0
                    else mismatch_count / coverage_count,
                    "mismatch_percent_coverage": 0
                    if coverage_count == 0
                    else mismatch_count * 100.0 / coverage_count,
                    "deletion_count": deletion_count,
                    "deletion_fraction_coverage": 0
                    if coverage_count == 0
                    else deletion_count / coverage_count,
                    "deletion_percent_coverage": 0
                    if coverage_count == 0
                    else deletion_count * 100.0 / coverage_count,
                    "insertion_count": insertion_count,
                    "insertion_fraction_coverage": 0
                    if coverage_count == 0
                    else insertion_count / coverage_count,
                    "insertion_percent_coverage": 0
                    if coverage_count == 0
                    else insertion_count * 100.0 / coverage_count,
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
    introns, anchor_index = load_reference(args.reference, args.anchor_window)

    counters = Counter()
    coverage_counts = Counter()
    mismatch_counts = Counter()
    deletion_counts = Counter()
    insertion_counts = Counter()
    position_counts = defaultdict(
        lambda: defaultdict(
            lambda: {
                "coverage_count": 0,
                "mismatch_count": 0,
                "deletion_count": 0,
                "insertion_count": 0,
            }
        )
    )
    site_counts = defaultdict(
        lambda: {
            "anchored_fragments": 0,
            "traversing_fragments": 0,
            "profiled_coverage_positions": 0,
            "mismatch_events": 0,
            "deletion_events": 0,
            "insertion_events": 0,
        }
    )

    with pysam.AlignmentFile(args.bam, "rb") as bam:
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

            anchor_hits = anchor_index.get((alignment.reference_name, strand, fragment_three_prime), [])
            if not anchor_hits:
                counters["fragments_without_branchpoint_anchor"] += 1
                continue
            if len(anchor_hits) > 1:
                counters["ambiguous_anchor_fragments"] += 1
                continue

            intron = introns[anchor_hits[0]]
            counters["anchored_fragments"] += 1
            site_counts[intron.intron_id]["anchored_fragments"] += 1

            coverage_offsets, deletion_offsets, insertion_offsets, branchpoint_covered = parse_structural_events(
                alignment,
                intron,
                args.profile_upstream,
                args.profile_downstream,
            )

            read1_five_prime = read1_five_prime_coordinate(alignment)
            if not branchpoint_covered or oriented_offset(read1_five_prime, intron.branchpoint_position, intron.strand) >= 0:
                counters["anchored_nontraversing_fragments"] += 1
                continue

            counters["traversing_fragments"] += 1
            site_counts[intron.intron_id]["traversing_fragments"] += 1

            mismatch_offsets = parse_mismatch_offsets(
                alignment,
                intron,
                args.profile_upstream,
                args.profile_downstream,
            )
            if mismatch_offsets is None:
                counters["traversing_fragments_missing_md"] += 1
                mismatch_offsets = []

            counters["profiled_coverage_positions"] += len(coverage_offsets)
            counters["mismatch_events"] += len(mismatch_offsets)
            counters["deletion_events"] += len(deletion_offsets)
            counters["insertion_events"] += len(insertion_offsets)

            site_counts[intron.intron_id]["profiled_coverage_positions"] += len(coverage_offsets)
            site_counts[intron.intron_id]["mismatch_events"] += len(mismatch_offsets)
            site_counts[intron.intron_id]["deletion_events"] += len(deletion_offsets)
            site_counts[intron.intron_id]["insertion_events"] += len(insertion_offsets)

            for offset in coverage_offsets:
                coverage_counts[offset] += 1
                position_counts[intron.intron_id][offset]["coverage_count"] += 1
            for offset in mismatch_offsets:
                mismatch_counts[offset] += 1
                position_counts[intron.intron_id][offset]["mismatch_count"] += 1
            for offset in deletion_offsets:
                deletion_counts[offset] += 1
                position_counts[intron.intron_id][offset]["deletion_count"] += 1
            for offset in insertion_offsets:
                insertion_counts[offset] += 1
                position_counts[intron.intron_id][offset]["insertion_count"] += 1

    summary_row = {
        "sample": args.sample,
        "condition": args.condition,
        "reference_introns": len(introns),
        "anchor_window_nt": args.anchor_window,
        "profile_upstream_nt": args.profile_upstream,
        "profile_downstream_nt": args.profile_downstream,
        "read1_records_examined": counters["read1_records_examined"],
        "library_fragments": counters["library_fragments"],
        "anchored_fragments": counters["anchored_fragments"],
        "anchored_fragments_cpm": 0
        if counters["library_fragments"] == 0
        else counters["anchored_fragments"] * 1_000_000.0 / counters["library_fragments"],
        "traversing_fragments": counters["traversing_fragments"],
        "traversing_fragments_cpm": 0
        if counters["library_fragments"] == 0
        else counters["traversing_fragments"] * 1_000_000.0 / counters["library_fragments"],
        "traversing_introns_with_reads": sum(1 for counts in site_counts.values() if counts["traversing_fragments"] > 0),
        "traversing_fraction_anchored": 0
        if counters["anchored_fragments"] == 0
        else counters["traversing_fragments"] / counters["anchored_fragments"],
        "traversing_percent_anchored": 0
        if counters["anchored_fragments"] == 0
        else counters["traversing_fragments"] * 100.0 / counters["anchored_fragments"],
        "profiled_coverage_positions": counters["profiled_coverage_positions"],
        "mismatch_events": counters["mismatch_events"],
        "mismatch_events_per_100_covered_positions": 0
        if counters["profiled_coverage_positions"] == 0
        else counters["mismatch_events"] * 100.0 / counters["profiled_coverage_positions"],
        "deletion_events": counters["deletion_events"],
        "deletion_events_per_100_covered_positions": 0
        if counters["profiled_coverage_positions"] == 0
        else counters["deletion_events"] * 100.0 / counters["profiled_coverage_positions"],
        "insertion_events": counters["insertion_events"],
        "insertion_events_per_100_covered_positions": 0
        if counters["profiled_coverage_positions"] == 0
        else counters["insertion_events"] * 100.0 / counters["profiled_coverage_positions"],
        "filtered_secondary_or_supplementary": counters["filtered_secondary_or_supplementary"],
        "filtered_unmapped_or_mate_unmapped": counters["filtered_unmapped_or_mate_unmapped"],
        "filtered_improper_pair": counters["filtered_improper_pair"],
        "filtered_nonunique": counters["filtered_nonunique"],
        "filtered_missing_fragment_end": counters["filtered_missing_fragment_end"],
        "fragments_without_branchpoint_anchor": counters["fragments_without_branchpoint_anchor"],
        "ambiguous_anchor_fragments": counters["ambiguous_anchor_fragments"],
        "anchored_nontraversing_fragments": counters["anchored_nontraversing_fragments"],
        "traversing_fragments_missing_md": counters["traversing_fragments_missing_md"],
    }

    write_site_counts(
        args.output_site_counts,
        args.sample,
        args.condition,
        introns,
        site_counts,
        counters["library_fragments"],
    )
    write_position_counts(args.output_position_counts, args.sample, args.condition, position_counts)
    write_metaprofile(
        args.output_metaprofile,
        args.sample,
        args.condition,
        counters["traversing_fragments"],
        coverage_counts,
        mismatch_counts,
        deletion_counts,
        insertion_counts,
        args.profile_upstream,
        args.profile_downstream,
    )
    write_summary(args.output_summary, summary_row)

    print(f"Sample: {args.sample}")
    print(f"Condition: {args.condition}")
    print(f"Library fragments: {counters['library_fragments']}")
    print(f"Anchored fragments: {counters['anchored_fragments']}")
    print(f"Traversing fragments: {counters['traversing_fragments']}")
    print(f"Profiled coverage positions: {counters['profiled_coverage_positions']}")
    print(f"Mismatch events: {counters['mismatch_events']}")
    print(f"Deletion events: {counters['deletion_events']}")
    print(f"Insertion events: {counters['insertion_events']}")


if __name__ == "__main__":
    main()
