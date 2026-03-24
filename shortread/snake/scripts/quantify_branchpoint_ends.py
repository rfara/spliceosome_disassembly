#!/usr/bin/env python3

import argparse
import csv
import gzip
from collections import Counter, defaultdict
from dataclasses import dataclass

import pysam


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
    parser.add_argument("--output-offset-counts", required=True)
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
    profile_counts = Counter()
    offset_counts = defaultdict(Counter)
    site_counts = defaultdict(
        lambda: {
            "anchored_fragments": 0,
            "exact_branchpoint_fragments": 0,
            "plus_one_branchpoint_fragments": 0,
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
            read1_five_prime = read1_five_prime_coordinate(alignment)
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

            if -args.profile_upstream <= offset <= args.profile_downstream:
                counters["profile_window_fragments"] += 1
                profile_counts[offset] += 1
                offset_counts[intron.intron_id][offset] += 1

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
        "fragments_without_branchpoint_anchor": counters["fragments_without_branchpoint_anchor"],
        "ambiguous_anchor_fragments": counters["ambiguous_anchor_fragments"],
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
    print(f"Library fragments: {counters['library_fragments']}")
    print(f"Anchored fragments: {counters['anchored_fragments']}")
    print(f"Exact branchpoint fragments: {counters['exact_branchpoint_fragments']}")
    print(f"Plus-one branchpoint fragments: {counters['plus_one_branchpoint_fragments']}")
    print(f"Zero-or-plus-one branchpoint fragments: {counters['zero_or_plus_one_branchpoint_fragments']}")


if __name__ == "__main__":
    main()
