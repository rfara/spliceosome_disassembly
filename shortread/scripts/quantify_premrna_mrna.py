#!/usr/bin/env python3

import argparse
import csv
import gzip
import pickle
import signal
import subprocess
from collections import Counter, defaultdict

import pysam
from intervaltree import IntervalTree


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("--bam", required=True)
    parser.add_argument("--reference", required=True)
    parser.add_argument("--sample", required=True)
    parser.add_argument("--condition", required=True)
    parser.add_argument("--threads", type=int, default=1)
    parser.add_argument("--metaprofile-upstream", type=int, required=True)
    parser.add_argument("--metaprofile-downstream", type=int, required=True)
    parser.add_argument("--output-gene-counts", required=True)
    parser.add_argument("--output-coverage", required=True)
    parser.add_argument("--output-summary", required=True)
    return parser.parse_args()


def open_text(path, mode="rt"):
    if path.endswith(".gz"):
        return gzip.open(path, mode)
    return open(path, mode)


def fragment_strand(read1_alignment):
    return "-" if read1_alignment.is_reverse else "+"


def iter_name_collated_groups(bam_path, threads):
    command = ["samtools", "collate", "-f", "-u", "-O", "-@", str(max(1, threads)), bam_path]
    process = subprocess.Popen(command, stdout=subprocess.PIPE)
    try:
        with pysam.AlignmentFile(process.stdout, "rb") as bam:
            current_name = None
            bucket = []
            for alignment in bam.fetch(until_eof=True):
                if current_name is None:
                    current_name = alignment.query_name
                if alignment.query_name != current_name:
                    yield current_name, bucket
                    current_name = alignment.query_name
                    bucket = [alignment]
                else:
                    bucket.append(alignment)
            if bucket:
                yield current_name, bucket
    finally:
        if process.stdout is not None:
            process.stdout.close()
        return_code = process.wait()
        if return_code not in (0, -signal.SIGPIPE):
            raise subprocess.CalledProcessError(return_code, command)


def window_exon_spans(gene_record, upstream, downstream):
    anchor = gene_record["three_prime_end"]
    if gene_record["strand"] == "+":
        window_start = anchor - upstream
        window_end = anchor + downstream
    else:
        window_start = anchor - downstream
        window_end = anchor + upstream

    spans = []
    for exon_start, exon_end in gene_record["exon_spans"]:
        overlap_start = max(exon_start, window_start)
        overlap_end = min(exon_end, window_end)
        if overlap_start <= overlap_end:
            spans.append((overlap_start, overlap_end))
    return tuple(spans)


def load_reference(path, upstream, downstream):
    with open(path, "rb") as handle:
        reference = pickle.load(handle)

    genes = {}
    exon_trees = defaultdict(IntervalTree)
    intron_trees = defaultdict(IntervalTree)
    small_ncrna_trees = defaultdict(IntervalTree)

    for gene_id, raw_gene in reference["genes"].items():
        gene = dict(raw_gene)
        gene["exon_spans"] = tuple(tuple(span) for span in raw_gene["exon_spans"])
        gene["intron_spans"] = tuple(tuple(span) for span in raw_gene["intron_spans"])
        gene["window_exon_spans"] = window_exon_spans(gene, upstream, downstream)
        genes[gene_id] = gene

        exon_tree = exon_trees[(gene["chrom"], gene["strand"])]
        intron_tree = intron_trees[(gene["chrom"], gene["strand"])]
        for exon_start, exon_end in gene["exon_spans"]:
            exon_tree.addi(exon_start - 1, exon_end, gene_id)
        for intron_start, intron_end in gene["intron_spans"]:
            intron_tree.addi(intron_start - 1, intron_end, gene_id)

    for raw_gene in reference["small_ncrna_genes"]:
        small_ncrna_trees[(raw_gene["chrom"], raw_gene["strand"])].addi(
            raw_gene["start"] - 1,
            raw_gene["end"],
            raw_gene["gene_id"],
        )

    return genes, exon_trees, intron_trees, small_ncrna_trees, len(reference["small_ncrna_genes"])


def overlap_ids(tree_map, chrom, strand, start0, end0):
    tree = tree_map.get((chrom, strand))
    if tree is None:
        return set()
    return {interval.data for interval in tree.overlap(start0, end0)}


def has_overlap(tree_map, chrom, strand, start0, end0):
    tree = tree_map.get((chrom, strand))
    if tree is None:
        return False
    return bool(tree.overlap(start0, end0))


def classify_fragment(primary_mapped_alignments, read1_alignment, exon_trees, intron_trees, small_ncrna_trees):
    strand = fragment_strand(read1_alignment)
    exon_gene_ids = set()
    intron_gene_ids = set()
    small_ncrna_overlap = False

    for alignment in primary_mapped_alignments:
        chrom = alignment.reference_name
        for block_start0, block_end0 in alignment.get_blocks():
            exon_gene_ids.update(overlap_ids(exon_trees, chrom, strand, block_start0, block_end0))
            intron_gene_ids.update(overlap_ids(intron_trees, chrom, strand, block_start0, block_end0))
            if not small_ncrna_overlap and has_overlap(small_ncrna_trees, chrom, strand, block_start0, block_end0):
                small_ncrna_overlap = True

    candidate_gene_ids = exon_gene_ids.union(intron_gene_ids)
    if not candidate_gene_ids:
        return "no_gene_overlap", None
    if len(candidate_gene_ids) > 1:
        return "ambiguous_same_strand_genes", None

    gene_id = next(iter(candidate_gene_ids))
    if small_ncrna_overlap:
        return "small_ncrna_overlap", gene_id
    if gene_id in intron_gene_ids:
        return "intronic", gene_id
    if gene_id in exon_gene_ids:
        return "mrna", gene_id
    return "unclassified", gene_id


def interval_offsets(interval_start, interval_end, anchor, strand):
    if strand == "+":
        return interval_start - anchor, interval_end - anchor
    return anchor - interval_end, anchor - interval_start


def update_coverage_deltas(coverage_deltas, gene_record, alignments, upstream, downstream):
    window_size = upstream + downstream + 1
    delta = coverage_deltas.setdefault(gene_record["gene_id"], [0] * (window_size + 1))
    anchor = gene_record["three_prime_end"]

    for alignment in alignments:
        for block_start0, block_end0 in alignment.get_blocks():
            block_start = block_start0 + 1
            block_end = block_end0
            for exon_start, exon_end in gene_record["window_exon_spans"]:
                overlap_start = max(block_start, exon_start)
                overlap_end = min(block_end, exon_end)
                if overlap_start > overlap_end:
                    continue

                start_offset, end_offset = interval_offsets(
                    overlap_start,
                    overlap_end,
                    anchor,
                    gene_record["strand"],
                )
                start_offset = max(start_offset, -upstream)
                end_offset = min(end_offset, downstream)
                if start_offset > end_offset:
                    continue

                start_index = start_offset + upstream
                end_index = end_offset + upstream
                delta[start_index] += 1
                delta[end_index + 1] -= 1


def ratio_or_na(numerator, denominator):
    if denominator == 0:
        return "NA"
    return numerator / denominator


def write_gene_counts(path, sample, condition, genes, gene_counts, library_fragments):
    fieldnames = [
        "sample",
        "condition",
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
        "transcript_length",
        "exonic_length",
        "exon_count",
        "intron_count",
        "total_gene_fragments",
        "total_gene_cpm",
        "mrna_fragments",
        "mrna_cpm",
        "mrna_fraction_gene_reads",
        "mrna_percent_gene_reads",
        "intronic_fragments",
        "intronic_cpm",
        "intronic_fraction_gene_reads",
        "intronic_percent_gene_reads",
        "intronic_to_mrna_ratio",
    ]

    ordered_gene_ids = sorted(
        gene_counts,
        key=lambda gene_id: (
            -sum(gene_counts[gene_id].values()),
            genes[gene_id]["gene_name"],
            gene_id,
        ),
    )

    with open_text(path, "wt") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        for gene_id in ordered_gene_ids:
            gene = genes[gene_id]
            mrna_fragments = gene_counts[gene_id]["mrna_fragments"]
            intronic_fragments = gene_counts[gene_id]["intronic_fragments"]
            total_gene_fragments = mrna_fragments + intronic_fragments
            writer.writerow(
                {
                    "sample": sample,
                    "condition": condition,
                    "gene_id": gene_id,
                    "gene_name": gene["gene_name"],
                    "transcript_id": gene["transcript_id"],
                    "chrom": gene["chrom"],
                    "strand": gene["strand"],
                    "transcript_start": gene["transcript_start"],
                    "transcript_end": gene["transcript_end"],
                    "three_prime_end": gene["three_prime_end"],
                    "terminal_exon_start": gene["terminal_exon_start"],
                    "terminal_exon_end": gene["terminal_exon_end"],
                    "transcript_length": gene["transcript_length"],
                    "exonic_length": gene["exonic_length"],
                    "exon_count": gene["exon_count"],
                    "intron_count": gene["intron_count"],
                    "total_gene_fragments": total_gene_fragments,
                    "total_gene_cpm": 0.0
                    if library_fragments == 0
                    else total_gene_fragments * 1_000_000.0 / library_fragments,
                    "mrna_fragments": mrna_fragments,
                    "mrna_cpm": 0.0 if library_fragments == 0 else mrna_fragments * 1_000_000.0 / library_fragments,
                    "mrna_fraction_gene_reads": 0.0 if total_gene_fragments == 0 else mrna_fragments / total_gene_fragments,
                    "mrna_percent_gene_reads": 0.0
                    if total_gene_fragments == 0
                    else mrna_fragments * 100.0 / total_gene_fragments,
                    "intronic_fragments": intronic_fragments,
                    "intronic_cpm": 0.0
                    if library_fragments == 0
                    else intronic_fragments * 1_000_000.0 / library_fragments,
                    "intronic_fraction_gene_reads": 0.0
                    if total_gene_fragments == 0
                    else intronic_fragments / total_gene_fragments,
                    "intronic_percent_gene_reads": 0.0
                    if total_gene_fragments == 0
                    else intronic_fragments * 100.0 / total_gene_fragments,
                    "intronic_to_mrna_ratio": ratio_or_na(intronic_fragments, mrna_fragments),
                }
            )


def write_coverage(path, sample, condition, genes, coverage_deltas, upstream, downstream):
    fieldnames = ["sample", "condition", "gene_id", "gene_name", "offset_nt", "coverage_count"]
    total_coverage_positions = 0
    window_size = upstream + downstream + 1

    with open_text(path, "wt") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        for gene_id in sorted(coverage_deltas, key=lambda value: (genes[value]["gene_name"], value)):
            running = 0
            gene_name = genes[gene_id]["gene_name"]
            deltas = coverage_deltas[gene_id]
            for index in range(window_size):
                running += deltas[index]
                if running == 0:
                    continue
                total_coverage_positions += running
                writer.writerow(
                    {
                        "sample": sample,
                        "condition": condition,
                        "gene_id": gene_id,
                        "gene_name": gene_name,
                        "offset_nt": index - upstream,
                        "coverage_count": running,
                    }
                )

    return total_coverage_positions


def write_summary(path, summary_row):
    with open(path, "w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(summary_row.keys()), delimiter="\t")
        writer.writeheader()
        writer.writerow(summary_row)


def main():
    args = parse_args()
    genes, exon_trees, intron_trees, small_ncrna_trees, small_ncrna_reference_genes = load_reference(
        args.reference,
        args.metaprofile_upstream,
        args.metaprofile_downstream,
    )

    counters = Counter()
    gene_counts = defaultdict(lambda: {"mrna_fragments": 0, "intronic_fragments": 0})
    coverage_deltas = {}

    with pysam.AlignmentFile(args.bam, "rb") as bam:
        for _, alignments in iter_name_collated_groups(args.bam, args.threads):
            counters["query_groups_examined"] += 1

            primary_mapped = [
                alignment
                for alignment in alignments
                if not alignment.is_secondary and not alignment.is_supplementary and not alignment.is_unmapped
            ]
            if not primary_mapped:
                counters["filtered_without_primary_mapped_alignment"] += 1
                continue

            read1_alignments = [alignment for alignment in primary_mapped if alignment.is_read1]
            if len(read1_alignments) != 1:
                counters["filtered_missing_or_multiple_read1"] += 1
                continue
            read1_alignment = read1_alignments[0]

            if read1_alignment.mate_is_unmapped:
                counters["filtered_unmapped_or_mate_unmapped"] += 1
                continue
            if not read1_alignment.is_proper_pair or read1_alignment.reference_id != read1_alignment.next_reference_id:
                counters["filtered_improper_pair"] += 1
                continue
            if read1_alignment.has_tag("NH") and read1_alignment.get_tag("NH") != 1:
                counters["filtered_nonunique"] += 1
                continue
            if not read1_alignment.has_tag("NH") and read1_alignment.mapping_quality <= 0:
                counters["filtered_nonunique"] += 1
                continue

            counters["library_fragments"] += 1

            classification, gene_id = classify_fragment(
                primary_mapped,
                read1_alignment,
                exon_trees,
                intron_trees,
                small_ncrna_trees,
            )
            counters[classification] += 1

            if classification == "mrna":
                gene_counts[gene_id]["mrna_fragments"] += 1
                update_coverage_deltas(
                    coverage_deltas,
                    genes[gene_id],
                    primary_mapped,
                    args.metaprofile_upstream,
                    args.metaprofile_downstream,
                )
            elif classification == "intronic":
                gene_counts[gene_id]["intronic_fragments"] += 1

    total_gene_fragments = counters["mrna"] + counters["intronic"]
    total_coverage_positions = write_coverage(
        args.output_coverage,
        args.sample,
        args.condition,
        genes,
        coverage_deltas,
        args.metaprofile_upstream,
        args.metaprofile_downstream,
    )
    write_gene_counts(
        args.output_gene_counts,
        args.sample,
        args.condition,
        genes,
        gene_counts,
        counters["library_fragments"],
    )

    summary_row = {
        "sample": args.sample,
        "condition": args.condition,
        "reference_genes": len(genes),
        "small_ncrna_reference_genes": small_ncrna_reference_genes,
        "metaprofile_upstream_nt": args.metaprofile_upstream,
        "metaprofile_downstream_nt": args.metaprofile_downstream,
        "query_groups_examined": counters["query_groups_examined"],
        "library_fragments": counters["library_fragments"],
        "filtered_without_primary_mapped_alignment": counters["filtered_without_primary_mapped_alignment"],
        "filtered_missing_or_multiple_read1": counters["filtered_missing_or_multiple_read1"],
        "filtered_unmapped_or_mate_unmapped": counters["filtered_unmapped_or_mate_unmapped"],
        "filtered_improper_pair": counters["filtered_improper_pair"],
        "filtered_nonunique": counters["filtered_nonunique"],
        "fragments_without_mane_gene_overlap": counters["no_gene_overlap"],
        "ambiguous_same_strand_gene_fragments": counters["ambiguous_same_strand_genes"],
        "fragments_excluded_by_small_ncrna_overlap": counters["small_ncrna_overlap"],
        "unclassified_gene_fragments": counters["unclassified"],
        "assigned_gene_fragments": total_gene_fragments,
        "assigned_gene_fragments_cpm": 0.0
        if counters["library_fragments"] == 0
        else total_gene_fragments * 1_000_000.0 / counters["library_fragments"],
        "genes_with_reads": len(gene_counts),
        "mrna_fragments": counters["mrna"],
        "mrna_fragments_cpm": 0.0
        if counters["library_fragments"] == 0
        else counters["mrna"] * 1_000_000.0 / counters["library_fragments"],
        "mrna_genes_with_reads": sum(1 for counts in gene_counts.values() if counts["mrna_fragments"] > 0),
        "mrna_fraction_gene_reads": 0.0 if total_gene_fragments == 0 else counters["mrna"] / total_gene_fragments,
        "mrna_percent_gene_reads": 0.0 if total_gene_fragments == 0 else counters["mrna"] * 100.0 / total_gene_fragments,
        "intronic_fragments": counters["intronic"],
        "intronic_fragments_cpm": 0.0
        if counters["library_fragments"] == 0
        else counters["intronic"] * 1_000_000.0 / counters["library_fragments"],
        "intronic_genes_with_reads": sum(1 for counts in gene_counts.values() if counts["intronic_fragments"] > 0),
        "intronic_fraction_gene_reads": 0.0
        if total_gene_fragments == 0
        else counters["intronic"] / total_gene_fragments,
        "intronic_percent_gene_reads": 0.0
        if total_gene_fragments == 0
        else counters["intronic"] * 100.0 / total_gene_fragments,
        "intronic_to_mrna_ratio": ratio_or_na(counters["intronic"], counters["mrna"]),
        "profiled_coverage_positions": total_coverage_positions,
    }
    write_summary(args.output_summary, summary_row)


if __name__ == "__main__":
    main()
