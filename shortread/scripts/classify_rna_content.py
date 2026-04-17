#!/usr/bin/env python3

import argparse
import collections
import gzip
import pickle
import shutil
import signal
import subprocess
import sys

import pysam
from intervaltree import IntervalTree


OUTPUT_ORDER = [
    "rRNA",
    "snRNA",
    "mitochondrial",
    "snoRNA",
    "scaRNA",
    "miRNA",
    "tRNA",
    "ribozyme",
    "vaultRNA",
    "misc_RNA",
    "lncRNA",
    "protein_coding_intron",
    "protein_coding_exon",
    "pseudogene",
    "other_annotated_genic",
    "intergenic",
    "ambiguous",
    "genome_unmapped",
]

PRIORITY_GROUPS = [
    ["mitochondrial"],
    ["snoRNA", "scaRNA", "miRNA", "tRNA", "ribozyme", "vaultRNA", "misc_RNA", "snRNA"],
    ["lncRNA"],
    ["protein_coding_intron"],
    ["protein_coding_exon"],
    ["pseudogene"],
]

MITO_CONTIGS = {"chrM", "MT", "M"}

SMALL_NCRNA_CATEGORIES = {
    "snorna": "snoRNA",
    "scarna": "scaRNA",
    "mirna": "miRNA",
    "trna": "tRNA",
    "ribozyme": "ribozyme",
    "vaultrna": "vaultRNA",
    "vault_rna": "vaultRNA",
    "misc_rna": "misc_RNA",
    "snrna": "snRNA",
}

LNC_RNA_TYPES = {
    "lncrna",
    "lincrna",
    "antisense",
    "sense_intronic",
    "sense_overlapping",
    "processed_transcript",
    "3prime_overlapping_ncrna",
    "3prime_overlapping_ncrna",
    "macro_lncrna",
    "non_coding",
    "bidirectional_promoter_lncrna",
}

ANNOTATION_CATEGORIES = (
    "snoRNA",
    "scaRNA",
    "miRNA",
    "tRNA",
    "ribozyme",
    "vaultRNA",
    "misc_RNA",
    "snRNA",
    "lncRNA",
    "protein_coding_intron",
    "protein_coding_exon",
    "pseudogene",
)


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("--umi-ready-r1", required=True)
    parser.add_argument("--star-input-r1", required=True)
    parser.add_argument("--rrna-bam", required=True)
    parser.add_argument("--snrna-bam", required=True)
    parser.add_argument("--star-bam", required=True)
    parser.add_argument("--gtf")
    parser.add_argument("--reference-cache")
    parser.add_argument("--threads", type=int, default=1)
    parser.add_argument("--output", required=True)
    args = parser.parse_args()
    if args.gtf is None and args.reference_cache is None:
        parser.error("one of --gtf or --reference-cache is required")
    return args


def open_text(path):
    if path.endswith(".gz"):
        return gzip.open(path, "rt")
    return open(path)


def count_fastq_reads(path, threads=1):
    if path.endswith(".gz"):
        pigz = shutil.which("pigz")
        if pigz is not None:
            command = [pigz, "-dc", "-p", str(max(1, threads)), path]
        else:
            command = ["gzip", "-dc", path]
        process = subprocess.Popen(command, stdout=subprocess.PIPE, text=False)
        try:
            if process.stdout is None:
                raise RuntimeError("Failed to open FASTQ decompression stream")
            counter = subprocess.run(
                ["wc", "-l"],
                stdin=process.stdout,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                text=True,
                check=True,
            )
        finally:
            if process.stdout is not None:
                process.stdout.close()
        return_code = process.wait()
        if return_code != 0:
            raise subprocess.CalledProcessError(return_code, command)
        return int(counter.stdout.strip().split()[0]) // 4

    count = 0
    with open_text(path) as handle:
        while True:
            header = handle.readline()
            if not header:
                break
            handle.readline()
            handle.readline()
            handle.readline()
            count += 1
    return count


def count_unique_query_names(path):
    names = set()
    with pysam.AlignmentFile(path, "rb") as bam:
        for aln in bam.fetch(until_eof=True):
            if aln.is_secondary or aln.is_supplementary:
                continue
            names.add(aln.query_name)
    return len(names)


def normalize_gene_type(raw_value):
    return raw_value.lower().replace("-", "_").replace(" ", "_")


def parse_attributes(raw_attributes):
    attributes = {}
    for part in raw_attributes.strip().split(";"):
        part = part.strip()
        if not part:
            continue
        key, _, value = part.partition(" ")
        attributes[key] = value.strip().strip('"')
    return attributes


def gene_type_to_category(normalized_gene_type):
    if normalized_gene_type == "protein_coding":
        return "protein_coding"
    if normalized_gene_type in SMALL_NCRNA_CATEGORIES:
        return SMALL_NCRNA_CATEGORIES[normalized_gene_type]
    if normalized_gene_type in LNC_RNA_TYPES:
        return "lncRNA"
    if "pseudogene" in normalized_gene_type:
        return "pseudogene"
    return None


def merge_spans(spans):
    merged = []
    for start, end in sorted(spans):
        if not merged or start > merged[-1][1]:
            merged.append([start, end])
        else:
            merged[-1][1] = max(merged[-1][1], end)
    return [(start, end) for start, end in merged]


def parse_gtf(gtf_path):
    category_trees = collections.defaultdict(lambda: collections.defaultdict(IntervalTree))
    annotated_gene_trees = collections.defaultdict(IntervalTree)
    protein_coding_genes = {}
    protein_coding_exons = collections.defaultdict(list)

    with open_text(gtf_path) as handle:
        for line in handle:
            if not line or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) != 9:
                continue
            chrom, _, feature, start, end, _, _, _, attributes_raw = fields
            attributes = parse_attributes(attributes_raw)
            gene_id = attributes.get("gene_id")
            gene_type = attributes.get("gene_type") or attributes.get("gene_biotype")
            if gene_id is None or gene_type is None:
                continue

            start0 = int(start) - 1
            end0 = int(end)

            if feature == "gene":
                annotated_gene_trees[chrom].addi(start0, end0)
                category = gene_type_to_category(normalize_gene_type(gene_type))
                if category == "protein_coding":
                    protein_coding_genes[gene_id] = (chrom, start0, end0)
                elif category is not None:
                    category_trees[category][chrom].addi(start0, end0)
            elif feature == "exon" and normalize_gene_type(gene_type) == "protein_coding":
                protein_coding_exons[gene_id].append((chrom, start0, end0))

    for gene_id, (chrom, gene_start, gene_end) in protein_coding_genes.items():
        exon_spans = merge_spans(
            [(start, end) for exon_chrom, start, end in protein_coding_exons.get(gene_id, []) if exon_chrom == chrom]
        )
        for exon_start, exon_end in exon_spans:
            category_trees["protein_coding_exon"][chrom].addi(exon_start, exon_end)

        cursor = gene_start
        for exon_start, exon_end in exon_spans:
            if cursor < exon_start:
                category_trees["protein_coding_intron"][chrom].addi(cursor, exon_start)
            cursor = max(cursor, exon_end)
        if cursor < gene_end:
            category_trees["protein_coding_intron"][chrom].addi(cursor, gene_end)

    for chrom_trees in category_trees.values():
        for tree in chrom_trees.values():
            tree.merge_overlaps(strict=False)
    for tree in annotated_gene_trees.values():
        tree.merge_overlaps(strict=False)

    frozen_category_trees = {
        category: dict(category_trees[category])
        for category in ANNOTATION_CATEGORIES
    }
    frozen_annotated_gene_trees = dict(annotated_gene_trees)

    return frozen_category_trees, frozen_annotated_gene_trees


def load_reference_cache(path):
    with open(path, "rb") as handle:
        category_trees, annotated_gene_trees = pickle.load(handle)

    missing_categories = [
        category for category in ANNOTATION_CATEGORIES
        if category not in category_trees
    ]
    if missing_categories:
        missing_text = ", ".join(sorted(missing_categories))
        raise ValueError(
            f"RNA-content reference cache at {path} is stale or incompatible; "
            f"missing categories: {missing_text}. Rebuild the cache."
        )

    return category_trees, annotated_gene_trees


def tree_has_overlap(tree_by_chrom, chrom, start, end):
    tree = tree_by_chrom.get(chrom)
    if tree is None:
        return False
    return bool(tree.overlap(start, end))


def classify_fragment(alignments, category_trees, annotated_gene_trees):
    overlaps = set()
    annotated = False

    for aln in alignments:
        chrom = aln.reference_name
        blocks = aln.get_blocks()
        if chrom in MITO_CONTIGS:
            overlaps.add("mitochondrial")
        for start, end in blocks:
            if not annotated and tree_has_overlap(annotated_gene_trees, chrom, start, end):
                annotated = True
            for category in ANNOTATION_CATEGORIES:
                if tree_has_overlap(category_trees[category], chrom, start, end):
                    overlaps.add(category)

    for group in PRIORITY_GROUPS:
        hits = [category for category in group if category in overlaps]
        if len(hits) == 1:
            return hits[0]
        if len(hits) > 1:
            return "ambiguous"

    if annotated:
        return "other_annotated_genic"
    return "intergenic"


def iter_name_collated_groups(star_bam_path, threads):
    # Fast collate keeps only primary alignments, which is all this classifier uses.
    command = ["samtools", "collate", "-f", "-u", "-O", "-@", str(threads), star_bam_path]
    process = subprocess.Popen(command, stdout=subprocess.PIPE)
    try:
        with pysam.AlignmentFile(process.stdout, "rb") as bam:
            current_name = None
            bucket = []
            for aln in bam.fetch(until_eof=True):
                if current_name is None:
                    current_name = aln.query_name
                if aln.query_name != current_name:
                    yield current_name, bucket
                    current_name = aln.query_name
                    bucket = [aln]
                else:
                    bucket.append(aln)
            if bucket:
                yield current_name, bucket
    finally:
        if process.stdout is not None:
            process.stdout.close()
        return_code = process.wait()
        if return_code not in (0, -signal.SIGPIPE):
            raise subprocess.CalledProcessError(return_code, command)


def classify_star_bam(star_bam_path, category_trees, annotated_gene_trees, threads):
    counts = collections.Counter()
    mapped_pairs = 0

    for _, alignments in iter_name_collated_groups(star_bam_path, threads):
        primary_mapped = [
            aln
            for aln in alignments
            if not aln.is_unmapped
        ]
        if not primary_mapped:
            continue
        category = classify_fragment(primary_mapped, category_trees, annotated_gene_trees)
        counts[category] += 1
        mapped_pairs += 1

    return counts, mapped_pairs


def main():
    args = parse_args()

    if args.reference_cache is not None:
        category_trees, annotated_gene_trees = load_reference_cache(args.reference_cache)
    else:
        category_trees, annotated_gene_trees = parse_gtf(args.gtf)

    total_input_pairs = count_fastq_reads(args.umi_ready_r1, args.threads)
    star_input_pairs = count_fastq_reads(args.star_input_r1, args.threads)
    rrna_pairs = count_unique_query_names(args.rrna_bam)
    snrna_pairs = count_unique_query_names(args.snrna_bam)

    genomic_counts, genomic_mapped_pairs = classify_star_bam(
        args.star_bam,
        category_trees,
        annotated_gene_trees,
        args.threads,
    )

    genome_unmapped_pairs = star_input_pairs - genomic_mapped_pairs
    if genome_unmapped_pairs < 0:
        raise ValueError(
            f"Computed a negative genome_unmapped count ({genome_unmapped_pairs})."
        )

    counts = collections.Counter()
    counts["rRNA"] = rrna_pairs
    counts["snRNA"] = snrna_pairs + genomic_counts.get("snRNA", 0)
    for category, count in genomic_counts.items():
        if category == "snRNA":
            continue
        counts[category] = count
    counts["genome_unmapped"] = genome_unmapped_pairs

    classified_pairs = total_input_pairs - counts["genome_unmapped"]

    with open(args.output, "w") as handle:
        handle.write("category\tcount\tfraction_of_total_input_pairs\tfraction_of_classified_pairs\n")
        for category in OUTPUT_ORDER:
            count = counts.get(category, 0)
            fraction_total = count / total_input_pairs if total_input_pairs else 0.0
            if category == "genome_unmapped":
                fraction_classified = "NA"
            else:
                fraction_classified = (
                    f"{count / classified_pairs:.6f}" if classified_pairs else "NA"
                )
            handle.write(
                f"{category}\t{count}\t{fraction_total:.6f}\t{fraction_classified}\n"
            )

    sys.stderr.write(f"total_input_pairs\t{total_input_pairs}\n")
    sys.stderr.write(f"rrna_pairs\t{rrna_pairs}\n")
    sys.stderr.write(f"snrna_pairs\t{snrna_pairs}\n")
    sys.stderr.write(f"star_input_pairs\t{star_input_pairs}\n")
    sys.stderr.write(f"genomic_mapped_pairs\t{genomic_mapped_pairs}\n")
    sys.stderr.write(f"genome_unmapped_pairs\t{genome_unmapped_pairs}\n")
    sys.stderr.write(
        f"classified_pair_sum\t{sum(counts.get(category, 0) for category in OUTPUT_ORDER)}\n"
    )


if __name__ == "__main__":
    main()
