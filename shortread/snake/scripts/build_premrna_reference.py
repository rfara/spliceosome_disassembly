#!/usr/bin/env python3

import argparse
import gzip
import pickle
from collections import Counter, defaultdict


SMALL_NCRNA_GENE_TYPES = {
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


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("--gtf", required=True)
    parser.add_argument("--output", required=True)
    return parser.parse_args()


def open_text(path):
    if path.endswith(".gz"):
        return gzip.open(path, "rt")
    return open(path)


def parse_attributes(raw_attributes):
    attributes = {}
    for part in raw_attributes.strip().split(";"):
        part = part.strip()
        if not part:
            continue
        key, _, value = part.partition(" ")
        attributes[key] = value.strip().strip('"')
    return attributes


def normalize_gene_type(raw_value):
    return raw_value.lower().replace("-", "_").replace(" ", "_")


def merge_spans(spans):
    merged = []
    for start, end in sorted(spans):
        if not merged or start > merged[-1][1] + 1:
            merged.append([start, end])
        else:
            merged[-1][1] = max(merged[-1][1], end)
    return [(start, end) for start, end in merged]


def parse_gtf(gtf_path):
    transcript_rows = {}
    exons_by_transcript = defaultdict(list)
    small_ncrna_genes = []
    small_ncrna_counts = Counter()

    with open_text(gtf_path) as handle:
        for line in handle:
            if not line or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) != 9:
                continue

            chrom, _, feature, start, end, _, strand, _, raw_attributes = fields
            attributes = parse_attributes(raw_attributes)
            gene_id = attributes.get("gene_id")
            gene_type = attributes.get("gene_type") or attributes.get("gene_biotype")
            if gene_id is None or gene_type is None:
                continue

            start_1based = int(start)
            end_1based = int(end)
            normalized_gene_type = normalize_gene_type(gene_type)

            if feature == "gene" and normalized_gene_type in SMALL_NCRNA_GENE_TYPES:
                category = SMALL_NCRNA_GENE_TYPES[normalized_gene_type]
                small_ncrna_genes.append(
                    {
                        "gene_id": gene_id,
                        "gene_name": attributes.get("gene_name", gene_id),
                        "chrom": chrom,
                        "strand": strand,
                        "start": start_1based,
                        "end": end_1based,
                        "category": category,
                    }
                )
                small_ncrna_counts[category] += 1

            if feature not in {"transcript", "exon"}:
                continue
            if normalized_gene_type != "protein_coding":
                continue
            if 'tag "MANE_Select"' not in raw_attributes:
                continue

            transcript_id = attributes.get("transcript_id")
            if transcript_id is None:
                continue

            if feature == "transcript":
                transcript_rows[transcript_id] = {
                    "gene_id": gene_id,
                    "gene_name": attributes.get("gene_name", gene_id),
                    "chrom": chrom,
                    "strand": strand,
                    "transcript_start": start_1based,
                    "transcript_end": end_1based,
                }
            else:
                exons_by_transcript[transcript_id].append((start_1based, end_1based))

    genes = {}
    for transcript_id, transcript_meta in transcript_rows.items():
        merged_exons = merge_spans(exons_by_transcript.get(transcript_id, []))
        if not merged_exons:
            continue

        introns = []
        for exon_index in range(len(merged_exons) - 1):
            intron_start = merged_exons[exon_index][1] + 1
            intron_end = merged_exons[exon_index + 1][0] - 1
            if intron_start <= intron_end:
                introns.append((intron_start, intron_end))

        gene_id = transcript_meta["gene_id"]
        if gene_id in genes:
            raise ValueError(
                f"Expected at most one MANE Select transcript per gene, found duplicate gene_id {gene_id}"
            )

        strand = transcript_meta["strand"]
        terminal_exon_start, terminal_exon_end = merged_exons[-1] if strand == "+" else merged_exons[0]
        three_prime_end = transcript_meta["transcript_end"] if strand == "+" else transcript_meta["transcript_start"]

        genes[gene_id] = {
            "gene_id": gene_id,
            "gene_name": transcript_meta["gene_name"],
            "transcript_id": transcript_id,
            "chrom": transcript_meta["chrom"],
            "strand": strand,
            "transcript_start": transcript_meta["transcript_start"],
            "transcript_end": transcript_meta["transcript_end"],
            "three_prime_end": three_prime_end,
            "terminal_exon_start": terminal_exon_start,
            "terminal_exon_end": terminal_exon_end,
            "transcript_length": transcript_meta["transcript_end"] - transcript_meta["transcript_start"] + 1,
            "exonic_length": sum(end - start + 1 for start, end in merged_exons),
            "exon_count": len(merged_exons),
            "intron_count": len(introns),
            "exon_spans": merged_exons,
            "intron_spans": introns,
        }

    return {
        "genes": genes,
        "small_ncrna_genes": small_ncrna_genes,
        "small_ncrna_counts": dict(sorted(small_ncrna_counts.items())),
    }


def main():
    args = parse_args()
    reference = parse_gtf(args.gtf)
    with open(args.output, "wb") as handle:
        pickle.dump(reference, handle, protocol=pickle.HIGHEST_PROTOCOL)

    print(f"MANE Select protein-coding genes loaded: {len(reference['genes'])}")
    print(f"Small ncRNA genes loaded: {len(reference['small_ncrna_genes'])}")
    for category, count in reference["small_ncrna_counts"].items():
        print(f"  {category}: {count}")


if __name__ == "__main__":
    main()
