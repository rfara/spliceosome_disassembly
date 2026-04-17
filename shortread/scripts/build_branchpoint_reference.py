#!/usr/bin/env python3

import argparse
import csv
import gzip
from collections import defaultdict
from dataclasses import dataclass, field


@dataclass
class BranchpointCandidate:
    position: int
    score: float
    intron_type: str


@dataclass
class IntronRecord:
    transcript_id: str
    intron_number: int
    gene_id: str
    gene_name: str
    chrom: str
    strand: str
    start: int
    end: int
    candidates: list[BranchpointCandidate] = field(default_factory=list)

    @property
    def intron_id(self):
        return f"{self.transcript_id}:intron_{self.intron_number}"

    @property
    def length(self):
        return self.end - self.start + 1

    @property
    def three_prime_ss(self):
        return self.end if self.strand == "+" else self.start

    def distance_to_three_prime(self, position):
        if self.strand == "+":
            return self.end - position
        return position - self.start


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("--gtf", required=True)
    parser.add_argument("--branchpoints", required=True)
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


def parse_mane_select_introns(gtf_path):
    transcript_rows = {}
    exons_by_transcript = defaultdict(list)

    with open_text(gtf_path) as handle:
        for line in handle:
            if not line or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) != 9:
                continue
            chrom, _, feature, start, end, _, strand, _, raw_attributes = fields
            if feature not in {"transcript", "exon"}:
                continue
            if 'gene_type "protein_coding"' not in raw_attributes:
                continue
            if 'tag "MANE_Select"' not in raw_attributes:
                continue

            attributes = parse_attributes(raw_attributes)
            transcript_id = attributes.get("transcript_id")
            gene_id = attributes.get("gene_id")
            if transcript_id is None or gene_id is None:
                continue

            if feature == "transcript":
                transcript_rows[transcript_id] = {
                    "chrom": chrom,
                    "strand": strand,
                    "gene_id": gene_id,
                    "gene_name": attributes.get("gene_name", gene_id),
                }
            else:
                exons_by_transcript[transcript_id].append((int(start), int(end)))

    introns_by_key = defaultdict(list)
    introns = []
    for transcript_id, transcript_meta in transcript_rows.items():
        exons = sorted(exons_by_transcript.get(transcript_id, []))
        for idx in range(len(exons) - 1):
            intron_start = exons[idx][1] + 1
            intron_end = exons[idx + 1][0] - 1
            if intron_start > intron_end:
                continue
            intron = IntronRecord(
                transcript_id=transcript_id,
                intron_number=idx + 1,
                gene_id=transcript_meta["gene_id"],
                gene_name=transcript_meta["gene_name"],
                chrom=transcript_meta["chrom"],
                strand=transcript_meta["strand"],
                start=intron_start,
                end=intron_end,
            )
            introns_by_key[(intron.chrom, intron.strand)].append(intron)
            introns.append(intron)

    for key in introns_by_key:
        introns_by_key[key].sort(key=lambda intron: (intron.start, intron.end))

    return introns, introns_by_key


def load_branchpoints(path):
    branchpoints_by_key = defaultdict(list)
    total_rows = 0
    with open_text(path) as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            position = int(row["start"])
            score = float(row["score"])
            branchpoints_by_key[(row["chr"], row["strand"])].append(
                (position, score, row["intron_type"])
            )
            total_rows += 1

    for key in branchpoints_by_key:
        branchpoints_by_key[key].sort()

    return branchpoints_by_key, total_rows


def assign_branchpoints(introns_by_key, branchpoints_by_key):
    assigned_rows = 0
    unassigned_rows = 0
    ambiguous_rows = 0

    for key, branchpoints in branchpoints_by_key.items():
        introns = introns_by_key.get(key, [])
        intron_idx = 0
        active_introns = []

        for position, score, intron_type in branchpoints:
            while intron_idx < len(introns) and introns[intron_idx].start <= position:
                active_introns.append(introns[intron_idx])
                intron_idx += 1
            active_introns = [intron for intron in active_introns if intron.end >= position]

            if not active_introns:
                unassigned_rows += 1
                continue

            assigned_rows += 1
            if len(active_introns) > 1:
                ambiguous_rows += 1

            for intron in active_introns:
                intron.candidates.append(
                    BranchpointCandidate(
                        position=position,
                        score=score,
                        intron_type=intron_type,
                    )
                )

    return assigned_rows, unassigned_rows, ambiguous_rows


def choose_best_candidate(intron):
    def ranking_key(candidate):
        distance = intron.distance_to_three_prime(candidate.position)
        strand_tiebreak = candidate.position if intron.strand == "+" else -candidate.position
        return (candidate.score, -distance, strand_tiebreak)

    return max(intron.candidates, key=ranking_key)


def write_reference(introns, output_path):
    selected_introns = [intron for intron in introns if intron.candidates]
    selected_introns.sort(key=lambda intron: (intron.chrom, intron.start, intron.end, intron.transcript_id))

    fieldnames = [
        "intron_id",
        "transcript_id",
        "intron_number",
        "gene_id",
        "gene_name",
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
    ]

    with open(output_path, "w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        for intron in selected_introns:
            best_candidate = choose_best_candidate(intron)
            writer.writerow(
                {
                    "intron_id": intron.intron_id,
                    "transcript_id": intron.transcript_id,
                    "intron_number": intron.intron_number,
                    "gene_id": intron.gene_id,
                    "gene_name": intron.gene_name,
                    "chrom": intron.chrom,
                    "strand": intron.strand,
                    "intron_start": intron.start,
                    "intron_end": intron.end,
                    "intron_length": intron.length,
                    "three_prime_ss": intron.three_prime_ss,
                    "branchpoint_position": best_candidate.position,
                    "branchpoint_score": best_candidate.score,
                    "branchpoint_to_3ss_nt": intron.distance_to_three_prime(best_candidate.position),
                    "branchpoint_type": best_candidate.intron_type,
                    "branchpoint_candidates": len(intron.candidates),
                }
            )

    return len(selected_introns), sum(1 for intron in selected_introns if len(intron.candidates) > 1)


def main():
    args = parse_args()
    introns, introns_by_key = parse_mane_select_introns(args.gtf)
    branchpoints_by_key, total_branchpoint_rows = load_branchpoints(args.branchpoints)
    assigned_rows, unassigned_rows, ambiguous_rows = assign_branchpoints(introns_by_key, branchpoints_by_key)
    selected_introns, multi_candidate_introns = write_reference(introns, args.output)

    print(f"MANE-select protein-coding introns: {len(introns)}")
    print(f"Branchpoint rows loaded: {total_branchpoint_rows}")
    print(f"Branchpoint rows assigned to at least one intron: {assigned_rows}")
    print(f"Branchpoint rows outside all MANE-select introns: {unassigned_rows}")
    print(f"Branchpoint rows overlapping multiple introns: {ambiguous_rows}")
    print(f"Introns with at least one branchpoint: {selected_introns}")
    print(f"Introns with multiple branchpoint candidates: {multi_candidate_introns}")


if __name__ == "__main__":
    main()
