#!/usr/bin/env python3
"""Build the long-read 3'SS-coverage and branchpoint-deletion source tables.

The implementation preserves the coordinate conventions of the exploratory R
analysis, but avoids expanding every CIGAR deletion to a persistent table and
avoids repeated all-by-all genomic-range joins.  Alignments are streamed twice;
sorted point indexes, interval trees, and coverage difference arrays keep the
runtime and memory use bounded.
"""

from __future__ import annotations

import argparse
import bisect
import csv
import gzip
import re
import sys
from collections import Counter, defaultdict
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, Iterator, Sequence

import pysam


@dataclass(frozen=True, slots=True)
class Intron:
    index: int
    chrom: str
    start: int
    end: int
    strand: str
    transcript_id: str
    gene_id: str
    gene_name: str


@dataclass(frozen=True, slots=True)
class Branchpoint:
    index: int
    chrom: str
    position: int
    strand: str
    score: float
    intron_type: str


@dataclass(frozen=True, slots=True)
class ReadAlignment:
    chrom: str
    start: int
    end: int
    strand: str
    qname: str
    cigar: str


class ContainmentTree:
    """Static interval tree returning intervals that contain a query range."""

    __slots__ = ("center", "left", "right", "crossing")

    def __init__(self, intervals: Sequence[tuple[int, int, int]]):
        if not intervals:
            self.center = None
            self.left = None
            self.right = None
            self.crossing = ()
            return

        self.center = sorted((start + end) // 2 for start, end, _ in intervals)[
            len(intervals) // 2
        ]
        left: list[tuple[int, int, int]] = []
        right: list[tuple[int, int, int]] = []
        crossing: list[tuple[int, int, int]] = []
        for interval in intervals:
            start, end, _ = interval
            if end < self.center:
                left.append(interval)
            elif start > self.center:
                right.append(interval)
            else:
                crossing.append(interval)
        self.crossing = tuple(crossing)
        self.left = ContainmentTree(left) if left else None
        self.right = ContainmentTree(right) if right else None

    def containers(self, query_start: int, query_end: int) -> list[int]:
        if self.center is None:
            return []
        hits: list[int] = []
        for start, end, interval_id in self.crossing:
            if start <= query_start and end >= query_end:
                hits.append(interval_id)
        if query_end < self.center and self.left is not None:
            hits.extend(self.left.containers(query_start, query_end))
        elif query_start > self.center and self.right is not None:
            hits.extend(self.right.containers(query_start, query_end))
        return hits


class PointIndex:
    """Strand-aware, chromosome-aware index for integer genomic positions."""

    def __init__(self, records: Iterable[tuple[str, str, int, int]]):
        grouped: dict[tuple[str, str], list[tuple[int, int]]] = defaultdict(list)
        for chrom, strand, position, record_id in records:
            grouped[(chrom, strand)].append((position, record_id))
        self._positions: dict[tuple[str, str], tuple[list[int], list[int]]] = {}
        for key, values in grouped.items():
            values.sort()
            self._positions[key] = (
                [position for position, _ in values],
                [record_id for _, record_id in values],
            )

    def between(self, chrom: str, strand: str, start: int, end: int) -> list[int]:
        positions_and_ids = self._positions.get((chrom, strand))
        if positions_and_ids is None:
            return []
        positions, record_ids = positions_and_ids
        left = bisect.bisect_left(positions, start)
        right = bisect.bisect_right(positions, end)
        return record_ids[left:right]


class CoverageAccumulator:
    def __init__(self, half_width: int, include_positive_edge: bool = True):
        self.half_width = half_width
        self.maximum_position = half_width if include_positive_edge else half_width - 1
        self._difference = [0] * (
            self.maximum_position + half_width + 2
        )

    def add(self, relative_start: int, relative_end: int, multiplicity: int = 1) -> None:
        relative_start = max(relative_start, -self.half_width)
        relative_end = min(relative_end, self.maximum_position)
        if relative_start > relative_end:
            return
        left = relative_start + self.half_width
        right = relative_end + self.half_width
        self._difference[left] += multiplicity
        self._difference[right + 1] -= multiplicity

    def scores(self) -> list[int]:
        running = 0
        result: list[int] = []
        for delta in self._difference[:-1]:
            running += delta
            result.append(running)
        return result


def open_text(path: Path):
    return gzip.open(path, "rt") if path.suffix == ".gz" else path.open()


def parse_gtf_attributes(text: str) -> dict[str, str]:
    attributes: dict[str, str] = {}
    for field in text.rstrip().split(";"):
        field = field.strip()
        if not field:
            continue
        if "=" in field:
            key, value = field.split("=", 1)
        else:
            key, _, value = field.partition(" ")
        attributes[key] = value.strip().strip('"')
    return attributes


def load_selected_transcripts(path: Path) -> dict[str, str]:
    """Return the first transcript listed for each gene, as in dplyr::slice(1L)."""
    selected: dict[str, str] = {}
    with open_text(path) as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            selected.setdefault(row["gene_id"], row["transcript_id"])
    return selected


def build_introns(gtf_path: Path, transcript_details_path: Path) -> list[Intron]:
    selected_by_gene = load_selected_transcripts(transcript_details_path)
    selected_transcripts = set(selected_by_gene.values())
    exons_by_gene: dict[str, list[tuple[int, int, str, str, str, str]]] = defaultdict(list)

    with open_text(gtf_path) as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) != 9 or fields[2] != "exon":
                continue
            attributes = parse_gtf_attributes(fields[8])
            transcript_id = attributes.get("transcript_id", "")
            if transcript_id not in selected_transcripts:
                continue
            gene_id = attributes.get("gene_id", "")
            if selected_by_gene.get(gene_id) != transcript_id:
                continue
            exons_by_gene[gene_id].append(
                (
                    int(fields[3]),
                    int(fields[4]),
                    fields[0],
                    fields[6],
                    transcript_id,
                    attributes.get("gene_name", ""),
                )
            )

    introns: list[Intron] = []
    for gene_id, exons in exons_by_gene.items():
        if len(exons) <= 1:
            continue
        exons.sort(key=lambda exon: exon[0])
        for left, right in zip(exons, exons[1:]):
            left_start, left_end, chrom, strand, transcript_id, gene_name = left
            right_start, _, right_chrom, right_strand, _, _ = right
            if chrom != right_chrom or strand != right_strand:
                raise ValueError(f"Inconsistent exon coordinates for {transcript_id}")
            introns.append(
                Intron(
                    index=len(introns),
                    chrom=chrom,
                    start=left_end + 1,
                    end=right_start - 1,
                    strand=strand,
                    transcript_id=transcript_id,
                    gene_id=gene_id,
                    gene_name=gene_name,
                )
            )
    return introns


def load_branchpoints(path: Path) -> list[Branchpoint]:
    branchpoints: list[Branchpoint] = []
    with open_text(path) as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            chrom = row.get("chr") or row.get("seqnames")
            branchpoints.append(
                Branchpoint(
                    index=len(branchpoints),
                    chrom=chrom,
                    position=int(row["start"]),
                    strand=row["strand"],
                    score=float(row["score"]),
                    intron_type=row.get("intron_type", ""),
                )
            )
    return branchpoints


def make_containment_indexes(
    introns: Sequence[Intron], selected_ids: Iterable[int] | None = None, stretch: int = 0
) -> dict[tuple[str, str], ContainmentTree]:
    ids = range(len(introns)) if selected_ids is None else selected_ids
    left_extension = stretch // 2
    right_extension = stretch - left_extension
    grouped: dict[tuple[str, str], list[tuple[int, int, int]]] = defaultdict(list)
    for intron_id in ids:
        intron = introns[intron_id]
        grouped[(intron.chrom, intron.strand)].append(
            (intron.start - left_extension, intron.end + right_extension, intron_id)
        )
    return {key: ContainmentTree(values) for key, values in grouped.items()}


def containing_ids(
    indexes: dict[tuple[str, str], ContainmentTree],
    chrom: str,
    strand: str,
    start: int,
    end: int,
) -> list[int]:
    tree = indexes.get((chrom, strand))
    return [] if tree is None else tree.containers(start, end)


def choose_top_branchpoints(
    introns: Sequence[Intron], branchpoints: Sequence[Branchpoint]
) -> list[Branchpoint]:
    intron_index = make_containment_indexes(introns)
    best: dict[int, Branchpoint] = {}
    for branchpoint in branchpoints:
        for intron_id in containing_ids(
            intron_index,
            branchpoint.chrom,
            branchpoint.strand,
            branchpoint.position,
            branchpoint.position,
        ):
            previous = best.get(intron_id)
            if previous is None or branchpoint.score > previous.score:
                best[intron_id] = branchpoint
    return [best[intron_id] for intron_id in sorted(best)]


def introns_with_one_top_branchpoint(
    introns: Sequence[Intron], top_branchpoints: Sequence[Branchpoint]
) -> set[int]:
    # The exploratory code used plyranges::count_overlaps(), whose undirected
    # form ignores strand, for this one filtering step.
    top_index = PointIndex(
        (bp.chrom, "*", bp.position, index)
        for index, bp in enumerate(top_branchpoints)
    )
    selected: set[int] = set()
    for intron in introns:
        if len(
            top_index.between(
                intron.chrom, "*", intron.start, intron.end
            )
        ) == 1:
            selected.add(intron.index)
    return selected


_CIGAR_TOKEN = re.compile(r"(\d+)([MIDNSHP=X])")


def aligned_span_like_r(cigar: str) -> int:
    return sum(
        int(length) for length, operation in _CIGAR_TOKEN.findall(cigar) if operation in "MD"
    )


def deletion_positions_like_r(read_start: int, cigar: str) -> Iterator[int]:
    """Match convert_cigar_to_vector() + which() + start from the R analysis."""
    offset = 0
    for length_text, operation in _CIGAR_TOKEN.findall(cigar):
        length = int(length_text)
        if operation == "M":
            offset += length
        elif operation == "D":
            for within_deletion in range(1, length + 1):
                yield read_start + offset + within_deletion
            offset += length


def iter_bam_reads(path: Path, quality_threshold: int) -> Iterator[ReadAlignment]:
    with pysam.AlignmentFile(path, "rb") as alignment_file:
        for record in alignment_file.fetch(until_eof=True):
            if record.is_unmapped or record.mapping_quality <= quality_threshold:
                continue
            if record.is_secondary or record.is_supplementary:
                continue
            cigar = record.cigarstring
            if cigar is None:
                continue
            start = record.reference_start + 1
            yield ReadAlignment(
                chrom=record.reference_name,
                start=start,
                end=start + aligned_span_like_r(cigar),
                strand="-" if record.is_reverse else "+",
                qname=record.query_name,
                cigar=cigar,
            )


def add_read_coverage(
    accumulator: CoverageAccumulator,
    read: ReadAlignment,
    anchor_ids: Sequence[int],
    anchors: Sequence[tuple[str, int, str]],
    multiplicity: int,
    symmetric_minus: bool = False,
) -> None:
    half_width = accumulator.half_width
    for anchor_id in anchor_ids:
        _, anchor, strand = anchors[anchor_id]
        clipped_start = max(read.start, anchor - half_width)
        clipped_end = min(read.end, anchor + half_width)
        if strand == "+":
            relative_start = clipped_start - anchor
            relative_end = clipped_end - anchor
        elif not symmetric_minus:
            relative_start = anchor - clipped_end
            relative_end = anchor - clipped_start
        else:
            # Preserve the sequential dplyr::mutate() semantics in the original
            # R function.  The second assignment sees the already-flipped
            # relative_start, so reverse-strand intervals become symmetric
            # around the anchor. This preserves the transformation used by the
            # analysis being published.
            original_relative_end = clipped_end - anchor
            relative_start = -original_relative_end
            relative_end = original_relative_end
        accumulator.add(relative_start, relative_end, multiplicity)


def first_pass(
    alignment_path: Path,
    quality_threshold: int,
    intron_candidates: dict[tuple[str, str], ContainmentTree],
    branchpoint_index: PointIndex,
    splice_site_index: PointIndex,
    splice_sites: Sequence[tuple[str, int, str]],
    coverage_half_width: int,
) -> tuple[set[int], CoverageAccumulator, int]:
    covered_introns: set[int] = set()
    # anchor_3p() creates zero-width anchors in the R code.  Stretching those
    # to width 200 gives the asymmetric integer interval -100..99.
    coverage = CoverageAccumulator(coverage_half_width, include_positive_edge=False)
    read_count = 0
    for read in iter_bam_reads(alignment_path, quality_threshold):
        read_count += 1
        splice_site_hits = splice_site_index.between(
            read.chrom, read.strand, read.start, read.end
        )
        add_read_coverage(
            coverage,
            read,
            splice_site_hits,
            splice_sites,
            1,
            symmetric_minus=False,
        )

        if branchpoint_index.between(read.chrom, read.strand, read.start, read.end):
            covered_introns.update(
                containing_ids(
                    intron_candidates,
                    read.chrom,
                    read.strand,
                    read.start,
                    read.end,
                )
            )
    return covered_introns, coverage, read_count


def process_deletion_runs(
    relative_positions: list[int],
    raw_counts: Counter[int],
    consecutive_counts: Counter[int],
) -> None:
    relative_positions.sort()
    raw_counts.update(relative_positions)
    run_start = 0
    for index in range(1, len(relative_positions) + 1):
        if index == len(relative_positions) or (
            relative_positions[index] - relative_positions[index - 1] != 1
        ):
            run_length = index - run_start
            weight = 1.0 / run_length
            for position in relative_positions[run_start:index]:
                consecutive_counts[position] += weight
            run_start = index


def second_pass(
    alignment_path: Path,
    quality_threshold: int,
    selected_introns: dict[tuple[str, str], ContainmentTree],
    branchpoints: Sequence[Branchpoint],
    branchpoint_index: PointIndex,
    top_branchpoints: Sequence[Branchpoint],
    top_branchpoint_index: PointIndex,
    coverage_half_width: int,
) -> tuple[CoverageAccumulator, Counter[int], Counter[int], int, dict[str, int]]:
    coverage = CoverageAccumulator(coverage_half_width)
    raw_counts: Counter[int] = Counter()
    consecutive_counts: Counter[int] = Counter()
    true_read_number = 0
    diagnostics = Counter()
    anchors = [(bp.chrom, bp.position, bp.strand) for bp in top_branchpoints]

    for read in iter_bam_reads(alignment_path, quality_threshold):
        containing_introns = containing_ids(
            selected_introns, read.chrom, read.strand, read.start, read.end
        )
        if containing_introns:
            top_hits = top_branchpoint_index.between(
                read.chrom, read.strand, read.start, read.end
            )
            add_read_coverage(
                coverage,
                read,
                top_hits,
                anchors,
                len(containing_introns),
                symmetric_minus=True,
            )
            diagnostics["high_coverage_read_hits"] += len(containing_introns)

        overlapping_branchpoints = branchpoint_index.between(
            read.chrom, read.strand, read.start, read.end
        )
        if not overlapping_branchpoints:
            continue

        positions_by_branchpoint: dict[int, list[int]] = defaultdict(list)
        retained_deletion = False
        for deletion_position in deletion_positions_like_r(read.start, read.cigar):
            deletion_intron_multiplicity = len(
                containing_ids(
                    selected_introns,
                    read.chrom,
                    read.strand,
                    deletion_position,
                    deletion_position,
                )
            )
            if deletion_intron_multiplicity == 0:
                continue
            retained_deletion = True
            # find_overlaps_directed() repeated a read for every branchpoint it
            # covered, and find_overlaps_within_directed() then repeated each
            # deletion for every containing intron. Preserve those plyranges
            # multiplicities without materializing the expanded table.
            event_multiplicity = (
                len(overlapping_branchpoints) * deletion_intron_multiplicity
            )
            nearby = branchpoint_index.between(
                read.chrom,
                read.strand,
                deletion_position - coverage_half_width,
                deletion_position + coverage_half_width,
            )
            for branchpoint_id in nearby:
                branchpoint = branchpoints[branchpoint_id]
                relative = deletion_position - branchpoint.position
                if read.strand == "-":
                    relative *= -1
                positions_by_branchpoint[branchpoint_id].extend(
                    [relative] * event_multiplicity
                )
                diagnostics["metaprofile_deletion_rows"] += event_multiplicity

        if retained_deletion:
            true_read_number += len(overlapping_branchpoints)
            diagnostics["deletion_reads"] += 1
        for positions in positions_by_branchpoint.values():
            process_deletion_runs(positions, raw_counts, consecutive_counts)

    return coverage, raw_counts, consecutive_counts, true_read_number, dict(diagnostics)


def write_coverage_table(
    path: Path,
    sample_order: Sequence[str],
    sample_coverages: dict[str, CoverageAccumulator],
) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(
            ["protein", "strand_corrected_position", "score", "normalised_coverage"]
        )
        for sample in sample_order:
            scores = sample_coverages[sample].scores()
            maximum = max(scores)
            half_width = sample_coverages[sample].half_width
            for index, score in enumerate(scores):
                writer.writerow(
                    [sample, index - half_width, score, score / maximum]
                )


def write_deletion_table(
    path: Path,
    sample_order: Sequence[str],
    sample_results: dict[
        str, tuple[CoverageAccumulator, Counter[int], Counter[int], int]
    ],
) -> dict[str, tuple[float, float]]:
    path.parent.mkdir(parents=True, exist_ok=True)
    rows_by_sample: dict[str, list[dict[str, float]]] = {}
    with path.open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(
            [
                "protein",
                "strand_corrected_position",
                "count_cons",
                "count",
                "normalised_count",
                "normalised_count_cons",
                "score",
                "normalised_coverage",
                "coverage_normalised",
                "coverage_normalised_cons",
            ]
        )
        for sample in sample_order:
            coverage, raw_counts, consecutive_counts, denominator = sample_results[sample]
            scores = coverage.scores()
            maximum = max(scores)
            half_width = coverage.half_width
            sample_rows: list[dict[str, float]] = []
            for position in sorted(raw_counts):
                score = scores[position + half_width]
                normalised_coverage = score / maximum
                normalised_count = raw_counts[position] / denominator
                normalised_count_cons = consecutive_counts[position] / denominator
                row = {
                    "position": position,
                    "coverage_normalised": normalised_count / normalised_coverage,
                }
                sample_rows.append(row)
                writer.writerow(
                    [
                        sample,
                        position,
                        consecutive_counts[position],
                        raw_counts[position],
                        normalised_count,
                        normalised_count_cons,
                        score,
                        normalised_coverage,
                        normalised_count / normalised_coverage,
                        normalised_count_cons / normalised_coverage,
                    ]
                )
            rows_by_sample[sample] = sample_rows

    auc_values: dict[str, tuple[float, float]] = {}
    for sample, rows in rows_by_sample.items():
        baseline_values = [
            row["coverage_normalised"] for row in rows if row["position"] <= 10
        ]
        baseline_values.sort()
        middle = len(baseline_values) // 2
        if len(baseline_values) % 2:
            baseline = baseline_values[middle]
        else:
            baseline = (baseline_values[middle - 1] + baseline_values[middle]) / 2
        auc = sum(
            row["coverage_normalised"] - baseline
            for row in rows
            if -4 <= row["position"] <= 2
        )
        auc_values[sample] = (baseline, auc)
    return auc_values


def write_auc_table(
    path: Path, sample_order: Sequence[str], values: dict[str, tuple[float, float]]
) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["protein", "baseline", "area_under_peak"])
        for sample in sample_order:
            writer.writerow([sample, *values[sample]])


def parse_sample_paths(values: Sequence[str]) -> tuple[list[str], dict[str, Path]]:
    order: list[str] = []
    paths: dict[str, Path] = {}
    for value in values:
        sample, separator, path_text = value.partition("=")
        if not separator or not sample or not path_text:
            raise ValueError(f"Expected SAMPLE=PATH, got {value!r}")
        order.append(sample)
        paths[sample] = Path(path_text)
    return order, paths


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--alignment", action="append", required=True, metavar="SAMPLE=PATH")
    parser.add_argument("--gtf", required=True, type=Path)
    parser.add_argument("--transcript-details", required=True, type=Path)
    parser.add_argument("--branchpoints", required=True, type=Path)
    parser.add_argument("--coverage-output", required=True, type=Path)
    parser.add_argument("--deletion-output", required=True, type=Path)
    parser.add_argument("--auc-output", required=True, type=Path)
    parser.add_argument("--diagnostics-output", required=True, type=Path)
    parser.add_argument("--mapq-threshold", type=int, default=20)
    parser.add_argument("--three-prime-half-width", type=int, default=100)
    parser.add_argument("--branchpoint-half-width", type=int, default=60)
    parser.add_argument("--intron-stretch", type=int, default=6)
    args = parser.parse_args()

    sample_order, alignment_paths = parse_sample_paths(args.alignment)
    if len(sample_order) != 2:
        parser.error("This experiment and shared-intron analysis require exactly two samples")
    introns = build_introns(args.gtf, args.transcript_details)
    branchpoints = load_branchpoints(args.branchpoints)
    top_branchpoints = choose_top_branchpoints(introns, branchpoints)
    one_bp_introns = introns_with_one_top_branchpoint(introns, top_branchpoints)

    branchpoint_index = PointIndex(
        (bp.chrom, bp.strand, bp.position, bp.index) for bp in branchpoints
    )
    splice_sites = [
        (
            intron.chrom,
            intron.end if intron.strand == "+" else intron.start,
            intron.strand,
        )
        for intron in introns
    ]
    splice_site_index = PointIndex(
        (chrom, strand, position, index)
        for index, (chrom, position, strand) in enumerate(splice_sites)
    )
    candidate_index = make_containment_indexes(
        introns, one_bp_introns, stretch=args.intron_stretch
    )

    covered_by_sample: dict[str, set[int]] = {}
    three_prime_coverage: dict[str, CoverageAccumulator] = {}
    read_counts: dict[str, int] = {}
    for sample in sample_order:
        print(f"[{sample}] first alignment pass", file=sys.stderr, flush=True)
        covered, coverage, read_count = first_pass(
            alignment_paths[sample],
            args.mapq_threshold,
            candidate_index,
            branchpoint_index,
            splice_site_index,
            splice_sites,
            args.three_prime_half_width,
        )
        covered_by_sample[sample] = covered
        three_prime_coverage[sample] = coverage
        read_counts[sample] = read_count

    shared_analysis_introns = set().union(*covered_by_sample.values())
    selected_intron_index = make_containment_indexes(
        introns, shared_analysis_introns, stretch=args.intron_stretch
    )
    top_branchpoint_index = PointIndex(
        (bp.chrom, bp.strand, bp.position, index)
        for index, bp in enumerate(top_branchpoints)
    )

    deletion_results: dict[
        str, tuple[CoverageAccumulator, Counter[int], Counter[int], int]
    ] = {}
    diagnostics_by_sample: dict[str, dict[str, int]] = {}
    for sample in sample_order:
        print(f"[{sample}] second alignment pass", file=sys.stderr, flush=True)
        coverage, raw, consecutive, denominator, diagnostics = second_pass(
            alignment_paths[sample],
            args.mapq_threshold,
            selected_intron_index,
            branchpoints,
            branchpoint_index,
            top_branchpoints,
            top_branchpoint_index,
            args.branchpoint_half_width,
        )
        deletion_results[sample] = (coverage, raw, consecutive, denominator)
        diagnostics.update(
            {
                "filtered_primary_alignments": read_counts[sample],
                "covered_candidate_introns": len(covered_by_sample[sample]),
                "deletion_normalisation_read_overlaps": denominator,
            }
        )
        diagnostics_by_sample[sample] = diagnostics

    write_coverage_table(args.coverage_output, sample_order, three_prime_coverage)
    auc_values = write_deletion_table(args.deletion_output, sample_order, deletion_results)
    write_auc_table(args.auc_output, sample_order, auc_values)

    args.diagnostics_output.parent.mkdir(parents=True, exist_ok=True)
    with args.diagnostics_output.open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["scope", "metric", "value"])
        writer.writerow(["reference", "introns", len(introns)])
        writer.writerow(["reference", "branchpoints", len(branchpoints)])
        writer.writerow(["reference", "top_branchpoints", len(top_branchpoints)])
        writer.writerow(["reference", "single_top_branchpoint_introns", len(one_bp_introns)])
        writer.writerow(["analysis", "union_covered_introns", len(shared_analysis_introns)])
        for sample in sample_order:
            for metric, value in sorted(diagnostics_by_sample[sample].items()):
                writer.writerow([sample, metric, value])


if __name__ == "__main__":
    main()
